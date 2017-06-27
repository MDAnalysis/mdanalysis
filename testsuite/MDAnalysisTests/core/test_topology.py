"""Tests for MDAnalysis.core.topology objects.

Should convert between indices (*ix)
Should work with both a single or an array of indices
"""
from __future__ import absolute_import
from six.moves import zip
import itertools

from numpy.testing import (
    assert_,
    assert_equal,
    assert_array_equal,
    assert_raises,
)
import numpy as np

from MDAnalysisTests import make_Universe

from MDAnalysis.core.topology import (
    Topology,
    TransTable,
    make_downshift_arrays,
)
from MDAnalysis.core import topologyattrs as ta
from MDAnalysis.core import groups
from MDAnalysis import NoDataError
import MDAnalysis


class TestTransTable(object):
    def setUp(self):
        # Reference data
        self.Ridx = np.array([0, 0, 2, 2, 1, 1, 3, 3, 1, 2])
        self.Sidx = np.array([0, 1, 1, 0])
        self.tt = TransTable(10, 4, 2, self.Ridx, self.Sidx)

    def tearDown(self):
        del self.tt

    def test_a2r(self):
        for aix, rix in zip(
                [np.array([0, 1, 2]),
                 np.array([9, 6, 2]),
                 np.array([3, 3, 3])],
                [np.array([0, 0, 2]),
                 np.array([2, 3, 2]),
                 np.array([2, 2, 2])]
        ):
            assert_array_equal(self.tt.atoms2residues(aix), rix)

    def test_r2a_1d(self):
        for rix, aix in zip(
                [[0, 1], [1, 1], [3, 1]],
                [[0, 1, 4, 5, 8], [4, 5, 8, 4, 5, 8], [6, 7, 4, 5, 8]]
        ):
            assert_array_equal(self.tt.residues2atoms_1d(rix), aix)

    def test_r2a_2d(self):
        for rix, aix in zip(
                [[0, 1],
                 [1, 1],
                 [3, 1]],
                [[[0, 1], [4, 5, 8]],
                 [[4, 5, 8], [4, 5, 8]],
                 [[6, 7], [4, 5, 8]]]
        ):
            answer = self.tt.residues2atoms_2d(rix)
            for a1, a2 in zip(answer, aix):
                assert_array_equal(a1, a2)

    def test_r2s(self):
        for rix, sidx in zip(
                [np.array([0, 1]),
                 np.array([2, 1, 0]),
                 np.array([1, 1, 1])],
                [np.array([0, 1]),
                 np.array([1, 1, 0]),
                 np.array([1, 1, 1])]
        ):
            assert_array_equal(self.tt.residues2segments(rix), sidx)

    def test_s2r_1d(self):
        for sidx, rix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[0, 3, 1, 2],
                 [1, 2, 0, 3],
                 [1, 2, 1, 2]]
                ):
            assert_array_equal(self.tt.segments2residues_1d(sidx), rix)

    def test_s2r_2d(self):
        for sidx, rix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[[0, 3], [1, 2]],
                 [[1, 2], [0, 3]],
                 [[1, 2], [1, 2]]]
        ):
            answer = self.tt.segments2residues_2d(sidx)
            for a1, a2 in zip(answer, rix):
                assert_array_equal(a1, a2)

    def test_s2a_1d(self):
        for sidx, aix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[0, 1, 6, 7, 4, 5, 8, 2, 3, 9],
                 [4, 5, 8, 2, 3, 9, 0, 1, 6, 7],
                 [4, 5, 8, 2, 3, 9, 4, 5, 8, 2, 3, 9]],
        ):
            assert_array_equal(self.tt.segments2atoms_1d(sidx), aix)

    def test_s2a_2d(self):
        for sidx, aix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[[0, 1, 6, 7], [4, 5, 8, 2, 3, 9]],
                 [[4, 5, 8, 2, 3, 9], [0, 1, 6, 7]],
                 [[4, 5, 8, 2, 3, 9], [4, 5, 8, 2, 3, 9]]],
        ):
            answer = self.tt.segments2atoms_2d(sidx)
            for a1, a2 in zip(answer, aix):
                assert_array_equal(a1, a2)

    def test_move_atom_simple(self):
        tt = self.tt
        assert_equal(tt.atoms2residues(1), 0)
        assert_equal(len(tt.residues2atoms_1d(0)), 2)
        assert_equal(len(tt.residues2atoms_1d(3)), 2)

        # move 2nd atom to 4th residue (atom1 -> res3)
        tt.move_atom(1, 3)

        assert_equal(tt.atoms2residues(1), 3)  # identity of parent changed
        assert_equal(len(tt.residues2atoms_1d(0)), 1)  # 1 fewer here
        assert_equal(len(tt.residues2atoms_1d(3)), 3)  # 1 more here

    def test_move_residue_simple(self):
        tt = self.tt
        assert_equal(tt.residues2segments(1), 1)
        assert_equal(len(tt.segments2residues_1d(0)), 2)
        assert_equal(len(tt.segments2residues_1d(1)), 2)

        # move 2nd residue to 1st segment (res1 -> seg0)
        tt.move_residue(1, 0)

        assert_equal(tt.residues2segments(1), 0)
        assert_equal(len(tt.segments2residues_1d(0)), 3)
        assert_equal(len(tt.segments2residues_1d(1)), 1)


class TestLevelMoves(object):
    """Tests for moving atoms/residues between residues/segments

    
    Atoms can move between residues by setting .residue with a Residue
    Residues can move between segments by setting .segment with a Segment

    Moves are performed by setting either [res/seg]indices or [res/seg]ids

    
    """
    def setUp(self):
        self.u = make_Universe(('resids', 'resnames', 'segids'))

    def tearDown(self):
        del self.u

    @staticmethod
    def assert_atoms_match_residue(atom, residue):
        # check that an Atom's residue level properties match the residue
        if isinstance(atom, groups.Atom):
            atom = [atom]  # can work with either Atom or AG
        if isinstance(residue, groups.Residue):
            residue = itertools.cycle((residue,))  # either R or RG or [R, R]
        for at, r in zip(atom, residue):
            assert_(at.resindex == r.resindex)
            assert_(at.resid == r.resid)
            assert_(at.resname == r.resname)

    def test_move_atom(self):
        # move a single atom by providing a new Residue object
        at = self.u.atoms[0]
        source = self.u.residues[0]
        dest = self.u.residues[4]

        assert_(at in source.atoms)
        assert_(not at in dest.atoms)
        self.assert_atoms_match_residue(at, source)
        assert_equal(len(source.atoms), 5)
        assert_equal(len(dest.atoms), 5)

        at.residue = dest

        assert_(not at in source.atoms)
        assert_(at in dest.atoms)
        self.assert_atoms_match_residue(at, dest)
        assert_equal(len(source.atoms), 4)
        assert_equal(len(dest.atoms), 6)

    def test_move_atomgroup_single_residue(self):
        # move contents of AtomGroup into Residue object
        ag = self.u.atoms[[1, 3]]
        source = self.u.residues[0]
        dest = self.u.residues[4]

        for at in ag:
            assert_(at in source.atoms)
            assert_(not at in dest.atoms)
        self.assert_atoms_match_residue(ag, source)
        assert_equal(len(source.atoms), 5)
        assert_equal(len(dest.atoms), 5)

        ag.residues = dest

        for at in ag:
            assert_(not at in source.atoms)
            assert_(at in dest.atoms)
        self.assert_atoms_match_residue(ag, dest)
        assert_equal(len(source.atoms), 3)
        assert_equal(len(dest.atoms), 7)

    def test_move_atomgroup_residuegroup(self):
        # move contents of AtomGroup into Residue object
        ag = self.u.atoms[[1, 3]]
        source = self.u.residues[0]
        dest = self.u.residues[4] + self.u.residues[5]

        for at in ag:
            assert_(at in source.atoms)
            assert_(not at in dest.atoms)
        self.assert_atoms_match_residue(ag, source)
        assert_equal(len(source.atoms), 5)
        assert_equal(len(dest[0].atoms), 5)
        assert_equal(len(dest[1].atoms), 5)

        ag.residues = dest

        for at in ag:
            assert_(not at in source.atoms)
            assert_(at in dest.atoms)
        self.assert_atoms_match_residue(ag, dest)
        assert_equal(len(source.atoms), 3)
        assert_equal(len(dest[0].atoms), 6)
        assert_equal(len(dest[1].atoms), 6)

    def test_move_atomgroup_residue_list(self):
        # move contents of AtomGroup into Residue object
        ag = self.u.atoms[[1, 3]]
        source = self.u.residues[0]
        dest = [self.u.residues[4], self.u.residues[5]]

        for at, d in zip(ag, dest):
            assert_(at in source.atoms)
            assert_(not at in d.atoms)
        self.assert_atoms_match_residue(ag, source)
        assert_equal(len(source.atoms), 5)
        assert_equal(len(dest[0].atoms), 5)
        assert_equal(len(dest[1].atoms), 5)

        ag.residues = dest

        for at, d in zip(ag, dest):
            assert_(not at in source.atoms)
            assert_(at in d.atoms)
        self.assert_atoms_match_residue(ag, dest)
        assert_equal(len(source.atoms), 3)
        assert_equal(len(dest[0].atoms), 6)
        assert_equal(len(dest[1].atoms), 6)

    # Wrong size argument for these operations
    def test_move_atom_residuegroup_TE(self):
        assert_raises(TypeError,
                      setattr, self.u.atoms[0], 'residue', self.u.atoms[1:3])

    def test_move_atom_residue_list_TE(self):
        dest = [self.u.residues[1], self.u.residues[3]]
        assert_raises(TypeError,
                      setattr, self.u.atoms[0], 'residue', dest)
        
    def test_move_atomgroup_residuegroup_VE(self):
        ag = self.u.atoms[:2]
        dest = self.u.residues[5:10]

        assert_raises(ValueError, setattr, ag, 'residues', dest)

    def test_move_atomgroup_residue_list_VE(self):
        ag = self.u.atoms[:2]
        dest = [self.u.residues[0], self.u.residues[10], self.u.residues[15]]

        assert_raises(ValueError, setattr, ag, 'residues', dest)

    # Setting to non-Residue/ResidueGroup raises TE
    def test_move_atom_TE(self):
        assert_raises(TypeError,
                      setattr, self.u.atoms[0], 'residue', 14)

    def test_move_atomgroup_TE(self):
        assert_raises(TypeError,
                      setattr, self.u.atoms[:5], 'residues', 15)

    def test_move_atomgroup_list_TE(self):
        assert_raises(TypeError,
                      setattr, self.u.atoms[:5], 'residues', [14, 12])
        
    # Test illegal moves - Atom.segment can't be changed
    def test_move_atom_segment_NIE(self):
        assert_raises(NotImplementedError,
                      setattr, self.u.atoms[0], 'segment', self.u.segments[1])

    def test_move_atomgroup_segment_NIE(self):
        assert_raises(NotImplementedError,
                      setattr, self.u.atoms[:3], 'segments', self.u.segments[1])

    @staticmethod
    def assert_residue_matches_segment(res, seg):
        if isinstance(res, groups.Residue):
            res = [res]
        if isinstance(seg, groups.Segment):
            seg = itertools.cycle((seg,))
        for r, s in zip(res, seg):
            assert_(r.segindex == s.segindex)
            assert_(r.segid == s.segid)

    def test_move_residue(self):
        res = self.u.residues[0]
        source = self.u.segments[0]
        dest = self.u.segments[2]

        assert_(res in source.residues)
        assert_(not res in dest.residues)
        self.assert_residue_matches_segment(res, source)
        assert_equal(len(source.residues), 5)
        assert_equal(len(dest.residues), 5)

        res.segment = dest

        assert_(not res in source.residues)
        assert_(res in dest.residues)
        self.assert_residue_matches_segment(res, dest)
        assert_equal(len(source.residues), 4)
        assert_equal(len(dest.residues), 6)

    def test_move_residuegroup_single_segment(self):
        res = self.u.residues[[1, 3]]
        source = self.u.segments[0]
        dest = self.u.segments[2]

        for r in res:
            assert_(r in source.residues)
            assert_(not r in dest.residues)
        self.assert_residue_matches_segment(res, source)
        assert_equal(len(source.residues), 5)
        assert_equal(len(dest.residues), 5)

        res.segments = dest

        for r in res:
            assert_(not r in source.residues)
            assert_(r in dest.residues)
        self.assert_residue_matches_segment(res, dest)
        assert_equal(len(source.residues), 3)
        assert_equal(len(dest.residues), 7)

    def test_move_residuegroup_segmentgroup(self):
        res = self.u.residues[[1, 3]]
        source = self.u.segments[0]
        dest = self.u.segments[2] + self.u.segments[3]

        for r in res:
            assert_(r in source.residues)
            assert_(not r in dest.residues)
        self.assert_residue_matches_segment(res, source)
        assert_equal(len(source.residues), 5)
        assert_equal(len(dest[0].residues), 5)
        assert_equal(len(dest[1].residues), 5)

        res.segments = dest

        for r in res:
            assert_(not r in source.residues)
            assert_(r in dest.residues)
        self.assert_residue_matches_segment(res, dest)
        assert_equal(len(source.residues), 3)
        assert_equal(len(dest[0].residues), 6)
        assert_equal(len(dest[1].residues), 6)

    def test_move_residuegroup_segment_list(self):
        res = self.u.residues[[1, 3]]
        source = self.u.segments[0]
        dest = [self.u.segments[2], self.u.segments[3]]

        for r, d in zip(res, dest):
            assert_(r in source.residues)
            assert_(not r in d.residues)
        self.assert_residue_matches_segment(res, source)
        assert_equal(len(source.residues), 5)
        assert_equal(len(dest[0].residues), 5)
        assert_equal(len(dest[1].residues), 5)

        res.segments = dest

        for r, d in zip(res, dest):
            assert_(not r in source.residues)
            assert_(r in d.residues)
        self.assert_residue_matches_segment(res, dest)
        assert_equal(len(source.residues), 3)
        assert_equal(len(dest[0].residues), 6)
        assert_equal(len(dest[1].residues), 6)

    def test_move_residue_segmentgroup_TE(self):
        assert_raises(TypeError,
                      setattr, self.u.residues[0], 'segment', self.u.segments[:4])

    def test_move_residue_list_TE(self):
        dest = [self.u.segments[3], self.u.segments[4]]
        assert_raises(TypeError,
                      setattr, self.u.residues[0], 'segment', dest)

    def test_move_residuegroup_segmentgroup_VE(self):
        rg = self.u.residues[:3]
        sg = self.u.segments[1:]

        assert_raises(ValueError, setattr, rg, 'segments', sg)

    def test_move_residuegroup_list_VE(self):
        rg = self.u.residues[:2]
        sg = [self.u.segments[1], self.u.segments[2], self.u.segments[3]]

        assert_raises(ValueError, setattr, rg, 'segments', sg)

    def test_move_residue_TE(self):
        assert_raises(TypeError,
                      setattr, self.u.residues[0], 'segment', 1)

    def test_move_residuegroup_TE(self):
        assert_raises(TypeError,
                      setattr, self.u.residues[:3], 'segments', 4)

    def test_move_residuegroup_list_TE(self):
        assert_raises(TypeError,
                      setattr, self.u.residues[:3], 'segments', [1, 2, 3])


class TestDownshiftArrays(object):
    def setUp(self):
        # test for square and ragged shapes
        # square shapes sometimes simplify to 2d array
        # which is bad!
        self.square = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
        self.square_size = 3
        self.square_result = np.array([[0, 3, 6], [1, 4, 7], [2, 5, 8]])
        self.ragged = np.array([0, 1, 2, 2, 0, 1, 2, 0, 1, 2])
        self.ragged_size = 3
        self.ragged_result = np.array([[0, 4, 7], [1, 5, 8], [2, 3, 6, 9]])

    def tearDown(self):
        del self.square
        del self.square_result
        del self.ragged
        del self.ragged_result

    @staticmethod
    def assert_rows_match(a, b):
        for row_a, row_b in zip(a, b):
            assert_array_equal(row_a, row_b)
        
    # The array as a whole must be dtype object
    # While the subarrays must be integers
    def test_downshift_dtype_square(self):
        out = make_downshift_arrays(self.square, self.square_size)
        assert_(out.dtype == object)
        assert_(out[0].dtype == np.intp)

    def test_downshift_dtype_ragged(self):
        out = make_downshift_arrays(self.ragged, self.ragged_size)
        assert_(out.dtype == object)
        assert_(out[0].dtype == np.intp)

    # Check shape and size
    # Shape should be size N+1 as None is appended
    def test_shape_square(self):
        out = make_downshift_arrays(self.square, self.square_size)
        assert_(out.shape == (4,))
        assert_(out[-1] is None)

    def test_shape_ragged(self):
        out = make_downshift_arrays(self.ragged, self.ragged_size)
        assert_(out.shape == (4,))
        assert_(out[-1] is None)

    def test_contents_square(self):
        out = make_downshift_arrays(self.square, self.square_size)
        self.assert_rows_match(out, self.square_result)

    def test_contents_ragged(self):
        out = make_downshift_arrays(self.ragged, self.ragged_size)
        self.assert_rows_match(out, self.ragged_result)

    def test_missing_intra_values(self):
        out = make_downshift_arrays(
            np.array([0, 0, 2, 2, 3, 3]), 4)
        self.assert_rows_match(out,
                               np.array([np.array([0, 1]),
                                         np.array([], dtype=np.int),
                                         np.array([2, 3]),
                                         np.array([4, 5]),
                                         None], dtype=object))

    def test_missing_intra_values_2(self):
        out = make_downshift_arrays(
            np.array([0, 0, 3, 3, 4, 4]), 5)
        self.assert_rows_match(out,
                               np.array([np.array([0, 1]),
                                         np.array([], dtype=np.int),
                                         np.array([], dtype=np.int),
                                         np.array([2, 3]),
                                         np.array([4, 5]),
                                         None], dtype=object))

    def test_missing_end_values(self):
        out = make_downshift_arrays(np.array([0, 0, 1, 1, 2, 2]), 4)
        self.assert_rows_match(out,
                               np.array([np.array([0, 1]),
                                         np.array([2, 3]),
                                         np.array([4, 5]),
                                         np.array([], dtype=np.int),
                                         None], dtype=object))

    def test_missing_end_values_2(self):
        out = make_downshift_arrays(np.array([0, 0, 1, 1, 2, 2]), 6)
        self.assert_rows_match(out,
                               np.array([np.array([0, 1]),
                                         np.array([2, 3]),
                                         np.array([4, 5]),
                                         np.array([], dtype=np.int),
                                         np.array([], dtype=np.int),
                                         None], dtype=object))


class TestAddingResidues(object):
    """Tests for adding residues and segments to a Universe

    Adding Residues
    Adding Segments

    Adding Residue without an attr
    Adding Residue with ambiguous segment
    Adding Segment without an attr

    Adding Residue and moving atoms to it
    Adding Segment and moving residues to it
    """
    def test_add_segment_no_attrs(self):
        u = make_Universe()

        assert_(len(u.segments) == 5)
        s = u.add_Segment()
        assert_(isinstance(s, MDAnalysis.core.groups.Segment))
        assert_(len(u.segments) == 6)

    def test_add_residue_no_attrs(self):
        u = make_Universe()

        assert_(len(u.residues) == 25)
        assert_(len(u.segments[1].residues) == 5)

        res = u.add_Residue(segment=u.segments[1])

        assert_(len(u.residues) == 26)
        assert_(len(u.segments[1].residues) == 6)
        assert_(res in u.segments[1].residues)

    def test_add_residue_no_attrs_one_segment(self):
        u = make_Universe(extras=[], size=(125, 25, 1))

        res = u.add_Residue()

        assert_(isinstance(res, MDAnalysis.core.groups.Residue))

    def test_add_Residue_ambiguous_segment_NDE(self):
        u = make_Universe()

        assert_raises(NoDataError, u.add_Residue)
    
    def test_add_Residue_missing_attr_NDE(self):
        u = make_Universe(('resids',))

        assert_raises(NoDataError, u.add_Residue, segment=u.segments[0])

    def test_add_Residue_NDE_message(self):
        # check error message asks for missing attr
        u = make_Universe(('resnames', 'resids'))

        try:
            u.add_Residue(segment=u.segments[0], resid=42)
        except NoDataError as e:
            assert_('resname' in e[0])
        else:
            raise AssertionError

    def test_add_Residue_NDE_message_2(self):
        # multiple missing attrs, check all get mentioned in error
        u = make_Universe(('resnames', 'resids'))

        try:
            u.add_Residue(segment=u.segments[0])
        except NoDataError as e:
            assert_('resname' in e[0])
            assert_('resid' in e[0])
        else:
            raise AssertionError

    def test_add_Residue_with_attrs(self):
        u = make_Universe(('resnames', 'resids'))

        r_new = u.add_Residue(segment=u.segments[0], resid=4321, resname='New')

        assert_(r_new.resid == 4321)
        assert_(r_new.resname == 'New')

    def test_missing_attr_NDE_Segment(self):
        u = make_Universe(('segids',))

        assert_raises(NoDataError, u.add_Segment)

    def test_add_Segment_NDE_message(self):
        u = make_Universe(('segids',))

        try:
            u.add_Segment()
        except NoDataError as e:
            assert_('segid' in e[0])
        else:
            raise AssertionError

    def test_add_Segment_with_attr(self):
        u = make_Universe(('segids',))

        new_seg = u.add_Segment(segid='New')

        assert_(new_seg.segid == 'New')


class TestTopologyGuessed(object):
    def setUp(self):
        names = self.names = ta.Atomnames(np.array(['A', 'B', 'C'], dtype=object))
        types = self.types = ta.Atomtypes(np.array(['X', 'Y', 'Z'], dtype=object),
                                          guessed=True)
        resids = self.resids = ta.Resids(np.array([1]))
        resnames = self.resnames = ta.Resnames(np.array(['ABC'], dtype=object),
                                               guessed=True)
        self.top = Topology(n_atoms=3, n_res=1,
                            attrs=[names, types, resids, resnames])

    def tearDown(self):
        del self.names
        del self.types
        del self.resids
        del self.resnames
        del self.top

    def test_guessed(self):
        guessed = self.top.guessed_attributes

        assert_(self.types in guessed)
        assert_(self.resnames in guessed)
        assert_(not self.names in guessed)
        assert_(not self.resids in guessed)

    def test_read(self):
        read = self.top.read_attributes

        assert_(self.names in read)
        assert_(self.resids in read)
        assert_(not self.types in read)
        assert_(not self.resnames in read)


class TestTopologyCreation(object):
    @staticmethod
    def test_make_topology_no_attrs():
        # should still make attrs list when attrs=None
        top = Topology()

        assert_(hasattr(top, 'attrs'))
        assert_(isinstance(top.attrs, list))

    @staticmethod
    def test_resindex_VE():
        # wrong sized atom to residue array
        AR = np.arange(10)
        assert_raises(ValueError,
                      Topology, n_atoms=5, atom_resindex=AR)

    @staticmethod
    def test_segindex_VE():
        # wrong sized residue to segment array
        AR = np.arange(5)
        RS = np.arange(10)
        assert_raises(ValueError,
                      Topology, n_atoms=5, n_res=5, atom_resindex=AR, residue_segindex=RS)
