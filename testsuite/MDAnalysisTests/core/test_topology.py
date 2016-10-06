"""Tests for MDAnalysis.core.topology objects.

Should convert between indices (*ix)
Should work with both a single or an array of indices
"""
import six
from six.moves import zip
import itertools

from numpy.testing import (
    assert_,
    assert_equal,
    assert_array_equal,
    assert_raises,
)
import numpy as np

from MDAnalysisTests.core.groupbase import make_Universe

from MDAnalysis.core.topology import (
    Topology,
    TransTable,
    make_downshift_arrays,
)
from MDAnalysis.core import groups


class TestTopology(object):
    """Tests for Topology object.

    """
    # Reference data
    Ridx = np.array([0, 0, 2, 2, 1, 1, 3, 3, 1, 2])
    Sidx = np.array([0, 1, 1, 0])

    def setUp(self):
        self.top = Topology(10, 4, 2, attrs=[],
                            atom_resindex=self.Ridx,
                            residue_segindex=self.Sidx)

    def tearDown(self):
        del self.top


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
        for rix, six in zip(
                [np.array([0, 1]),
                 np.array([2, 1, 0]),
                 np.array([1, 1, 1])],
                [np.array([0, 1]),
                 np.array([1, 1, 0]),
                 np.array([1, 1, 1])]
        ):
            assert_array_equal(self.tt.residues2segments(rix), six)

    def test_s2r_1d(self):
        for six, rix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[0, 3, 1, 2],
                 [1, 2, 0, 3],
                 [1, 2, 1, 2]]
                ):
            assert_array_equal(self.tt.segments2residues_1d(six), rix)

    def test_s2r_2d(self):
        for six, rix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[[0, 3], [1, 2]],
                 [[1, 2], [0, 3]],
                 [[1, 2], [1, 2]]]
        ):
            answer = self.tt.segments2residues_2d(six)
            for a1, a2 in zip(answer, rix):
                assert_array_equal(a1, a2)

    def test_s2a_1d(self):
        for six, aix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[0, 1, 6, 7, 4, 5, 8, 2, 3, 9],
                 [4, 5, 8, 2, 3, 9, 0, 1, 6, 7],
                 [4, 5, 8, 2, 3, 9, 4, 5, 8, 2, 3, 9]],
        ):
            assert_array_equal(self.tt.segments2atoms_1d(six), aix)

    def test_s2a_2d(self):
        for six, aix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[[0, 1, 6, 7], [4, 5, 8, 2, 3, 9]],
                 [[4, 5, 8, 2, 3, 9], [0, 1, 6, 7]],
                 [[4, 5, 8, 2, 3, 9], [4, 5, 8, 2, 3, 9]]],
        ):
            answer = self.tt.segments2atoms_2d(six)
            for a1, a2 in zip(answer, aix):
                assert_array_equal(a1, a2)

    # Moving within transtable without resizes
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
        self.u = make_Universe('resids', 'resnames', 'segids')

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
                      self.u.residues[0], 'segment', 1)

    def test_move_residuegroup_TE(self):
        assert_raises(TypeError,
                      self.u.residues[:3], 'segments', 4)

    def test_move_residuegroup_list_TE(self):
        assert_raises(TypeError,
                      self.u.residues[:3], 'segments', [1, 2, 3])


class TestDownshiftArrays(object):
    def setUp(self):
        # test for square and ragged shapes
        # square shapes sometimes simplify to 2d array
        # which is bad!
        self.square = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
        self.square_result = np.array([[0, 3, 6], [1, 4, 7], [2, 5, 8]])
        self.ragged = np.array([0, 1, 2, 2, 0, 1, 2, 0, 1, 2])
        self.ragged_result = np.array([[0, 4, 7], [1, 5, 8], [2, 3, 6, 9]])

    def tearDown(self):
        del self.square
        del self.square_result
        del self.ragged
        del self.ragged_result

    # The array as a whole must be dtype object
    # While the subarrays must be integers
    def test_downshift_dtype_square(self):
        out = make_downshift_arrays(self.square)
        assert_(out.dtype == object)
        assert_(out[0].dtype == np.int64)

    def test_downshift_dtype_ragged(self):
        out = make_downshift_arrays(self.ragged)
        assert_(out.dtype == object)
        assert_(out[0].dtype == np.int64)

    # Check shape and size
    # Shape should be size N+1 as None is appended
    def test_shape_square(self):
        out = make_downshift_arrays(self.square)
        assert_(out.shape == (4,))
        assert_(out[-1] is None)

    def test_shape_ragged(self):
        out = make_downshift_arrays(self.ragged)
        assert_(out.shape == (4,))
        assert_(out[-1] is None)

    def test_contents_square(self):
        out = make_downshift_arrays(self.square)
        for row, ref in zip(out, self.square_result):
            assert_array_equal(row, ref)

    def test_contents_ragged(self):
        out = make_downshift_arrays(self.ragged)
        for row, ref in zip(out, self.ragged_result):
            assert_array_equal(row, ref)
