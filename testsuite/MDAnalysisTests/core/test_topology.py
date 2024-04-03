"""Tests for MDAnalysis.core.topology objects.

Should convert between indices (*ix)
Should work with both a single or an array of indices
"""
import itertools
from numpy.testing import (
    assert_equal,
)
import pytest
import numpy as np
import pickle

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


def assert_rows_match(a, b):
    for row_a, row_b in zip(a, b):
        assert_equal(row_a, row_b)


class TestTransTable(object):
    @pytest.fixture()
    def tt(self):
        Ridx = np.array([0, 0, 2, 2, 1, 1, 3, 3, 1, 2])
        Sidx = np.array([0, 1, 1, 0])
        return TransTable(10, 4, 2, Ridx, Sidx)

    def test_a2r(self, tt):
        for aix, rix in zip(
                [np.array([0, 1, 2]),
                 np.array([9, 6, 2]),
                 np.array([3, 3, 3])],
                [np.array([0, 0, 2]),
                 np.array([2, 3, 2]),
                 np.array([2, 2, 2])]
        ):
            assert_equal(tt.atoms2residues(aix), rix)

    def test_r2a_1d(self, tt):
        for rix, aix in zip(
                [[0, 1], [1, 1], [3, 1]],
                [[0, 1, 4, 5, 8], [4, 5, 8, 4, 5, 8], [6, 7, 4, 5, 8]]
        ):
            assert_equal(tt.residues2atoms_1d(rix), aix)

    def test_r2a_2d(self, tt):
        for rix, aix in zip(
                [[0, 1],
                 [1, 1],
                 [3, 1]],
                [[[0, 1], [4, 5, 8]],
                 [[4, 5, 8], [4, 5, 8]],
                 [[6, 7], [4, 5, 8]]]
        ):
            answer = tt.residues2atoms_2d(rix)
            for a1, a2 in zip(answer, aix):
                assert_equal(a1, a2)

    def test_r2s(self, tt):
        for rix, sidx in zip(
                [np.array([0, 1]),
                 np.array([2, 1, 0]),
                 np.array([1, 1, 1])],
                [np.array([0, 1]),
                 np.array([1, 1, 0]),
                 np.array([1, 1, 1])]
        ):
            assert_equal(tt.residues2segments(rix), sidx)

    def test_s2r_1d(self, tt):
        for sidx, rix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[0, 3, 1, 2],
                 [1, 2, 0, 3],
                 [1, 2, 1, 2]]
        ):
            assert_equal(tt.segments2residues_1d(sidx), rix)

    def test_s2r_2d(self, tt):
        for sidx, rix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[[0, 3], [1, 2]],
                 [[1, 2], [0, 3]],
                 [[1, 2], [1, 2]]]
        ):
            answer = tt.segments2residues_2d(sidx)
            for a1, a2 in zip(answer, rix):
                assert_equal(a1, a2)

    def test_s2a_1d(self, tt):
        for sidx, aix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[0, 1, 6, 7, 4, 5, 8, 2, 3, 9],
                 [4, 5, 8, 2, 3, 9, 0, 1, 6, 7],
                 [4, 5, 8, 2, 3, 9, 4, 5, 8, 2, 3, 9]],
        ):
            assert_equal(tt.segments2atoms_1d(sidx), aix)

    def test_s2a_2d(self, tt):
        for sidx, aix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[[0, 1, 6, 7], [4, 5, 8, 2, 3, 9]],
                 [[4, 5, 8, 2, 3, 9], [0, 1, 6, 7]],
                 [[4, 5, 8, 2, 3, 9], [4, 5, 8, 2, 3, 9]]],
        ):
            answer = tt.segments2atoms_2d(sidx)
            for a1, a2 in zip(answer, aix):
                assert_equal(a1, a2)

    def test_move_atom_simple(self, tt):
        assert_equal(tt.atoms2residues(1), 0)
        assert_equal(len(tt.residues2atoms_1d(0)), 2)
        assert_equal(len(tt.residues2atoms_1d(3)), 2)

        # move 2nd atom to 4th residue (atom1 -> res3)
        tt.move_atom(1, 3)

        assert_equal(tt.atoms2residues(1), 3)  # identity of parent changed
        assert_equal(len(tt.residues2atoms_1d(0)), 1)  # 1 fewer here
        assert_equal(len(tt.residues2atoms_1d(3)), 3)  # 1 more here

    def test_move_residue_simple(self, tt):
        assert_equal(tt.residues2segments(1), 1)
        assert_equal(len(tt.segments2residues_1d(0)), 2)
        assert_equal(len(tt.segments2residues_1d(1)), 2)

        # move 2nd residue to 1st segment (res1 -> seg0)
        tt.move_residue(1, 0)

        assert_equal(tt.residues2segments(1), 0)
        assert_equal(len(tt.segments2residues_1d(0)), 3)
        assert_equal(len(tt.segments2residues_1d(1)), 1)

    def test_lazy_building_RA(self, tt):
        assert_equal(tt._RA, None)
        RA = tt.RA
        assert_rows_match(tt.RA,
                          np.array([np.array([0, 1]),
                                    np.array([4, 5, 8]),
                                    np.array([2, 3, 9]),
                                    np.array([6, 7]),
                                    None], dtype=object))

        tt.move_atom(1, 3)
        assert_equal(tt._RA, None)

    def test_lazy_building_SR(self, tt):
        assert_equal(tt._SR, None)
        SR = tt.SR
        assert_rows_match(tt.SR,
                          np.array([np.array([0, 3]),
                                    np.array([1, 2]),
                                    None], dtype=object))

        tt.move_residue(1, 0)
        assert_equal(tt._SR, None)

    def test_serialization(self, tt):
        _ = tt.RA
        _ = tt.SR
        tt_loaded = pickle.loads(pickle.dumps(tt))
        assert_equal(tt_loaded._RA, None)
        assert_equal(tt_loaded._SR, None)
        assert_rows_match(tt_loaded.RA, tt.RA)
        assert_rows_match(tt_loaded.SR, tt.SR)


class TestLevelMoves(object):
    """Tests for moving atoms/residues between residues/segments

    
    Atoms can move between residues by setting .residue with a Residue
    Residues can move between segments by setting .segment with a Segment

    Moves are performed by setting either [res/seg]indices or [res/seg]ids

    
    """

    @pytest.fixture()
    def u(self):
        return make_Universe(('resids', 'resnames', 'segids'))

    @staticmethod
    def assert_atoms_match_residue(atom, residue):
        # check that an Atom's residue level properties match the residue
        if isinstance(atom, groups.Atom):
            atom = [atom]  # can work with either Atom or AG
        if isinstance(residue, groups.Residue):
            residue = itertools.cycle((residue,))  # either R or RG or [R, R]
        for at, r in zip(atom, residue):
            assert at.resindex == r.resindex
            assert at.resid == r.resid
            assert at.resname == r.resname

    def test_move_atom(self, u):
        # move a single atom by providing a new Residue object
        at = u.atoms[0]
        source = u.residues[0]
        dest = u.residues[4]

        assert at in source.atoms
        assert not at in dest.atoms
        self.assert_atoms_match_residue(at, source)
        assert_equal(len(source.atoms), 5)
        assert_equal(len(dest.atoms), 5)

        at.residue = dest

        assert not at in source.atoms
        assert at in dest.atoms
        self.assert_atoms_match_residue(at, dest)
        assert_equal(len(source.atoms), 4)
        assert_equal(len(dest.atoms), 6)

    def test_move_atomgroup_single_residue(self, u):
        # move contents of AtomGroup into Residue object
        ag = u.atoms[[1, 3]]
        source = u.residues[0]
        dest = u.residues[4]

        for at in ag:
            assert at in source.atoms
            assert not at in dest.atoms
        self.assert_atoms_match_residue(ag, source)
        assert_equal(len(source.atoms), 5)
        assert_equal(len(dest.atoms), 5)

        ag.residues = dest

        for at in ag:
            assert not at in source.atoms
            assert at in dest.atoms
        self.assert_atoms_match_residue(ag, dest)
        assert_equal(len(source.atoms), 3)
        assert_equal(len(dest.atoms), 7)

    def test_move_atomgroup_residuegroup(self, u):
        # move contents of AtomGroup into Residue object
        ag = u.atoms[[1, 3]]
        source = u.residues[0]
        dest = u.residues[4] + u.residues[5]

        for at in ag:
            assert at in source.atoms
            assert not at in dest.atoms
        self.assert_atoms_match_residue(ag, source)
        assert_equal(len(source.atoms), 5)
        assert_equal(len(dest[0].atoms), 5)
        assert_equal(len(dest[1].atoms), 5)

        ag.residues = dest

        for at in ag:
            assert not at in source.atoms
            assert at in dest.atoms
        self.assert_atoms_match_residue(ag, dest)
        assert_equal(len(source.atoms), 3)
        assert_equal(len(dest[0].atoms), 6)
        assert_equal(len(dest[1].atoms), 6)

    def test_move_atomgroup_residue_list(self, u):
        # move contents of AtomGroup into Residue object
        ag = u.atoms[[1, 3]]
        source = u.residues[0]
        dest = [u.residues[4], u.residues[5]]

        for at, d in zip(ag, dest):
            assert at in source.atoms
            assert not at in d.atoms
        self.assert_atoms_match_residue(ag, source)
        assert_equal(len(source.atoms), 5)
        assert_equal(len(dest[0].atoms), 5)
        assert_equal(len(dest[1].atoms), 5)

        ag.residues = dest

        for at, d in zip(ag, dest):
            assert not at in source.atoms
            assert at in d.atoms
        self.assert_atoms_match_residue(ag, dest)
        assert_equal(len(source.atoms), 3)
        assert_equal(len(dest[0].atoms), 6)
        assert_equal(len(dest[1].atoms), 6)

    # Wrong size argument for these operations
    def test_move_atom_residuegroup_TE(self, u):
        with pytest.raises(TypeError):
            setattr(u.atoms[0], 'residue', u.atoms[1:3])

    def test_move_atom_residue_list_TE(self, u):
        dest = [u.residues[1], u.residues[3]]
        with pytest.raises(TypeError):
            setattr(u.atoms[0], 'residue', dest)

    def test_move_atomgroup_residuegroup_VE(self, u):
        ag = u.atoms[:2]
        dest = u.residues[5:10]

        with pytest.raises(ValueError):
            setattr(ag, 'residues', dest)

    def test_move_atomgroup_residue_list_VE(self, u):
        ag = u.atoms[:2]
        dest = [u.residues[0], u.residues[10], u.residues[15]]

        with pytest.raises(ValueError):
            setattr(ag, 'residues', dest)

    # Setting to non-Residue/ResidueGroup raises TE
    def test_move_atom_TE(self, u):
        with pytest.raises(TypeError):
            setattr(u.atoms[0], 'residue', 14)

    def test_move_atomgroup_TE(self, u):
        with pytest.raises(TypeError):
            setattr(u.atoms[:5], 'residues', 15)

    def test_move_atomgroup_list_TE(self, u):
        with pytest.raises(TypeError):
            setattr(u.atoms[:5], 'residues', [14, 12])

    # Test illegal moves - Atom.segment can't be changed
    def test_move_atom_segment_NIE(self, u):
        with pytest.raises(NotImplementedError):
            setattr(u.atoms[0], 'segment', u.segments[1])

    def test_move_atomgroup_segment_NIE(self, u):
        with pytest.raises(NotImplementedError):
            setattr(u.atoms[:3], 'segments', u.segments[1])

    @staticmethod
    def assert_residue_matches_segment(res, seg):
        if isinstance(res, groups.Residue):
            res = [res]
        if isinstance(seg, groups.Segment):
            seg = itertools.cycle((seg,))
        for r, s in zip(res, seg):
            assert r.segindex == s.segindex
            assert r.segid == s.segid

    def test_move_residue(self, u):
        res = u.residues[0]
        source = u.segments[0]
        dest = u.segments[2]

        assert res in source.residues
        assert not res in dest.residues
        self.assert_residue_matches_segment(res, source)
        assert_equal(len(source.residues), 5)
        assert_equal(len(dest.residues), 5)

        res.segment = dest

        assert not res in source.residues
        assert res in dest.residues
        self.assert_residue_matches_segment(res, dest)
        assert_equal(len(source.residues), 4)
        assert_equal(len(dest.residues), 6)

    def test_move_residuegroup_single_segment(self, u):
        res = u.residues[[1, 3]]
        source = u.segments[0]
        dest = u.segments[2]

        for r in res:
            assert r in source.residues
            assert not r in dest.residues
        self.assert_residue_matches_segment(res, source)
        assert_equal(len(source.residues), 5)
        assert_equal(len(dest.residues), 5)

        res.segments = dest

        for r in res:
            assert not r in source.residues
            assert r in dest.residues
        self.assert_residue_matches_segment(res, dest)
        assert_equal(len(source.residues), 3)
        assert_equal(len(dest.residues), 7)

    def test_move_residuegroup_segmentgroup(self, u):
        res = u.residues[[1, 3]]
        source = u.segments[0]
        dest = u.segments[2] + u.segments[3]

        for r in res:
            assert r in source.residues
            assert not r in dest.residues
        self.assert_residue_matches_segment(res, source)
        assert_equal(len(source.residues), 5)
        assert_equal(len(dest[0].residues), 5)
        assert_equal(len(dest[1].residues), 5)

        res.segments = dest

        for r in res:
            assert not r in source.residues
            assert r in dest.residues
        self.assert_residue_matches_segment(res, dest)
        assert_equal(len(source.residues), 3)
        assert_equal(len(dest[0].residues), 6)
        assert_equal(len(dest[1].residues), 6)

    def test_move_residuegroup_segment_list(self, u):
        res = u.residues[[1, 3]]
        source = u.segments[0]
        dest = [u.segments[2], u.segments[3]]

        for r, d in zip(res, dest):
            assert r in source.residues
            assert not r in d.residues
        self.assert_residue_matches_segment(res, source)
        assert_equal(len(source.residues), 5)
        assert_equal(len(dest[0].residues), 5)
        assert_equal(len(dest[1].residues), 5)

        res.segments = dest

        for r, d in zip(res, dest):
            assert not r in source.residues
            assert r in d.residues
        self.assert_residue_matches_segment(res, dest)
        assert_equal(len(source.residues), 3)
        assert_equal(len(dest[0].residues), 6)
        assert_equal(len(dest[1].residues), 6)

    def test_move_residue_segmentgroup_TE(self, u):
        with pytest.raises(TypeError):
            setattr(u.residues[0], 'segment', u.segments[:4])

    def test_move_residue_list_TE(self, u):
        dest = [u.segments[3], u.segments[4]]
        with pytest.raises(TypeError):
            setattr(u.residues[0], 'segment', dest)

    def test_move_residuegroup_segmentgroup_VE(self, u):
        rg = u.residues[:3]
        sg = u.segments[1:]

        with pytest.raises(ValueError):
            setattr(rg, 'segments', sg)

    def test_move_residuegroup_list_VE(self, u):
        rg = u.residues[:2]
        sg = [u.segments[1], u.segments[2], u.segments[3]]

        with pytest.raises(ValueError):
            setattr(rg, 'segments', sg)

    def test_move_residue_TE(self, u):
        with pytest.raises(TypeError):
            setattr(u.residues[0], 'segment', 1)

    def test_move_residuegroup_TE(self, u):
        with pytest.raises(TypeError):
            setattr(u.residues[:3], 'segments', 4)

    def test_move_residuegroup_list_TE(self, u):
        with pytest.raises(TypeError):
            setattr(u.residues[:3], 'segments', [1, 2, 3])


class TestDownshiftArrays(object):
    @pytest.fixture()
    def square(self):
        return np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])

    @pytest.fixture()
    def square_size(self):
        return 3

    @pytest.fixture()
    def square_result(self):
        return np.array([[0, 3, 6], [1, 4, 7], [2, 5, 8]])

    @pytest.fixture()
    def ragged(self):
        return np.array([0, 1, 2, 2, 0, 1, 2, 0, 1, 2])

    @pytest.fixture()
    def ragged_size(self):
        return 3

    @pytest.fixture()
    def ragged_result(self):
        return np.array([[0, 4, 7], [1, 5, 8], [2, 3, 6, 9]],
                        dtype=object)

    # The array as a whole must be dtype object
    # While the subarrays must be integers
    def test_downshift_dtype_square(self, square, square_size):
        out = make_downshift_arrays(square, square_size)
        assert out.dtype == object
        assert out[0].dtype == np.intp

    def test_downshift_dtype_ragged(self, ragged, ragged_size):
        out = make_downshift_arrays(ragged, ragged_size)
        assert out.dtype == object
        assert out[0].dtype == np.intp

    # Check shape and size
    # Shape should be size N+1 as None is appended
    def test_shape_square(self, square, square_size):
        out = make_downshift_arrays(square, square_size)
        assert out.shape == (4,)
        assert out[-1] is None

    def test_shape_ragged(self, ragged, ragged_size):
        out = make_downshift_arrays(ragged, ragged_size)
        assert out.shape == (4,)
        assert out[-1] is None

    def test_contents_square(self, square, square_size, square_result):
        out = make_downshift_arrays(square, square_size)
        assert_rows_match(out, square_result)

    def test_contents_ragged(self, ragged, ragged_size, ragged_result):
        out = make_downshift_arrays(ragged, ragged_size)
        assert_rows_match(out, ragged_result)

    def test_missing_intra_values(self):
        out = make_downshift_arrays(
            np.array([0, 0, 2, 2, 3, 3]), 4)
        assert_rows_match(out,
                               np.array([np.array([0, 1]),
                                         np.array([], dtype=int),
                                         np.array([2, 3]),
                                         np.array([4, 5]),
                                         None], dtype=object))

    def test_missing_intra_values_2(self):
        out = make_downshift_arrays(
            np.array([0, 0, 3, 3, 4, 4]), 5)
        assert_rows_match(out,
                               np.array([np.array([0, 1]),
                                         np.array([], dtype=int),
                                         np.array([], dtype=int),
                                         np.array([2, 3]),
                                         np.array([4, 5]),
                                         None], dtype=object))

    def test_missing_end_values(self):
        out = make_downshift_arrays(np.array([0, 0, 1, 1, 2, 2]), 4)
        assert_rows_match(out,
                               np.array([np.array([0, 1]),
                                         np.array([2, 3]),
                                         np.array([4, 5]),
                                         np.array([], dtype=int),
                                         None], dtype=object))

    def test_missing_end_values_2(self):
        out = make_downshift_arrays(np.array([0, 0, 1, 1, 2, 2]), 6)
        assert_rows_match(out,
                               np.array([np.array([0, 1]),
                                         np.array([2, 3]),
                                         np.array([4, 5]),
                                         np.array([], dtype=int),
                                         np.array([], dtype=int),
                                         None], dtype=object))

    def test_missing_start_values_2(self):
        out = make_downshift_arrays(np.array([1, 1, 2, 2, 3, 3]), 4)
        assert_rows_match(out,
                          np.array([np.array([], dtype=int),
                                    np.array([0, 1]),
                                    np.array([2, 3]),
                                    np.array([4, 5]),
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

        assert len(u.segments) == 5
        s = u.add_Segment()
        assert isinstance(s, MDAnalysis.core.groups.Segment)
        assert len(u.segments) == 6

    def test_add_residue_no_attrs(self):
        u = make_Universe()

        assert len(u.residues) == 25
        assert len(u.segments[1].residues) == 5

        res = u.add_Residue(segment=u.segments[1])

        assert len(u.residues) == 26
        assert len(u.segments[1].residues) == 6
        assert res in u.segments[1].residues

    def test_add_residue_no_attrs_one_segment(self):
        u = make_Universe(extras=[], size=(125, 25, 1))

        res = u.add_Residue()

        assert isinstance(res, MDAnalysis.core.groups.Residue)

    def test_add_Residue_ambiguous_segment_NDE(self):
        u = make_Universe()

        with pytest.raises(NoDataError):
            u.add_Residue()

    def test_add_Residue_missing_attr_NDE(self):
        u = make_Universe(('resids',))

        with pytest.raises(NoDataError):
            u.add_Residue(segment=u.segments[0])

    def test_add_Residue_NDE_message(self):
        # check error message asks for missing attr
        u = make_Universe(('resnames', 'resids'))

        try:
            u.add_Residue(segment=u.segments[0], resid=42)
        except NoDataError as e:
            assert 'resname' in str(e)
        else:
            raise AssertionError

    def test_add_Residue_NDE_message_2(self):
        # multiple missing attrs, check all get mentioned in error
        u = make_Universe(('resnames', 'resids'))

        try:
            u.add_Residue(segment=u.segments[0])
        except NoDataError as e:
            assert 'resname' in str(e)
            assert 'resid' in str(e)
        else:
            raise AssertionError

    def test_add_Residue_with_attrs(self):
        u = make_Universe(('resnames', 'resids'))

        r_new = u.add_Residue(segment=u.segments[0], resid=4321, resname='New')

        assert r_new.resid == 4321
        assert r_new.resname == 'New'

    def test_missing_attr_NDE_Segment(self):
        u = make_Universe(('segids',))
        with pytest.raises(NoDataError):
            u.add_Segment()

    def test_add_Segment_NDE_message(self):
        u = make_Universe(('segids',))

        try:
            u.add_Segment()
        except NoDataError as e:
            assert 'segid' in str(e)
        else:
            raise AssertionError

    def test_add_Segment_with_attr(self):
        u = make_Universe(('segids',))

        new_seg = u.add_Segment(segid='New')

        assert new_seg.segid == 'New'


class TestTopologyGuessed(object):
    @pytest.fixture()
    def names(self):
        return ta.Atomnames(np.array(['A', 'B', 'C'], dtype=object))

    @pytest.fixture()
    def types(self):
        return ta.Atomtypes(np.array(['X', 'Y', 'Z'], dtype=object),
                            guessed=True)

    @pytest.fixture()
    def resids(self):
        return ta.Resids(np.array([1]))

    @pytest.fixture()
    def resnames(self):
        return ta.Resnames(np.array(['ABC'], dtype=object),
                           guessed=True)

    @pytest.fixture()
    def bonds(self):
        return ta.Bonds([(1, 2), (2, 3)], guessed=False)

    @pytest.fixture()
    def top(self, names, types, resids, resnames, bonds):
        return Topology(n_atoms=3, n_res=1,
                        attrs=[names, types, resids, resnames, bonds])

    def test_guessed(self, names, types, resids, resnames, bonds, top):
        guessed = top.guessed_attributes

        assert types in guessed
        assert resnames in guessed
        assert not names in guessed
        assert not resids in guessed
        assert bonds not in guessed

    def test_read(self, names, types, resids, resnames, bonds, top):
        read = top.read_attributes
        assert names in read
        assert resids in read
        assert bonds in read
        assert not types in read
        assert not resnames in read

class TestTopologyCreation(object):
    def test_make_topology_no_attrs(self):
        # should still make attrs list when attrs=None
        top = Topology()

        assert hasattr(top, 'attrs')
        assert isinstance(top.attrs, list)

    def test_resindex_VE(self):
        # wrong sized atom to residue array
        AR = np.arange(10)
        with pytest.raises(ValueError):
            Topology(n_atoms=5, atom_resindex=AR)

    def test_segindex_VE(self):
        # wrong sized residue to segment array
        AR = np.arange(5)
        RS = np.arange(10)
        with pytest.raises(ValueError):
            Topology(n_atoms=5, n_res=5, atom_resindex=AR,
                     residue_segindex=RS)
