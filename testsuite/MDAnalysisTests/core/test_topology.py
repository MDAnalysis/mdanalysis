"""Tests for MDAnalysis.core.topology objects.

Should convert between indices (*ix)
Should work with both a single or an array of indices
"""
from numpy.testing import (
    assert_,
    assert_equal,
    assert_array_equal,
)
import numpy as np

from MDAnalysisTests.core.groupbase import make_Universe

from MDAnalysis.core.topology import Topology, TransTable


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
    def setUp(self):
        self.u = make_Universe()

    def tearDown(self):
        del self.u

    def test_move_atom_via_resid(self):
        # move an atom between residues based on resid
        # residue must already exist!
        at = self.u.atoms[0]
        
        assert_equal(at.resindex, 0)
        assert_equal(at.resid, 1)
        assert_equal(at.resname, 'AAA')
        assert_equal(len(self.u.residues[0].atoms), 5)
        assert_equal(len(self.u.residues[4].atoms), 5)

        at.resid = 5

        assert_equal(at.resindex, 4)
        assert_equal(at.resid, 5)
        assert_equal(at.resname, 'EEE')
        assert_equal(len(self.u.residues[0].atoms), 4)
        assert_equal(len(self.u.residues[4].atoms), 6)

    def test_move_atom_via_resindex(self):
        # move an atom between residues based on resindex
        # residue must already exist!
        at = self.u.atoms[0]
        
        assert_equal(at.resindex, 0)
        assert_equal(at.resname, 'AAA')
        assert_equal(len(self.u.residues[0].atoms), 5)
        assert_equal(len(self.u.residues[4].atoms), 5)

        at.resindex = 4

        assert_equal(at.resindex, 4)
        assert_equal(at.resname, 'EEE')
        assert_equal(len(self.u.residues[0].atoms), 4)
        assert_equal(len(self.u.residues[4].atoms), 6)

    def test_move_atomgroup_via_resid(self):
        # move some atoms between residues based on resid
        # residue must already exist!
        ag = self.u.atoms[[1, 3]]
        
        assert_equal(ag.resindices, 0)
        assert_equal(ag.resids, 1)
        assert_equal(ag.resnames, 'AAA')
        assert_equal(len(self.u.residues[0].atoms), 5)
        assert_equal(len(self.u.residues[4].atoms), 5)

        # TODO!
        # Should also be able to supply an array of len(ag)
        # to set atoms to different residues each
        # ie ag.resids = [4, 5]
        ag.resids = 5

        assert_equal(ag.resindices, 4)
        assert_equal(ag.resids, 5)
        assert_equal(ag.resnames, 'EEE')
        assert_equal(len(self.u.residues[0].atoms), 3)
        assert_equal(len(self.u.residues[4].atoms), 7)

    def test_move_atomgroup_via_resindex(self):
        # move some atoms between residues based on resindex
        # residue must already exist!
        ag = self.u.atoms[[1, 3]]
        
        assert_equal(ag.resindices, 0)
        assert_equal(ag.resnames, 'AAA')
        assert_equal(len(self.u.residues[0].atoms), 5)
        assert_equal(len(self.u.residues[4].atoms), 5)

        ag.resindices = 4

        assert_equal(ag.resindices, 4)
        assert_equal(ag.resnames, 'EEE')
        assert_equal(len(self.u.residues[0].atoms), 3)
        assert_equal(len(self.u.residues[4].atoms), 7)

    def test_move_residue_via_segindex(self):
        res = self.u.residues[0]

        assert_equal(res.segindex, 0)
        assert_equal(res.segid, 'SegA')
        assert_equal(len(self.u.segments[0].residues), 5)
        assert_equal(len(self.u.segments[1].residues), 5)

        res.segindex = 1

        assert_equal(res.segindex, 1)
        assert_equal(res.segid, 'SegB')
        assert_equal(len(self.u.segments[0].residues), 4)
        assert_equal(len(self.u.segments[1].residues), 6)

    def test_move_residue_via_segid(self):
        res = self.u.residues[0]

        assert_equal(res.segindex, 0)
        assert_equal(res.segid, 'SegA')
        assert_equal(len(self.u.segments[0].residues), 5)
        assert_equal(len(self.u.segments[1].residues), 5)

        res.segid = 'SegB'

        assert_equal(res.segindex, 1)
        assert_equal(res.segid, 'SegB')
        assert_equal(len(self.u.segments[0].residues), 4)
        assert_equal(len(self.u.segments[1].residues), 6)

    def test_move_residuegroup_via_segindex(self):
        rg = self.u.residues[[1, 3]]

        assert_equal(rg.segindices, 0)
        assert_equal(rg.segids, 'SegA')
        assert_equal(len(self.u.segments[0].residues), 5)
        assert_equal(len(self.u.segments[1].residues), 5)

        # TODO: Provide array ie [1, 2]
        rg.segindices = 1

        assert_equal(rg.segindices, 1)
        assert_equal(rg.segids, 'SegB')
        assert_equal(len(self.u.segments[0].residues), 3)
        assert_equal(len(self.u.segments[1].residues), 7)


    def test_move_residuegroup_via_segid(self):
        rg = self.u.residues[[1, 3]]

        assert_equal(rg.segindices, 0)
        assert_equal(rg.segids, 'SegA')
        assert_equal(len(self.u.segments[0].residues), 5)
        assert_equal(len(self.u.segments[1].residues), 5)

        rg.segids = 'SegB'

        assert_equal(rg.segindices, 1)
        assert_equal(rg.segids, 'SegB')
        assert_equal(len(self.u.segments[0].residues), 3)
        assert_equal(len(self.u.segments[1].residues), 7)


