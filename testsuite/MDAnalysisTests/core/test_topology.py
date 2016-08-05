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

from MDAnalysisTests.datafiles import PSF, DCD

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


class TestResidueMoves(object):
    def setUp(self):
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    def move_atom(self):
        at = self.u.atoms[0]
        
        assert_equal(at.resid, 1)
        assert_equal(at.resname, 'MET')
        assert_equal(len(self.u.residues[0]), 19)
        assert_equal(len(self.u.residues[4]), 19)

        at.resid = 5

        assert_equal(at.resid, 5)
        assert_equal(at.resname, 'LEU')
        assert_equal(len(self.u.residues[0]), 18)
        assert_equal(len(self.u.residues[4]), 20)
