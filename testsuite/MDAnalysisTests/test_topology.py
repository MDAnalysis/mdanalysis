"""Tests for MDAnalysis.core.topology objects.

Should convert between indices (*ix)
Should work with both a single or an array of indices
"""
from numpy.testing import (
    assert_,
    assert_array_equal,
)
import numpy as np

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
    """Tests for Atom/Residue/Segment index translation.
    
    Should convert between indices (*ix)
    Should work with both a single or an array of indices
    """
    # Reference data
    Ridx = np.array([0, 0, 2, 2, 1, 1, 3, 3, 1, 2])
    Sidx = np.array([0, 1, 1, 0])

    def setUp(self):
        self.tt = TransTable(10, 4, 2, self.Ridx, self.Sidx)

    def tearDown(self):
        del self.tt

    def test_a2r_single(self):
        for ai, ri in enumerate(self.Ridx):
            assert_(self.tt.a2r(ai) == ri)

    def test_a2r_many(self):
        for aix, rix in zip(
                [np.array([0, 1, 2]),
                 np.array([9, 6, 2]),
                 np.array([3, 3, 3])],
                [np.array([0, 0, 2]),
                 np.array([2, 3, 2]),
                 np.array([2, 2, 2])]
        ):
            assert_array_equal(self.tt.a2r(aix), rix)

    def test_r2a_1d_single(self):
        for ri, aix in zip(
                [0, 1, 2, 3],
                [[0, 1], [4, 5, 8], [2, 3, 9], [6, 7]]
        ):
            assert_array_equal(self.tt.r2a_1d(ri), aix)

    def test_r2a_1d_many(self):
        for rix, aix in zip(
                [[0, 1], [1, 1], [3, 1]],
                [[0, 1, 4, 5, 8], [4, 5, 8, 4, 5, 8], [6, 7, 4, 5, 8]]
        ):
            assert_array_equal(self.tt.r2a_1d(rix), aix)

    def test_r2a_2d_many(self):
        for rix, aix in zip(
                [[0, 1],
                 [1, 1],
                 [3, 1]],
                [[[0, 1], [4, 5, 8]],
                 [[4, 5, 8], [4, 5, 8]],
                 [[6, 7], [4, 5, 8]]]
        ):
            answer = self.tt.r2a_2d(rix)
            for a1, a2 in zip(answer, aix):
                assert_array_equal(a1, a2)


    def test_r2s_single(self):
        for ri, si in enumerate(self.Sidx):
            assert_(self.tt.r2s(ri) == si)

    def test_r2s_many(self):
        for rix, six in zip(
                [np.array([0, 1]),
                 np.array([2, 1, 0]),
                 np.array([1, 1, 1])],
                [np.array([0, 1]),
                 np.array([1, 1, 0]),
                 np.array([1, 1, 1])]
        ):
            assert_array_equal(self.tt.r2s(rix), six)

    def test_s2r_1d_single(self):
        for si, rix in zip(
                [0, 1],
                [[0, 3], [1, 2]]
        ):
            assert_array_equal(self.tt.s2r_1d(si), rix)

    def test_s2r_1d_many(self):
        for six, rix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[0, 3, 1, 2],
                 [1, 2, 0, 3],
                 [1, 2, 1, 2]]
                ):
            assert_array_equal(self.tt.s2r_1d(six), rix)

    def test_s2r_2d_single(self):
        for si, rix in zip(
                [0, 1],
                [[[0, 3]], [[1, 2]]]
        ):
            answer = self.tt.s2r_2d(si)
            for a1, a2 in zip(answer, rix):
                assert_array_equal(a1, a2)

    def test_s2r_2d_many(self):
        for six, rix in zip(
                [[0, 1],
                 [1, 0],
                 [1, 1]],
                [[[0, 3], [1, 2]],
                 [[1, 2], [0, 3]],
                 [[1, 2], [1, 2]]]
        ):
            answer = self.tt.s2r_2d(six)
            for a1, a2 in zip(answer, rix):
                assert_array_equal(a1, a2)

