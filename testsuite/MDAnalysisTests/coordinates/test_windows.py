# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""Tests for DOS style line endings

The files here are written from a Windows system, so have DOS style line endings.

This was previously problematic for ascii files.

These tests are just basic repeats of the tests from elsewhere, but on these DOS
line ending files.
"""

import numpy as np
import pytest
from numpy.testing import assert_almost_equal, assert_equal
from MDAnalysisTests.datafiles import (
    WIN_PDB_multiframe,
    WIN_ARC,
    WIN_DLP_HISTORY,
    WIN_TRJ,
    WIN_LAMMPSDUMP,
)
from MDAnalysisTests.coordinates.test_trj import TestTRJReader
from MDAnalysisTests.coordinates.test_dlpoly import TestDLPolyHistory
from MDAnalysisTests.coordinates.test_lammps import TestLammpsDumpReader
import MDAnalysis as mda


class TestWinLammpsDump(TestLammpsDumpReader):
    @pytest.fixture
    def u(self):
        return mda.Universe(WIN_LAMMPSDUMP, format='LAMMPSDUMP')


class TestWinPDB(object):
    @staticmethod
    @pytest.fixture(scope='class')
    def multiverse():
        return mda.Universe(WIN_PDB_multiframe, guess_bonds=True)

    def test_n_frames(self, multiverse):
        assert_equal(multiverse.trajectory.n_frames, 24,
                     "Wrong number of frames read from PDB muliple model file")
    def test_rewind(self, multiverse):
        u = multiverse
        u.trajectory[11]
        assert_equal(u.trajectory.ts.frame, 11,
                     "Failed to forward to 11th frame (frame index 11)")
        u.trajectory.rewind()
        assert_equal(u.trajectory.ts.frame, 0,
                     "Failed to rewind to 0th frame (frame index 0)")

    def test_iteration(self, multiverse):
        u = multiverse
        frames = []
        for frame in u.trajectory:
            pass
        # should rewind after previous test
        # problem was: the iterator is NoneType and next() cannot be called
        for ts in u.trajectory:
            frames.append(ts)
        assert_equal(
            len(frames), u.trajectory.n_frames,
            "iterated number of frames %d is not the expected number %d; "
            "trajectory iterator fails to rewind" %
            (len(frames), u.trajectory.n_frames))

    def test_slice_iteration(self, multiverse):
        u = multiverse
        frames = []
        for ts in u.trajectory[4:-2:4]:
            frames.append(ts.frame)
        assert_equal(np.array(frames),
                     np.arange(u.trajectory.n_frames)[4:-2:4],
                     err_msg="slicing did not produce the expected frames")


class TestWinDLPolyHistory(TestDLPolyHistory):
    f = WIN_DLP_HISTORY


class TestWinTRJ(TestTRJReader):
    trajectory_file = WIN_TRJ


class TestWinARC(object):
    @staticmethod
    @pytest.fixture
    def WIN_ARC_U():
        return mda.Universe(WIN_ARC)

    def test_n_frames(self, WIN_ARC_U):
        assert len(WIN_ARC_U.trajectory) == 2

    def test_positions(self, WIN_ARC_U):
        assert_almost_equal(WIN_ARC_U.atoms.positions[0],
                            [-6.553398, -1.854369, 0.000000])

    def test_positions_2(self, WIN_ARC_U):
        WIN_ARC_U.trajectory[1]

        assert_almost_equal(WIN_ARC_U.atoms.positions[0],
                            [-0.231579, -0.350841, -0.037475])
