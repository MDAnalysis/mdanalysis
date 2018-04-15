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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import division, absolute_import

from six.moves import zip

import numpy as np

import pytest

from numpy.testing import (assert_equal, assert_almost_equal, assert_array_almost_equal)

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (PDB, PSF, CRD, DCD,
                                       GRO, XTC, TRR, PDB_small, PDB_closed)


class TestChainReader(object):
    prec = 3

    @pytest.fixture()
    def universe(self):
        return mda.Universe(PSF,
                            [DCD, CRD, DCD, CRD, DCD, CRD, CRD])

    def test_next_trajectory(self, universe):
        universe.trajectory.rewind()
        universe.trajectory.next()
        assert_equal(universe.trajectory.ts.frame, 1, "loading frame 2")

    def test_n_atoms(self, universe):
        assert_equal(universe.trajectory.n_atoms, 3341,
                     "wrong number of atoms")

    def test_n_frames(self, universe):
        assert_equal(universe.trajectory.n_frames, 3 * 98 + 4,
                     "wrong number of frames in chained dcd")

    def test_iteration(self, universe):
        for ts in universe.trajectory:
            pass  # just forward to last frame
        assert_equal(
            universe.trajectory.n_frames - 1, ts.frame,
            "iteration yielded wrong number of frames ({0:d}), "
            "should be {1:d}".format(ts.frame, universe.trajectory.n_frames))

    def test_jump_lastframe_trajectory(self, universe):
        universe.trajectory[-1]
        assert_equal(universe.trajectory.ts.frame + 1, universe.trajectory.n_frames,
                     "indexing last frame with trajectory[-1]")

    def test_slice_trajectory(self, universe):
        frames = [ts.frame for ts in universe.trajectory[5:17:3]]
        assert_equal(frames, [5, 8, 11, 14], "slicing dcd [5:17:3]")

    def test_full_slice(self, universe):
        trj_iter = universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(universe.trajectory.n_frames))

    def test_frame_numbering(self, universe):
        universe.trajectory[98]  # index is 0-based and frames are 0-based
        assert_equal(universe.trajectory.frame, 98, "wrong frame number")

    def test_frame(self, universe):
        universe.trajectory[0]
        coord0 = universe.atoms.positions.copy()
        # forward to frame where we repeat original dcd again:
        # dcd:0..97 crd:98 dcd:99..196
        universe.trajectory[99]
        assert_equal(
            universe.atoms.positions, coord0,
            "coordinates at frame 1 and 100 should be the same!")

    def test_time(self, universe):
        universe.trajectory[30]  # index and frames 0-based
        assert_almost_equal(universe.trajectory.time,
                            30.0,
                            5,
                            err_msg="Wrong time of frame")

    def test_write_dcd(self, universe, tmpdir):
        """test that ChainReader written dcd (containing crds) is correct
        (Issue 81)"""
        outfile = str(tmpdir) + "chain-reader.dcd"
        with mda.Writer(outfile, universe.atoms.n_atoms) as W:
            for ts in universe.trajectory:
                W.write(universe)
        universe.trajectory.rewind()
        u = mda.Universe(PSF, outfile)
        for (ts_orig, ts_new) in zip(universe.trajectory,
                                     u.trajectory):
            assert_almost_equal(
                ts_orig._pos,
                ts_new._pos,
                self.prec,
                err_msg="Coordinates disagree at frame {0:d}".format(
                    ts_orig.frame))       

class TestChainReaderCommonDt(object):
    common_dt = 100.0
    prec = 3

    @pytest.fixture()
    def trajectory(self):
        universe = mda.Universe(PSF,
                                [DCD, CRD, DCD, CRD, DCD, CRD, CRD],
                                dt=self.common_dt)
        return universe.trajectory

    def test_time(self, trajectory):
        # We test this for the beginning, middle and end of the trajectory.
        for frame_n in (0, trajectory.n_frames // 2, -1):
            trajectory[frame_n]
            assert_almost_equal(trajectory.time,
                                trajectory.frame * self.common_dt,
                                5,
                                err_msg="Wrong time for frame {0:d}".format(
                                    frame_n))


class TestChainReaderFormats(object):
    """Test of ChainReader with explicit formats (Issue 76)."""

    def test_set_all_format_tuples(self):
        universe = mda.Universe(GRO, [(PDB, 'pdb'), (XTC, 'xtc'),
                                      (TRR, 'trr')])
        assert_equal(universe.trajectory.n_frames, 21)

    def test_set_one_format_tuple(self):
        universe = mda.Universe(PSF, [(PDB_small, 'pdb'), DCD])
        assert_equal(universe.trajectory.n_frames, 99)

    def test_set_all_formats(self):
        universe = mda.Universe(PSF, [PDB_small, PDB_closed], format='pdb')
        assert_equal(universe.trajectory.n_frames, 2)
