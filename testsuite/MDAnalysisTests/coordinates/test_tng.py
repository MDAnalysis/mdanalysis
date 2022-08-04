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
import numpy as np
import pytest

from numpy.testing import (
    assert_equal,
    assert_almost_equal
)

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (TNG_traj, TNG_traj_gro)


class TestTNGReaderTest(object):
    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(TNG_traj_gro, TNG_traj)

    def test_n_atoms(self, universe):
        assert_equal(universe.trajectory.n_atoms, 1000)

    def test_n_frames(self, universe):
        assert_equal(universe.trajectory.n_frames, 101,
                     "wrong number of frames in xyz")

    def test_initial_frame_is_0(self, universe):
        assert_equal(universe.trajectory.ts.frame, 0,
                     "initial frame is not 0 but {0}".format(
                         universe.trajectory.ts.frame))

    def test_starts_with_first_frame(self, universe):
        """Test that coordinate arrays are filled as soon as the trajectory
        has been opened."""
        assert np.any(universe.atoms.positions > 0), "Reader does not " \
                                                     "populate positions right away."

    def test_rewind(self, universe):
        trj = universe.trajectory
        trj.next()
        trj.next()  # for readers that do not support indexing
        assert_equal(trj.ts.frame, 2,
                     "failed to forward to frame 2 (frameindex 2)")
        trj.rewind()
        assert_equal(trj.ts.frame, 0, "failed to rewind to first frame")
        assert np.any(universe.atoms.positions > 0), "Reader does not " \
                                                     "populate positions after rewinding."

    def test_full_slice(self, universe):
        trj_iter = universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(universe.trajectory.n_frames))

    def test_random_access(self, universe):
        pos1 = universe.atoms[0].position
        universe.trajectory.next()
        universe.trajectory.next()
        pos3 = universe.atoms[0].position

        universe.trajectory[0]

        assert_equal(universe.atoms[0].position, pos1)

        universe.trajectory[2]

        assert_equal(universe.atoms[0].position, pos3)