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
import MDAnalysis as mda
import numpy as np
import pytest

from numpy.testing import assert_equal

from MDAnalysisTests.datafiles import (DMS)


class TestDMSReader(object):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(DMS)

    @pytest.fixture()
    def ts(self, universe):
        return universe.trajectory.ts

    def test_global_cell(self, ts):
        assert ts.dimensions is None

    # cythonised class can no longer raise AttributeError
    # so changed to test of has_velocities 
    def test_velocities(self, ts):
        assert_equal(ts.has_velocities, False)

    def test_number_of_coords(self, universe):
        # Desired value taken from VMD
        #      Info)    Atoms: 3341
        assert_equal(len(universe.atoms), 3341)

    def test_coords_atom_0(self, universe):
        # Desired coordinates taken directly from the SQLite file. Check unit
        # conversion
        coords_0 = np.array([-11.0530004501343,
                             26.6800003051758,
                             12.7419996261597, ],
                            dtype=np.float32)
        assert_equal(universe.atoms[0].position, coords_0)

    def test_n_frames(self, universe):
        assert_equal(universe.trajectory.n_frames, 1,
                     "wrong number of frames in pdb")

    def test_time(self, universe):
        assert_equal(universe.trajectory.time, 0.0,
                     "wrong time of the frame")

    def test_frame(self, universe):
        assert_equal(universe.trajectory.frame, 0, "wrong frame number "
                                                   "(0-based, should be 0 for single frame readers)")

    def test_frame_index_0(self, universe):
        universe.trajectory[0]
        assert_equal(universe.trajectory.ts.frame, 0,
                     "frame number for frame index 0 should be 0")

    def test_frame_index_1_raises_IndexError(self, universe):
        with pytest.raises(IndexError):
            universe.trajectory[1]
