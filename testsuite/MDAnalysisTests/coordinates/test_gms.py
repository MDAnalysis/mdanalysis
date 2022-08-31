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
import pytest

import numpy as np
from numpy.testing import (assert_equal, assert_almost_equal)

import MDAnalysis as mda
from MDAnalysis.coordinates.GMS import GMSReader
from MDAnalysisTests.datafiles import (GMS_ASYMOPT, GMS_ASYMSURF, GMS_SYMOPT)


class _GMSBase(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(self.filename)

    def test_n_frames(self, u):
        assert_equal(u.trajectory.n_frames,
                     self.n_frames,
                     err_msg="Wrong number of frames read from {}".format(
                         self.flavour))

    def test_random_access(self, u):
        u = u
        pos1 = u.atoms[-1].position

        u.trajectory.next()
        u.trajectory.next()

        pos3 = u.atoms[-1].position

        u.trajectory[0]
        assert_equal(u.atoms[-1].position, pos1)

        u.trajectory[2]
        assert_equal(u.atoms[-1].position, pos3)

    @staticmethod
    def _calcFD(u):
        u.trajectory.rewind()
        pp = (u.trajectory.ts._pos[0] - u.trajectory.ts._pos[3])
        z1 = np.sqrt(sum(pp ** 2))
        for i in range(5):
            u.trajectory.next()
        pp = (u.trajectory.ts._pos[0] - u.trajectory.ts._pos[3])
        z2 = np.sqrt(sum(pp ** 2))
        return z1 - z2

    def test_rewind(self, u):
        u.trajectory.rewind()
        assert_equal(u.trajectory.ts.frame, 0, "rewinding to frame 0")

    def test_next(self, u):
        u.trajectory.rewind()
        u.trajectory.next()
        assert_equal(u.trajectory.ts.frame, 1, "loading frame 1")

    def test_dt(self, u):
        assert_almost_equal(u.trajectory.dt,
                            1.0,
                            4,
                            err_msg="wrong timestep dt")

    def test_step5distances(self, u):
        assert_almost_equal(self._calcFD(u), self.step5d, decimal=5,
                            err_msg="Wrong 1-4 atom distance change after "
                                    "5 steps for {}".format(self.flavour))


class TestGMSReader(_GMSBase):
    n_frames = 21
    flavour = "GAMESS C1 optimization"
    step5d = -0.0484664
    filename = GMS_ASYMOPT


class TestGMSReaderSO(_GMSBase):
    n_frames = 8
    flavour = "GAMESS D4H optimization"
    step5d = 0.227637
    filename = GMS_SYMOPT


class TestGMSReaderASS(_GMSBase):
    n_frames = 10
    flavour = "GAMESS C1 surface"
    step5d = -0.499996
    filename = GMS_ASYMSURF
