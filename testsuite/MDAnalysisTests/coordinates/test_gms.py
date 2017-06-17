# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
from __future__ import absolute_import
from six.moves import range

import MDAnalysis as mda
import numpy as np

from numpy.testing import (assert_equal, assert_almost_equal)

from MDAnalysisTests.datafiles import (GMS_ASYMOPT, GMS_ASYMSURF, GMS_SYMOPT)


class _GMSBase(object):
    def tearDown(self):
        del self.u
        del self.n_frames
        del self.flavour
        del self.step5d

    def test_n_frames(self):
        assert_equal(self.u.trajectory.n_frames,
                     self.n_frames,
                     err_msg="Wrong number of frames read from {}".format(self.flavour))

    def test_random_access(self):
        u = self.u
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

    def test_rewind(self):
        self.u.trajectory.rewind()
        assert_equal(self.u.trajectory.ts.frame, 0, "rewinding to frame 0")

    def test_next(self):
        self.u.trajectory.rewind()
        self.u.trajectory.next()
        assert_equal(self.u.trajectory.ts.frame, 1, "loading frame 1")

    def test_dt(self):
        assert_almost_equal(self.u.trajectory.dt,
                            1.0,
                            4,
                            err_msg="wrong timestep dt")

    def test_step5distances(self):
        assert_almost_equal(self._calcFD(self.u), self.step5d, decimal=5,
                            err_msg="Wrong 1-4 atom distance change after "
                            "5 steps for {}".format(self.flavour))


class TestGMSReader(_GMSBase):
    def setUp(self):
        self.u =  mda.Universe(GMS_ASYMOPT)
        self.n_frames = 21
        self.flavour = "GAMESS C1 optimization"
        self.step5d = -0.0484664

class TestGMSReaderSO(_GMSBase):
    def setUp(self):
        self.u = mda.Universe(GMS_SYMOPT)
        self.n_frames = 8
        self.flavour = "GAMESS D4H optimization"
        self.step5d = 0.227637

class TestGMSReaderASS(_GMSBase):
    def setUp(self):
        self.u = mda.Universe(GMS_ASYMSURF)
        self.n_frames = 10
        self.flavour = "GAMESS C1 surface"
        self.step5d = -0.499996
