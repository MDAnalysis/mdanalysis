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
from __future__ import division, absolute_import
import numpy as np
from numpy.testing import (
    assert_,
    assert_array_almost_equal,
)

from MDAnalysisTests.datafiles import MMTF, MMTF_gz

from MDAnalysis.coordinates.MMTF import MMTFReader


class TestMMTFReader(object):
    def setUp(self):
        self.r = MMTFReader(MMTF)

    def tearDown(self):
        del self.r

    def test_read_frame_size(self):
        assert_(self.r.ts.n_atoms == 512)

    def test_read_positions(self):
        assert_array_almost_equal(self.r.ts.positions[0],
                                  np.array([-0.798, 12.632, 23.231]),
                                  decimal=4)
        assert_array_almost_equal(self.r.ts.positions[-1],
                                  np.array([10.677, 15.517, 11.1]),
                                  decimal=4)

    def test_velocities(self):
        assert_(not self.r.ts.has_velocities)

    def test_forces(self):
        assert_(not self.r.ts.has_forces)

    def test_len(self):
        # should be single frame
        assert_(len(self.r) == 1)

class TestMMTFReaderGZ(object):
    def setUp(self):
        self.r = MMTFReader(MMTF_gz)

    def tearDown(self):
        del self.r

    def test_read_frame_size(self):
        assert_(self.r.ts.n_atoms == 1140)

    def test_read_positions(self):
        assert_array_almost_equal(self.r.ts.positions[0],
                                  np.array([38.428, 16.440, 28.841]),
                                  decimal=4)
        assert_array_almost_equal(self.r.ts.positions[-1],
                                  np.array([36.684, 27.024, 20.468]),
                                  decimal=4)

    def test_velocities(self):
        assert_(not self.r.ts.has_velocities)

    def test_forces(self):
        assert_(not self.r.ts.has_forces)

    def test_len(self):
        # should be single frame
        assert_(len(self.r) == 1)
