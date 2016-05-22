# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function

import MDAnalysis
import MDAnalysis.analysis.gnm

from numpy.testing import (TestCase, assert_equal, assert_almost_equal)
import numpy as np

from nose.plugins.attrib import attr

from MDAnalysisTests.datafiles import GRO, XTC
from MDAnalysisTests import tempdir

class TestGNM(TestCase):
    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.universe = MDAnalysis.Universe(GRO, XTC)

    def tearDown(self):
        del self.universe
        del self.tmpdir

    def test_gnm(self):
        gnm = MDAnalysis.analysis.gnm.GNMAnalysis(self.universe, ReportVector="output.txt")
        gnm.run()
        result = gnm.results
        assert_equal(len(result), 10)
        time, eigenvalues, eigenvectors = zip(*result)
        assert_almost_equal(time, range(0, 1000, 100), decimal=4)
        assert_almost_equal(eigenvalues,
          [ 2.0287113e-15, 4.1471575e-15, 1.8539533e-15, 4.3810359e-15,
            3.9607304e-15, 4.1289113e-15, 2.5501084e-15, 4.0498182e-15,
            4.2058769e-15, 3.9839431e-15])

    def test_gnm_run_skip(self):
        gnm = MDAnalysis.analysis.gnm.GNMAnalysis(self.universe)
        gnm.run(skip=3)
        result = gnm.results
        assert_equal(len(result), 4)
        time, eigenvalues, eigenvectors = zip(*result)
        assert_almost_equal(time, range(0, 1200, 300), decimal=4)
        assert_almost_equal(eigenvalues,
          [ 2.0287113e-15, 4.3810359e-15, 2.5501084e-15, 3.9839431e-15])

    def test_generate_kirchoff(self):
        gnm = MDAnalysis.analysis.gnm.GNMAnalysis(self.universe)
        gnm.run()
        gen = gnm.generate_kirchoff()
        assert_almost_equal(gen[0],
          [ 7,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    @attr('slow')
    def test_closeContactGNMAnalysis(self):
        gnm = MDAnalysis.analysis.gnm.closeContactGNMAnalysis(self.universe)
        gnm.run()

        result = gnm.results
        assert_equal(len(result), 10)
        time, eigenvalues, eigenvectors = zip(*result)
        assert_almost_equal(time, range(0, 1000, 100), decimal=4)
        assert_almost_equal(eigenvalues,
          [ 0.1502614,  0.1426407,  0.1412389,  0.1478305,  0.1425449,
            0.1563304,  0.156915 ,  0.1503619,  0.1572592,  0.1542063])

        gen = gnm.generate_kirchoff()
        assert_almost_equal(gen[0],
          [ 16.326744128018923, -2.716098853586913, -1.94736842105263, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, -0.05263157894736842, 0.0, 0.0, 0.0, -3.3541953679557905, 0.0, -1.4210526315789465, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -1.0423368771244421, -1.3006649542861801, -0.30779350562554625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.927172649945531, -0.7509392614826383,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, -2.263157894736841, -0.24333213169614382])

    @attr('slow')
    def test_closeContactGNMAnalysis_noMassWeight(self):
        gnm = MDAnalysis.analysis.gnm.closeContactGNMAnalysis(self.universe, MassWeight=False)
        gnm.run()

        result = gnm.results
        assert_equal(len(result), 10)
        time, eigenvalues, eigenvectors = zip(*result)
        assert_almost_equal(time, range(0, 1000, 100), decimal=4)
        assert_almost_equal(eigenvalues,
          [ 2.4328739,  2.2967251,  2.2950061,  2.4110916,  2.3271343,
            2.5213111,  2.5189955,  2.4481649,  2.5224835,  2.4824345])

        gen = gnm.generate_kirchoff()
        assert_almost_equal(gen[0],
          [ 303.0, -58.0, -37.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
            0.0, 0.0, 0.0, -67.0, 0.0, -27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -17.0, -15.0,
            -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -14.0, -15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -43.0, -3.0])
