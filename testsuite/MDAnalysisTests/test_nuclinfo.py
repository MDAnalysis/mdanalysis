# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
:Author:   Elizabeth Denning
:Contact:  denniej0@gmail.com

Sample code to use the routine for nucleic acid analysis
For the example provided below, the backbone dihedrals and WC distances

"""

import numpy as np
import MDAnalysis

from MDAnalysis.analysis import nuclinfo
from MDAnalysis.tests.datafiles import NUCL

from numpy.testing import assert_almost_equal ,assert_array_almost_equal, TestCase

from nose.plugins.attrib import attr


class TestNuclinfo(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(NUCL)
        self.prec = 4

    def tearDown(self):
        del self.universe

    def test_wc_pair(self):
        seg1 = self.universe.residues[3].segids[0]
        seg2 = self.universe.residues[19].segids[0]
        wc = nuclinfo.wc_pair(self.universe, 4, 20, seg1, seg2)
        assert_almost_equal(wc, 2.9810174, err_msg="Watson-Crick distance does not match expected value.")

    def test_major_pair(self):
        seg1 = self.universe.residues[3].segids[0]
        seg2 = self.universe.residues[19].segids[0]
        maj = nuclinfo.major_pair(self.universe, 4, 20, seg1, seg2)
        assert_almost_equal(maj, 2.9400151, err_msg="Watson-Crick distance does not match expected value.")

    def test_minor_pair(self):
        seg1 = self.universe.residues[3].segids[0]
        seg2 = self.universe.residues[19].segids[0]

        minor = nuclinfo.minor_pair(self.universe, 4, 20, seg1, seg2)
        assert_almost_equal(minor, 3.7739358, err_msg="Watson-Crick distance does not match expected value.")

    def test_torsions(self):
        nucl_acid = np.array(nuclinfo.tors(self.universe, "RNAA", 4), dtype=np.float32)
        expected_nucl_acid = np.array(
            [296.45596313, 177.79353333, 48.67910767, 81.81109619, 205.58882141, 286.37353516, 198.09187317],
            dtype=np.float32)
        assert_almost_equal(nucl_acid, expected_nucl_acid, self.prec,
                            err_msg="Backbone torsion does not have expected values for "
                                    "alpha, beta, gamma, epsilon, zeta, chi.")

    def test_hydroxyl(self):
        hydroxyls = np.array([nuclinfo.hydroxyl(self.universe,
                                                   self.universe.atoms.segids[0], resid)
                                 for resid in (7, 10, 11, 22)])
        expected_hydroxyls = np.array(
            [ 122.73991394,  123.34986115,  123.20658112,  122.57156372],
            dtype=np.float32)
        assert_array_almost_equal(hydroxyls, expected_hydroxyls, self.prec,
                                  err_msg="RNA hydroxyl dihedrals do not match")

    def test_pseudo_dihe_baseflip(self):
        seg1 = self.universe.residues[3].segids[0]
        seg2 = self.universe.residues[19].segids[0]

        # There is not really a baseflip, just testing the code...
        flip = nuclinfo.pseudo_dihe_baseflip(self.universe, 4, 20, 5, seg1=seg1, seg2=seg2, seg3=seg1)
        assert_almost_equal(flip, 322.0826, self.prec,
                            err_msg="pseudo_dihedral for resid 5 against 4--20 do not match")
