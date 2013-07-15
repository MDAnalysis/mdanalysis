# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
:Author:   Elizabeth Denning
:Contact:  denniej0@gmail.com

Sample code to use the routine for nucleic acid analysis
For the example provided below, the backbone dihedrals and WC distances

"""

import numpy
import MDAnalysis

from MDAnalysis.analysis import nuclinfo
from MDAnalysis.tests.datafiles import NUCL

from numpy.testing import *
del test
from nose.plugins.attrib import attr

class TestNuclinfo(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(NUCL)
        self.prec = 4
    def tearDown(self):
        del self.universe
    def test_wc_pair(self):
        wc = numpy.array(nuclinfo.wc_pair(self.universe,4,20),dtype=numpy.float32)
        assert_almost_equal(wc, 2.9810174, err_msg="Watson-Crick distance does not match expected value.")
    def test_major_pair(self):
        maj = numpy.array(nuclinfo.major_pair(self.universe,4,20),dtype=numpy.float32)
        assert_almost_equal(maj, 2.9400151, err_msg="Watson-Crick distance does not match expected value.")
    def test_minor_pair(self):
        minor = numpy.array(nuclinfo.minor_pair(self.universe,4,20),dtype=numpy.float32)
        assert_almost_equal(minor, 3.7739358, err_msg="Watson-Crick distance does not match expected value.")
    def test_torsions(self):
        nucl_acid = numpy.array(nuclinfo.tors(self.universe,"RNAA",4),dtype=numpy.float32)
        expected_nucl_acid = numpy.array([296.45596313,  177.79353333,   48.67910767,   81.81109619,  205.58882141,  286.37353516,  198.09187317], dtype=numpy.float32)
        assert_almost_equal(nucl_acid, expected_nucl_acid, self.prec, err_msg="Backbone torsion does not have expected values for alpha, beta, gamma, epsilon, zeta, chi.")


