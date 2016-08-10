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
from __future__ import absolute_import, print_function

import MDAnalysis
from MDAnalysis.analysis.MSD import MSD

from numpy.testing import assert_almost_equal, assert_allclose
import numpy as np

from MDAnalysisTests.datafiles import waterPSF, waterDCD

class TestMSD(object):
    def setUp(self):
        self.u = MDAnalysis.Universe(waterPSF, waterDCD)
        self.this_MSD = MSD(self.u,['type OT'],0,9,5,write_output=False)

    def tearDown(self):
        del self.u
        del self.this_MSD

    def test_MSD(self):
        answer = [array([[ 0.        ,  0.01054136,  0.04591889,  0.08961838,  0.10439679,
                  0.12196799,  0.11978941,  0.1291698 ,  0.14870591,  0.19049508]]), 
                  array([[10, 10, 10, 10, 10,  5,  5,  5,  5,  5]])]
        
        msd_out = self.this_MSD.run()
        print(msd_out)
	    assert_allclose(answer,msd_out,10**-4,10**-4)
