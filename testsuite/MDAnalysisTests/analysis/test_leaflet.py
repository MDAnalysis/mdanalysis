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
import MDAnalysis
from MDAnalysisTests import module_not_found

from numpy.testing import TestCase, assert_equal, dec
import numpy as np

from MDAnalysisTests.datafiles import Martini_membrane_gro

class TestLeafletFinder(TestCase):
    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def setUp(self):
        self.universe = MDAnalysis.Universe(Martini_membrane_gro, Martini_membrane_gro)
        self.lipid_heads = self.universe.select_atoms("name PO4")

    def tearDown(self):
        del self.universe

    def test_leaflet_finder(self):
        from MDAnalysis.analysis.leaflet import LeafletFinder
        lfls = LeafletFinder(self.universe, self.lipid_heads, pbc=True)
        top_heads, bottom_heads = lfls.groups()
        # Make top be... on top.
        if top_heads.center_of_geometry()[2] < bottom_heads.center_of_geometry()[2]:
            top_heads,bottom_heads = (bottom_heads,top_heads)
        assert_equal(top_heads.indices, np.arange(1,2150,12), err_msg="Found wrong leaflet lipids")
        assert_equal(bottom_heads.indices, np.arange(2521,4670,12), err_msg="Found wrong leaflet lipids")


