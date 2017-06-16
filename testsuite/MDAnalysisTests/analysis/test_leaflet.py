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
from __future__ import  absolute_import
import MDAnalysis
from MDAnalysisTests import module_not_found

from numpy.testing import TestCase, assert_equal, assert_almost_equal, dec
import numpy as np

from MDAnalysisTests.datafiles import Martini_membrane_gro

class TestLeafletFinder(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(Martini_membrane_gro, Martini_membrane_gro)
        self.lipid_heads = self.universe.select_atoms("name PO4")
        self.lipid_head_string = "name PO4"

    def tearDown(self):
        del self.universe
        del self.lipid_heads
        del self.lipid_head_string

    def test_leaflet_finder(self):
        from MDAnalysis.analysis.leaflet import LeafletFinder
        lfls = LeafletFinder(self.universe, self.lipid_heads, pbc=True)
        top_heads, bottom_heads = lfls.groups()
        # Make top be... on top.
        if top_heads.center_of_geometry()[2] < bottom_heads.center_of_geometry()[2]:
            top_heads,bottom_heads = (bottom_heads,top_heads)
        assert_equal(top_heads.indices, np.arange(1,2150,12), err_msg="Found wrong leaflet lipids")
        assert_equal(bottom_heads.indices, np.arange(2521,4670,12), err_msg="Found wrong leaflet lipids")


    def test_string_vs_atomgroup_proper(self):
        from MDAnalysis.analysis.leaflet import LeafletFinder
        lfls_ag = LeafletFinder(self.universe, self.lipid_heads, pbc=True)
        lfls_string = LeafletFinder(self.universe, self.lipid_head_string, pbc=True)
        groups_ag = lfls_ag.groups()
        groups_string = lfls_string.groups()
        assert_equal(groups_string[0].indices, groups_ag[0].indices)
        assert_equal(groups_string[1].indices, groups_ag[1].indices)

    def test_optimize_cutoff(self):
        from MDAnalysis.analysis.leaflet import optimize_cutoff
        cutoff, N = optimize_cutoff(self.universe, self.lipid_heads, pbc=True)
        assert_equal(N, 2)
        assert_almost_equal(cutoff, 10.5, decimal=4)

