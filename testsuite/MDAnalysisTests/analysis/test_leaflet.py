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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
import MDAnalysis
import pytest

from numpy.testing import assert_equal, assert_almost_equal
import numpy as np

from MDAnalysis.analysis.leaflet import LeafletFinder, optimize_cutoff
from MDAnalysisTests.datafiles import Martini_membrane_gro


LIPID_HEAD_STRING = "name PO4"


@pytest.fixture()
def universe():
    return MDAnalysis.Universe(Martini_membrane_gro)


@pytest.fixture()
def lipid_heads(universe):
    return universe.select_atoms(LIPID_HEAD_STRING)


def test_leaflet_finder(universe, lipid_heads):
    lfls = LeafletFinder(universe, lipid_heads, pbc=True)
    top_heads, bottom_heads = lfls.groups()
    # Make top be... on top.
    if top_heads.center_of_geometry()[2] < bottom_heads.center_of_geometry()[2]:
        top_heads,bottom_heads = (bottom_heads,top_heads)
    assert_equal(top_heads.indices, np.arange(1,2150,12),
                 err_msg="Found wrong leaflet lipids")
    assert_equal(bottom_heads.indices, np.arange(2521,4670,12),
                 err_msg="Found wrong leaflet lipids")


def test_string_vs_atomgroup_proper(universe, lipid_heads):
    lfls_ag = LeafletFinder(universe, lipid_heads, pbc=True)
    lfls_string = LeafletFinder(universe, LIPID_HEAD_STRING, pbc=True)
    groups_ag = lfls_ag.groups()
    groups_string = lfls_string.groups()
    assert_equal(groups_string[0].indices, groups_ag[0].indices)
    assert_equal(groups_string[1].indices, groups_ag[1].indices)


def test_optimize_cutoff(universe, lipid_heads):
    cutoff, N = optimize_cutoff(universe, lipid_heads, pbc=True)
    assert N == 2
    assert_almost_equal(cutoff, 10.5, decimal=4)
