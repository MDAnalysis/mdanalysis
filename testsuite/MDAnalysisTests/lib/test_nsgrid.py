# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2018 The MDAnalysis Development Team and contributors
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

from __future__ import print_function, absolute_import

import pytest

from numpy.testing import assert_equal
import numpy as np

import MDAnalysis as mda
from MDAnalysisTests.datafiles import GRO
from MDAnalysis.lib import nsgrid


@pytest.fixture
def universe():
    u = mda.Universe(GRO)
    return u



def run_grid_search(u, ids, cutoff=3):
    coords = u.atoms.positions

    # Run grid search
    searcher = nsgrid.FastNS(u, cutoff, coords, debug=True)


    return searcher.search(ids, debug=False)


def test_gridsearch(universe):
    """Check that pkdtree and grid search return the same results (No PBC needed)"""

    ref_id = 0

    cutoff = 3
    results = np.array([2, 3, 4, 5, 6, 7, 8, 9, 18, 19, 1211, 10862, 10865, 17582, 17585, 38342,
                        38345]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!

    results_grid = run_grid_search(universe, ref_id, cutoff).get_indices()[0]

    assert_equal(results, results_grid)


def test_gridsearch_PBC(universe):
    """Check that pkdtree and grid search return the same results (No PBC needed)"""

    ref_id = 13937
    results = np.array([4398, 4401, 13939, 13940, 13941, 17987, 23518, 23519, 23521, 23734,
                        47451]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!

    results_grid = run_grid_search(universe, ref_id).get_indices()[0]

    assert_equal(results, results_grid)



# def test_gridsearch_PBC(universe):
#     """Check that pkdtree and grid search return the same results (PBC needed)"""
#
#     ref_id = 13937
#     results_pkdtree, results_grid = run_search(universe, ref_id)
#     assert_equal(results_pkdtree, results_grid)
#
#
# def test_gridsearch_arraycoord(universe):
#     """Check the NS routine accepts a single bead coordinate as well as array of coordinates"""
#     cutoff = 2
#     ref_pos = universe.atoms.positions[:5]
#
#     results = [
#         np.array([2, 1, 4, 3]),
#         np.array([2, 0, 3]),
#         np.array([0, 1, 3]),
#         np.array([    2,     0,     1, 38341]),
#         np.array([ 6,  0,  5, 17])
#     ]
#
#     results_grid = run_grid_search(universe, ref_pos, cutoff).get_indices()
#
#     assert_equal(results_grid, results)
#
#
# def test_gridsearch_search_coordinates(grid_results):
#     """Check the NS routine can return coordinates instead of ids"""
#
#     results = np.array(
#         [
#             [40.32, 34.25, 55.9],
#             [0.61, 76.33, -0.56],
#             [0.48999998, 75.9, 0.19999999],
#             [-0.11, 76.19, 0.77]
#         ])
#
#     assert_allclose(grid_results.get_coordinates()[0], results)
#
#
# def test_gridsearch_search_distances(grid_results):
#     """Check the NS routine can return PBC distances from neighbors"""
#     results = np.array([0.096, 0.096, 0.015, 0.179]) * 10  # These distances were obtained using gmx distance
#     results.sort()
#
#     rounded_results = np.round(grid_results.get_distances()[0], 2)
#
#     assert_allclose(sorted(rounded_results), results)
