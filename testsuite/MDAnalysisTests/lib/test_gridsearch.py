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
from numpy.testing import assert_equal, assert_allclose
import numpy as np

import MDAnalysis as mda
from MDAnalysis.lib import grid
from MDAnalysis.lib.pkdtree import PeriodicKDTree

from MDAnalysisTests.datafiles import GRO


@pytest.fixture
def universe():
    u = mda.Universe(GRO)
    return u


@pytest.fixture
def grid_results():
    u = mda.Universe(GRO)
    cutoff = 2
    ref_pos = u.atoms.positions[13937]
    return run_grid_search(u, ref_pos, cutoff)


def run_grid_search(u, ref_pos, cutoff):
    coords = u.atoms.positions

    # Run grid search
    searcher = grid.FastNS(u)
    searcher.set_cutoff(cutoff)
    searcher.set_coords(coords)
    searcher.prepare()

    return searcher.search(ref_pos)


def run_search(universe, ref_id):
    cutoff = 3
    coords = universe.atoms.positions
    ref_pos = coords[ref_id]

    # Run pkdtree search
    pkdt = PeriodicKDTree(universe.atoms.dimensions, bucket_size=10)
    pkdt.set_coords(coords)
    pkdt.search(ref_pos, cutoff)

    results_pkdtree = pkdt.get_indices()
    results_pkdtree.remove(ref_id)
    results_pkdtree = np.array(results_pkdtree)
    results_pkdtree.sort()

    # Run grid search
    results_grid = run_grid_search(universe, ref_pos, cutoff)
    results_grid = results_grid.get_indices()[0]
    results_grid.sort()

    return results_pkdtree, results_grid


def test_gridsearch(universe):
    """Check that pkdtree and grid search return the same results (No PBC needed)"""

    ref_id = 0
    results_pkdtree, results_grid = run_search(universe, ref_id)
    assert_equal(results_pkdtree, results_grid)


def test_gridsearch_PBC(universe):
    """Check that pkdtree and grid search return the same results (PBC needed)"""

    ref_id = 13937
    results_pkdtree, results_grid = run_search(universe, ref_id)
    assert_equal(results_pkdtree, results_grid)


def test_gridsearch_arraycoord(universe):
    """Check the NS routine accepts a single bead coordinate as well as array of coordinates"""
    cutoff = 2
    ref_pos = universe.atoms.positions[:5]

    results = [
        np.array([2, 1, 4, 3]),
        np.array([2, 0, 3]),
        np.array([0, 1, 3]),
        np.array([    2,     0,     1, 38341]),
        np.array([ 6,  0,  5, 17])
    ]

    results_grid = run_grid_search(universe, ref_pos, cutoff).get_indices()

    assert_equal(results_grid, results)


def test_gridsearch_search_coordinates(grid_results):
    """Check the NS routine can return coordinates instead of ids"""

    results = np.array(
        [
            [40.32, 34.25, 55.9],
            [0.61, 76.33, -0.56],
            [0.48999998, 75.9, 0.19999999],
            [-0.11, 76.19, 0.77]
        ])

    assert_allclose(grid_results.get_coordinates()[0], results)


def test_gridsearch_search_distances(grid_results):
    """Check the NS routine can return PBC distances from neighbors"""
    results = np.array([0.096, 0.096, 0.015, 0.179]) * 10  # These distances were obtained using gmx distance
    results.sort()

    rounded_results = np.round(grid_results.get_distances()[0], 2)

    assert_allclose(sorted(rounded_results), results)