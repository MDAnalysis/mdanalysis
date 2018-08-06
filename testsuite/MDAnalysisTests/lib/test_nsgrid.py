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
from MDAnalysisTests.datafiles import GRO, Martini_membrane_gro
from MDAnalysis.lib import nsgrid


@pytest.fixture
def universe():
    u = mda.Universe(GRO)
    return u

def run_grid_search(u, ref_id, cutoff=3):
    coords = u.atoms.positions
    searchcoords = u.atoms.positions[ref_id]
    if searchcoords.shape == (3, ):
        searchcoords = searchcoords[None, :]
    # Run grid search
    searcher = nsgrid.FastNS(cutoff, coords, box=u.dimensions)

    return searcher.search(searchcoords)


def test_pbc_box():
    """Check that PBC box accepts only well-formated boxes"""
    pbc = True
    with pytest.raises(TypeError):
        nsgrid.PBCBox([])

    with pytest.raises(ValueError):
        nsgrid.PBCBox(np.zeros((3)), pbc)  # Bad shape
        nsgrid.PBCBox(np.zeros((3, 3)), pbc)  # Collapsed box
        nsgrid.PBCBOX(np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1]]), pbc)  # 2D box
        nsgrid.PBCBOX(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), pbc)  # Box provided as array of integers
        nsgrid.PBCBOX(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float), pbc)  # Box provided as array of double


def test_nsgrid_badcutoff(universe):
    with pytest.raises(ValueError):
        run_grid_search(universe, 0, -4)
        run_grid_search(universe, 0, 100000)


def test_ns_grid_noneighbor(universe):
    """Check that grid search returns empty lists/arrays when there is no neighbors"""
    ref_id = 0
    cutoff = 0.5

    results_grid = run_grid_search(universe, ref_id, cutoff)

    # same indices will be selected as neighbour here
    assert len(results_grid.get_distances()[0]) == 1
    assert len(results_grid.get_indices()[0]) == 1
    assert len(results_grid.get_pairs()) == 1
    assert len(results_grid.get_pair_distances()) == 1


def test_nsgrid_PBC_rect():
    """Check that nsgrid works with rect boxes and PBC"""
    ref_id = 191
    results = np.array([191, 192, 672, 682, 683, 684, 995, 996, 2060, 2808, 3300, 3791,
                        3792]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!

    universe = mda.Universe(Martini_membrane_gro)
    cutoff = 7

    # FastNS is called differently to max coverage
    searcher = nsgrid.FastNS(cutoff, universe.atoms.positions, box=universe.dimensions)

    results_grid = searcher.search(universe.atoms.positions[ref_id][None, :]).get_indices()[0]

    results_grid2 = searcher.search(universe.atoms.positions).get_indices()  # call without specifying any ids, should do NS for all beads

    assert_equal(np.sort(results), np.sort(results_grid))
    assert_equal(len(universe.atoms), len(results_grid2))
    assert searcher.cutoff == 7
    assert_equal(np.sort(results_grid), np.sort(results_grid2[ref_id]))


def test_nsgrid_PBC(universe):
    """Check that grid search works when PBC is needed"""

    ref_id = 13937
    results = np.array([4398, 4401, 13938, 13939, 13940, 13941, 17987, 23518, 23519, 23521, 23734,
                        47451]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!

    results_grid = run_grid_search(universe, ref_id).get_indices()[0]

    assert_equal(np.sort(results), np.sort(results_grid))


def test_nsgrid_pairs(universe):
    """Check that grid search returns the proper pairs"""

    ref_id = 13937
    neighbors = np.array([4398, 4401, 13938, 13939, 13940, 13941, 17987, 23518, 23519, 23521, 23734,
                          47451]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!
    results = []

    results = np.array(results)

    results_grid = run_grid_search(universe, ref_id).get_pairs()

    assert_equal(np.sort(neighbors, axis=0), np.sort(results_grid[:, 1], axis=0))


def test_nsgrid_pair_distances(universe):
    """Check that grid search returns the proper pair distances"""

    ref_id = 13937
    results = np.array([0.0, 0.270, 0.285, 0.096, 0.096, 0.015, 0.278, 0.268, 0.179, 0.259, 0.290,
                        0.270]) * 10  # These distances where obtained by gmx distance so they are in nm

    results_grid = run_grid_search(universe, ref_id).get_pair_distances()

    assert_allclose(np.sort(results), np.sort(results_grid), atol=1e-2)


def test_nsgrid_distances(universe):
    """Check that grid search returns the proper distances"""

    ref_id = 13937
    results = np.array([0.0, 0.270, 0.285, 0.096, 0.096, 0.015, 0.278, 0.268, 0.179, 0.259, 0.290,
                        0.270]) * 10  # These distances where obtained by gmx distance so they are in nm

    results_grid = run_grid_search(universe, ref_id).get_distances()[0]

    assert_allclose(np.sort(results), np.sort(results_grid), atol=1e-2)


@pytest.mark.parametrize('box, results',
                         ((None, [3, 13, 24]),
                          (np.array([10., 10., 10., 90., 90., 90.]), [3, 13, 24, 39, 67]),
                          (np.array([10., 10., 10., 60., 75., 90.]), [3, 13, 24, 39, 60, 79])))
def test_nsgrid_search(box, results):
    np.random.seed(90003)
    points = (np.random.uniform(low=0, high=1.0,
                        size=(100, 3))*(10.)).astype(np.float32)
    cutoff = 2.0
    query = np.array([1., 1., 1.], dtype=np.float32).reshape((1, 3))

    if box is None:
        pseudobox = np.zeros(6, dtype=np.float32)
        all_coords = np.concatenate([points, query])
        lmax = all_coords.max(axis=0)
        lmin = all_coords.min(axis=0)
        pseudobox[:3] = 1.1*(lmax - lmin)
        pseudobox[3:] = 90.
        shiftpoints, shiftquery = points.copy(), query.copy()
        shiftpoints -= lmin
        shiftquery -= lmin
        searcher = nsgrid.FastNS(cutoff, shiftpoints, box=pseudobox, pbc=False)
        searchresults = searcher.search(shiftquery)
    else:
        searcher = nsgrid.FastNS(cutoff, points, box)
        searchresults = searcher.search(query)
    indices = searchresults.get_indices()[0]
    assert_equal(np.sort(indices), results)


@pytest.mark.parametrize('box, result',
                         ((None, 21),
                          (np.array([0., 0., 0., 90., 90., 90.]), 21),
                          (np.array([10., 10., 10., 90., 90., 90.]), 26),
                          (np.array([10., 10., 10., 60., 75., 90.]), 33)))
def test_nsgrid_selfsearch(box, result):
    np.random.seed(90003)
    points = (np.random.uniform(low=0, high=1.0,
                        size=(100, 3))*(10.)).astype(np.float32)
    cutoff = 1.0
    if box is None or np.allclose(box[:3], 0):
        # create a pseudobox
        # define the max range
        # and supply the pseudobox
        # along with only one set of coordinates
        pseudobox = np.zeros(6, dtype=np.float32)
        lmax = points.max(axis=0)
        lmin = points.min(axis=0)
        pseudobox[:3] = 1.1*(lmax - lmin)
        pseudobox[3:] = 90.
        shiftref = points.copy()
        shiftref -= lmin
        searcher = nsgrid.FastNS(cutoff, shiftref, box=pseudobox, pbc=False)
        searchresults = searcher.self_search()
    else:
        searcher = nsgrid.FastNS(cutoff, points, box=box)
        searchresults = searcher.self_search()
    pairs = searchresults.get_pairs()
    assert_equal(len(pairs)//2, result)
