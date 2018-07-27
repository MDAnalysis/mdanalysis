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



def run_grid_search(u, ids, cutoff=3):
    coords = u.atoms.positions

    # Run grid search
    searcher = nsgrid.FastNS(u.dimensions, cutoff, coords, debug=True)

    return searcher.search(ids)


def test_pbc_badbox():
    """Check that PBC box accepts only well-formated boxes"""
    with pytest.raises(TypeError):
        nsgrid.PBCBox([])

    with pytest.raises(ValueError):
        nsgrid.PBCBox(np.zeros((3)))  # Bad shape
        nsgrid.PBCBox(np.zeros((3, 3)))  # Collapsed box
        nsgrid.PBCBOX(np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1]]))  # 2D box
        nsgrid.PBCBOX(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))  # Box provided as array of integers
        nsgrid.PBCBOX(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float))  # Box provided as array of double


def test_pbc_distances():
    """Check that PBC box computes distances"""
    box = np.identity(3, dtype=np.float32)
    bad = np.array([0.1, 0.2], dtype=np.float32)
    a = np.array([0.1, 0.1, 0.1], dtype=np.float32)
    b = np.array([1.1, -0.1, 0.2], dtype=np.float32)
    dx = np.array([0, -0.2, 0.1], dtype=np.float32)
    pbcbox = nsgrid.PBCBox(box)

    with pytest.raises(ValueError):
        pbcbox.distance(a, bad)
        pbcbox.distance(bad, a)

        pbcbox.distance2(a, bad)
        pbcbox.distance2(bad, a)

        pbcbox.dx(bad, a)
        pbcbox.dx(a, bad)

    assert_equal(pbcbox.dx(a, b), dx)
    assert_allclose(pbcbox.distance(a, b), np.sqrt(np.sum(dx*dx)), atol=1e-5)
    assert_allclose(pbcbox.distance2(a, b), np.sum(dx*dx), atol=1e-5)


def test_pbc_put_in_bbox():
    "Check that PBC put beads in brick-shaped box"
    box = np.identity(3, dtype=np.float32)
    coords = np.array(
        [
            [0.1, 0.1, 0.1],
            [-0.1, 1.1, 0.9]
        ],
        dtype=np.float32
    )
    results = np.array(
        [
            [0.1, 0.1, 0.1],
            [0.9, 0.1, 0.9]
        ],
        dtype=np.float32
    )

    pbcbox = nsgrid.PBCBox(box)

    assert_allclose(pbcbox.put_atoms_in_bbox(coords), results, atol=1e-5)


def test_nsgrid_badinit():
    with pytest.raises(TypeError):
        nsgrid.FastNS(None, 1)

def test_nsgrid_badcutoff(universe):
    with pytest.raises(ValueError):
        run_grid_search(universe, 0, -4)
        run_grid_search(universe, 0, 100000)

def test_ns_grid_noneighbor(universe):
    """Check that grid search returns empty lists/arrays when there is no neighbors"""
    ref_id = 0
    cutoff = 0.5

    results_grid = run_grid_search(universe, ref_id, cutoff)

    assert len(results_grid.get_coordinates()[0]) == 0
    assert len(results_grid.get_distances()[0]) == 0
    assert len(results_grid.get_indices()[0]) == 0
    assert len(results_grid.get_pairs()) == 0
    assert len(results_grid.get_pair_distances()) == 0
    assert len(results_grid.get_pair_coordinates()) == 0


def test_nsgrid_noPBC(universe):
    """Check that grid search works when no PBC is needed"""

    ref_id = 0

    cutoff = 3
    results = np.array([2, 3, 4, 5, 6, 7, 8, 9, 18, 19, 1211, 10862, 10865, 17582, 17585, 38342,
                        38345]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!

    results_grid = run_grid_search(universe, ref_id, cutoff).get_indices()[0]

    assert_equal(results, results_grid)


def test_nsgrid_PBC_rect():
    """Check that nsgrid works with rect boxes and PBC"""
    ref_id = 191
    results = np.array([191, 672, 682, 683, 684, 995, 996, 2060, 2808, 3300, 3791,
                        3792]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!

    universe = mda.Universe(Martini_membrane_gro)
    cutoff = 7

    # FastNS is called differently to max coverage
    searcher = nsgrid.FastNS(universe.dimensions, cutoff, universe.atoms.positions, prepare=False)

    results_grid = searcher.search([ref_id,]).get_indices()[0] # pass the id as a list for test+coverage purpose

    searcher.prepare()  # Does nothing, called here for coverage
    results_grid2 = searcher.search().get_indices() # call without specifying any ids, should do NS for all beads

    assert_equal(results, results_grid)
    assert_equal(len(universe.atoms), len(results_grid2))
    assert searcher.cutoff == 7
    assert_equal(results_grid, results_grid2[ref_id])


def test_nsgrid_PBC(universe):
    """Check that grid search works when PBC is needed"""

    ref_id = 13937
    results = np.array([4398, 4401, 13939, 13940, 13941, 17987, 23518, 23519, 23521, 23734,
                        47451]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!

    results_grid = run_grid_search(universe, ref_id).get_indices()[0]

    assert_equal(results, results_grid)



def test_nsgrid_pairs(universe):
    """Check that grid search returns the proper pairs"""

    ref_id = 13937
    neighbors = np.array([4398, 4401, 13939, 13940, 13941, 17987, 23518, 23519, 23521, 23734,
                          47451]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!
    results = []
    for nid in neighbors:
        if nid < ref_id:
            results.append([nid, ref_id])
        else:
            results.append([ref_id, nid])
    results = np.array(results)

    results_grid = run_grid_search(universe, ref_id).get_pairs()

    assert_equal(np.sort(results, axis=0), np.sort(results_grid, axis=0))


def test_nsgrid_pair_distances(universe):
    """Check that grid search returns the proper pair distances"""

    ref_id = 13937
    results = np.array([0.270, 0.285, 0.096, 0.096, 0.015, 0.278, 0.268, 0.179, 0.259, 0.290,
                        0.270]) * 10  # These distances where obtained by gmx distance so they are in nm

    results_grid = run_grid_search(universe, ref_id).get_pair_distances()

    assert_allclose(np.sort(results), np.sort(results_grid), atol=1e-2)


def test_nsgrid_pair_coordinates(universe):
    """Check that grid search return the proper pair coordinates"""

    ref_id = 13937
    neighbors = np.array([4398, 4401, 13939, 13940, 13941, 17987, 23518, 23519, 23521, 23734,
                          47451]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!
    coords = universe.atoms.positions

    results = []
    for nid in neighbors:
        if nid < ref_id:
            results.append([coords[nid], coords[ref_id]])
        else:
            results.append([coords[ref_id], coords[nid]])
    results = np.array(results)

    results_grid = run_grid_search(universe, ref_id).get_pair_coordinates()

    assert_allclose(np.sort(results, axis=0), np.sort(results_grid, axis=0), atol=1e-5)


def test_nsgrid_distances(universe):
    """Check that grid search returns the proper distances"""

    ref_id = 13937
    results = np.array([0.270, 0.285, 0.096, 0.096, 0.015, 0.278, 0.268, 0.179, 0.259, 0.290,
                        0.270]) * 10  # These distances where obtained by gmx distance so they are in nm

    results_grid = run_grid_search(universe, ref_id).get_distances()[0]

    assert_allclose(np.sort(results), np.sort(results_grid), atol=1e-2)


def test_nsgrid_coordinates(universe):
    """Check that grid search return the proper coordinates"""

    ref_id = 13937
    neighbors = np.array([4398, 4401, 13939, 13940, 13941, 17987, 23518, 23519, 23521, 23734,
                          47451]) - 1  # Atomid are from gmx select so there start from 1 and not 0. hence -1!

    results = universe.atoms.positions[neighbors]

    results_grid = run_grid_search(universe, ref_id).get_coordinates()[0]

    assert_allclose(np.sort(results, axis=0), np.sort(results_grid, axis=0), atol=1e-5)

