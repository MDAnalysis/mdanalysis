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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import os

import pytest

from collections import defaultdict, Counter
from numpy.testing import assert_equal, assert_allclose
import numpy as np

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (
    GRO, Martini_membrane_gro, PDB, PDB_xvf, SURFACE_PDB, SURFACE_TRR
)
from MDAnalysis.lib import nsgrid
from MDAnalysis.transformations.translate import center_in_box


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

@pytest.mark.parametrize('box', [
    np.zeros(3),  # Bad shape
    np.zeros((3, 3)),  # Collapsed box
    np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1]]),  # 2D box
    np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),  # Box provided as array of integers
    np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64),  # Box provided as array of double
])
def test_pbc_box(box):
    """Check that PBC box accepts only well-formated boxes"""
    coords = np.array([[1.0, 1.0, 1.0]], dtype=np.float32)

    with pytest.raises(ValueError):
        nsgrid.FastNS(4.0, coords, box=box)


@pytest.mark.parametrize('cutoff, match', ((-4, "Cutoff must be positive"),
                                           (100000,
                                            "Cutoff 100000 too large for box")))
def test_nsgrid_badcutoff(universe, cutoff, match):
    with pytest.raises(ValueError, match=match):
        run_grid_search(universe, 0, cutoff)


def test_ns_grid_noneighbor(universe):
    """Check that grid search returns empty lists/arrays when there is no neighbors"""
    ref_id = 0
    cutoff = 0.5

    results_grid = run_grid_search(universe, ref_id, cutoff)

    # same indices will be selected as neighbour here
    assert len(results_grid.get_pairs()) == 1
    assert len(results_grid.get_pair_distances()) == 1


def test_nsgrid_PBC_rect():
    """Check that nsgrid works with rect boxes and PBC"""
    ref_id = 191
    # Atomid are from gmx select so there start from 1 and not 0. hence -1!
    results = np.array([191, 192, 672, 682, 683, 684, 995, 996, 2060, 2808, 3300, 3791,
                        3792]) - 1

    universe = mda.Universe(Martini_membrane_gro)
    cutoff = 7

    # FastNS is called differently to max coverage
    searcher = nsgrid.FastNS(cutoff, universe.atoms.positions, box=universe.dimensions)

    results_grid = searcher.search(universe.atoms.positions[ref_id][None, :]).get_pairs()
    other_ix = sorted(i for (_, i) in results_grid)

    assert len(results) == len(results_grid)
    assert other_ix == sorted(results)


def test_nsgrid_PBC(universe):
    """Check that grid search works when PBC is needed"""
    # Atomid are from gmx select so there start from 1 and not 0. hence -1!
    ref_id = 13937
    results = np.array([4398, 4401, 13938, 13939, 13940, 13941, 17987, 23518, 23519, 23521, 23734,
                        47451]) - 1

    results_grid = run_grid_search(universe, ref_id).get_pairs()

    other_ix = sorted(i for (_, i) in results_grid)

    assert len(results) == len(other_ix)
    assert other_ix == sorted(results)


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
    # These distances where obtained by gmx distance so they are in nm
    ref_id = 13937
    results = np.array([0.0, 0.270, 0.285, 0.096, 0.096, 0.015, 0.278, 0.268, 0.179, 0.259, 0.290,
                        0.270]) * 10

    results_grid = run_grid_search(universe, ref_id).get_pair_distances()

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
    indices = searchresults.get_pairs()[:, 1]
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
    assert_equal(len(pairs), result)

def test_nsgrid_probe_close_to_box_boundary():
    # FastNS.search used to segfault with this box, cutoff and reference
    # coordinate prior to PR #2136, so we ensure that this remains fixed.
    # See Issue #2132 for further information.
    ref = np.array([[55.783722, 44.190044, -54.16671]], dtype=np.float32)
    box = np.array([53.785854, 43.951054, 57.17597, 90., 90., 90.], dtype=np.float32)
    cutoff = 3.0
    # search within a configuration where we know the expected outcome:
    conf = np.ones((1, 3), dtype=np.float32)
    searcher = nsgrid.FastNS(cutoff, conf, box)
    results = searcher.search(ref)
    # check if results are as expected:
    expected_pairs = np.zeros((1, 2), dtype=np.int64)
    expected_dists = np.array([2.3689647], dtype=np.float64)
    assert_equal(results.get_pairs(), expected_pairs)
    assert_allclose(results.get_pair_distances(), expected_dists, rtol=1.e-6)


def test_zero_max_dist():
    # see issue #2656
    # searching with max_dist = 0.0 shouldn't cause segfault (and infinite subboxes)
    ref = np.array([1.0, 1.0, 1.0], dtype=np.float32)
    conf = np.array([2.0, 1.0, 1.0], dtype=np.float32)

    box = np.array([10., 10., 10., 90., 90., 90.], dtype=np.float32)

    res = mda.lib.distances._nsgrid_capped(ref, conf, box=box, max_cutoff=0.0)


@pytest.fixture()
def u_pbc_triclinic():
    u = mda.Universe(PDB)
    u.dimensions = [10, 10, 10, 60, 60, 60]
    return u


def test_around_res(u_pbc_triclinic):
    # sanity check for issue 2656, shouldn't segfault (obviously)
    ag = u_pbc_triclinic.select_atoms('around 0.0 resid 3')
    assert len(ag) == 0


def test_around_overlapping():
    # check that around 0.0 catches when atoms *are* superimposed
    u = mda.Universe.empty(60, trajectory=True)
    xyz = np.zeros((60, 3))
    x = np.tile(np.arange(12), (5,))+np.repeat(np.arange(5)*100, 12)
    # x is 5 images of 12 atoms

    xyz[:, 0] = x  # y and z are 0
    u.load_new(xyz)

    u.dimensions = [100, 100, 100, 60, 60, 60]
    # Technically true but not what we're testing:
    # dist = mda.lib.distances.distance_array(u.atoms[:12].positions,
    #                                         u.atoms[12:].positions,
    #                                         box=u.dimensions)
    # assert np.count_nonzero(np.any(dist <= 0.0, axis=0)) == 48
    assert u.select_atoms('around 0.0 index 0:11').n_atoms == 48


def test_issue_2229_part1():
    # reproducing first case in GH issue 2229
    u = mda.Universe.empty(2, trajectory=True)

    u.dimensions = [57.45585, 50.0000, 50.0000, 90, 90, 90]

    u.atoms[0].position = [0, 0, 0]
    u.atoms[1].position = [55.00, 0, 0]

    g = mda.lib.nsgrid.FastNS(3.0, u.atoms[[0]].positions, box=u.dimensions)
    assert len(g.search(u.atoms[[1]].positions).get_pairs()) == 1

    g = mda.lib.nsgrid.FastNS(3.0, u.atoms[[1]].positions, box=u.dimensions)
    assert len(g.search(u.atoms[[0]].positions).get_pairs()) == 1


def test_issue_2229_part2():
    u = mda.Universe.empty(2, trajectory=True)

    u.dimensions = [45.0000, 55.0000, 109.8375, 90, 90, 90]

    u.atoms[0].position = [0, 0, 29.29]
    u.atoms[1].position = [0, 0, 28.23]

    g = mda.lib.nsgrid.FastNS(3.0, u.atoms[[0]].positions, box=u.dimensions, pbc=False)
    assert len(g.search(u.atoms[[1]].positions).get_pairs()) == 1

    g = mda.lib.nsgrid.FastNS(3.0, u.atoms[[1]].positions, box=u.dimensions)
    assert len(g.search(u.atoms[[0]].positions).get_pairs()) == 1


def test_issue_2919():
    # regression test reported in issue 2919
    # other methods will also give 1115 or 2479 results
    u = mda.Universe(PDB_xvf)
    ag = u.select_atoms('index 0')
    u.trajectory.ts = center_in_box(ag)(u.trajectory.ts)

    box = u.dimensions
    reference = u.select_atoms('protein')
    configuration = u.select_atoms('not protein')

    for cutoff, expected in [(2.8, 1115), (3.2, 2497)]:
        pairs, distances = mda.lib.distances.capped_distance(
            reference.positions,
            configuration.positions,
            max_cutoff=cutoff,
            box=box,
            method='nsgrid',
            return_distances=True,
        )
        assert len(pairs) == expected


def test_issue_2345():
    # another example of NSGrid being wrong
    # this is a 111 FCC slab
    # coordination numbers for atoms should be either 9 or 12, 50 of each
    u = mda.Universe(SURFACE_PDB, SURFACE_TRR)

    g = mda.lib.nsgrid.FastNS(2.9, u.atoms.positions, box=u.dimensions)

    cn = defaultdict(list)

    idx = g.self_search().get_pairs()
    # count number of contacts for each atom
    for (i, j) in idx:
        cn[i].append(j)
        cn[j].append(i)
    c = Counter(len(v) for v in cn.values())

    assert c == {9: 50, 12: 50}


def test_issue_2670():
    # Tests that NSGrid no longer crashes when using small box sizes
    u = mda.Universe(PDB)
    u.dimensions = [1e-3, 1e-3, 1e-3, 90, 90, 90]

    # PDB files only have a coordinate precision of 1.0e-3, so we need to scale
    # the coordinates for this test to make any sense:
    u.atoms.positions = u.atoms.positions * 1.0e-3

    ag1 = u.select_atoms('resid 2 3')
    # should return nothing as nothing except resid 3 is within 0.0 or resid 3
    assert len(ag1.select_atoms('around 0.0 resid 3')) == 0

    # force atom 0 of resid 1 to overlap with atom 0 of resid 3
    u.residues[0].atoms[0].position = u.residues[2].atoms[0].position
    ag2 = u.select_atoms('resid 1 3')

    # should return the one atom overlap
    assert len(ag2.select_atoms('around 0.0 resid 3')) == 1


def high_mem_tests_enabled():
    """ Returns true if ENABLE_HIGH_MEM_UNIT_TESTS is set to true."""
    env = os.getenv("ENABLE_HIGH_MEM_UNIT_TESTS", default="false").lower()
    if env == 'true':
        return True
    return False


reason = ("Turned off by default. The test can be enabled by setting "
          "the ENABLE_HIGH_MEM_UNIT_TESTS "
          "environment variable. Make sure you have at least 10GB of RAM.")


# Tests that with a tiny cutoff to box ratio, the number of grids is capped
# to avoid indexing overflow. Expected results copied from test_nsgrid_search
# with no box.
@pytest.mark.skipif(not high_mem_tests_enabled(), reason=reason)
def test_issue_3183():
    np.random.seed(90003)
    points = (np.random.uniform(low=0, high=1.0,
                                size=(100, 3)) * (10.)).astype(np.float32)
    cutoff = 2.0
    query = np.array([1., 1., 1.], dtype=np.float32).reshape((1, 3))
    box = np.array([10000., 10000., 10000., 90., 90., 90.])

    searcher = nsgrid.FastNS(cutoff, points, box)
    searchresults = searcher.search(query)
    indices = searchresults.get_pairs()[:, 1]
    want_results = [3, 13, 24]
    assert_equal(np.sort(indices), want_results)
