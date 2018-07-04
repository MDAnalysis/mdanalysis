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
from MDAnalysis.lib import grid
from MDAnalysis.lib.pkdtree import PeriodicKDTree

from MDAnalysis.lib.mdamath import triclinic_vectors

from MDAnalysisTests.datafiles import GRO


@pytest.fixture
def universe():
    u = mda.Universe(GRO)
    return u


def run_search(universe, ref_id):
    cutoff = 3

    coords = universe.atoms.positions
    ref_pos = coords[ref_id]
    triclinic_box = triclinic_vectors(universe.dimensions)


    # Run pkdtree search
    pkdt = PeriodicKDTree(universe.atoms.dimensions, bucket_size=10)
    pkdt.set_coords(coords)
    pkdt.search(ref_pos, cutoff)

    results_pkdtree = pkdt.get_indices()
    results_pkdtree.remove(ref_id)
    results_pkdtree = np.array(results_pkdtree)
    results_pkdtree.sort()

    # Run grid search

    searcher = grid.FastNS(triclinic_box)
    searcher.set_cutoff(cutoff)
    searcher.set_coords(coords)
    searcher.prepare()

    results_grid = searcher.search(np.array([ref_pos, ]), return_ids=True)[0]
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