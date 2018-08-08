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

from __future__ import print_function, absolute_import

import pytest
import numpy as np
from numpy.testing import assert_equal


from MDAnalysis.lib.pkdtree import PeriodicKDTree
from MDAnalysis.lib.distances import transform_StoR


# fractional coordinates for data points
f_dataset = np.array([[0.2, 0.2, 0.2],  # center of the box
                      [0.5, 0.5, 0.5],
                      [0.11, 0.11, 0.11],
                      [1.1, -1.1, 1.1],  # wrapped to [1, 9, 1]
                      [2.1, 2.1, 0.3]],  # wrapped to [1, 1, 3]
                     dtype=np.float32)


@pytest.mark.parametrize('b, cut, result', (
                         (None, 1.0,
                          'Donot provide cutoff distance'
                          ' for non PBC aware calculations'),
                         ([10, 10, 10, 90, 90, 90], None,
                          'Provide a cutoff distance with'
                          ' tree.set_coords(...)')))
def test_setcoords(b, cut, result):
    coords = np.array([[1, 1, 1], [2, 2, 2]], dtype=np.float32)
    if b is not None:
        b = np.array(b, dtype=np.float32)
    tree = PeriodicKDTree(box=b)
    print(b, tree.box, cut, result)
    with pytest.raises(RuntimeError, match=result):
        tree.set_coords(coords, cutoff=cut)


def test_searchfail():
    coords = np.array([[1, 1, 1], [2, 2, 2]], dtype=np.float32)
    b = np.array([10, 10, 10, 90, 90, 90], dtype=np.float32)
    cutoff = 1.0
    search_radius = 2.0
    query = np.array([1, 1, 1], dtype=np.float32)
    tree = PeriodicKDTree(box=b)
    tree.set_coords(coords, cutoff=cutoff)
    match = 'Set cutoff greater or equal to the radius.'
    with pytest.raises(RuntimeError, match=match):
        tree.search(query, search_radius)


@pytest.mark.parametrize('b, q, result', (
                         ([10, 10, 10, 90, 90, 90], [0.5, -0.1, 1.1], []),
                         ([10, 10, 10, 90, 90, 90], [2.1, -3.1, 0.1], [2, 3, 4]),
                         ([10, 10, 10, 45, 60, 90], [2.1, -3.1, 0.1], [2, 3])
                         ))
def test_search(b, q, result):
    b = np.array(b, dtype=np.float32)
    q = transform_StoR(np.array(q, dtype=np.float32), b)
    cutoff = 3.0
    coords = transform_StoR(f_dataset, b)
    tree = PeriodicKDTree(box=b)
    tree.set_coords(coords, cutoff=cutoff)
    indices = tree.search(q, cutoff)
    assert_equal(indices, result)


def test_nopbc():
    cutoff = 0.3
    q = np.array([0.2, 0.3, 0.1])
    coords = f_dataset.copy()
    tree = PeriodicKDTree(box=None)
    tree.set_coords(coords)
    indices = tree.search(q, cutoff)
    assert_equal(indices, [0, 2])


@pytest.mark.parametrize('b, radius, result', (
                         ([10, 10, 10, 90, 90, 90], 2.0,  [[0, 2],
                                                           [0, 4],
                                                           [2, 4]]),
                         ([10, 10, 10, 45, 60, 90], 2.0,  [[0, 4],
                                                           [2, 4]]),
                         ([10, 10, 10, 45, 60, 90], 4.5,
                          'Set cutoff greater or equal to the radius.'),
                         ([10, 10, 10, 45, 60, 90], 0.1, [])
                         ))
def test_searchpairs(b, radius, result):
    b = np.array(b, dtype=np.float32)
    cutoff = 2.0
    coords = transform_StoR(f_dataset, b)
    tree = PeriodicKDTree(box=b)
    tree.set_coords(coords, cutoff=cutoff)
    if cutoff < radius:
        with pytest.raises(RuntimeError, match=result):
            indices = tree.search_pairs(radius)
    else:
        indices = tree.search_pairs(radius)
        assert_equal(len(indices), len(result))

@pytest.mark.parametrize('radius, result', ((0.1, []),
                                            (0.3, [[0, 2]])))
def test_ckd_searchpairs_nopbc(radius, result):
    coords = f_dataset.copy()
    tree = PeriodicKDTree()
    tree.set_coords(coords)
    indices = tree.search_pairs(radius)
    assert_equal(indices, result)
