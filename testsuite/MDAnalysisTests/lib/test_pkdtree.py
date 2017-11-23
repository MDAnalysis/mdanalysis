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
from six.moves import zip
from itertools import product

import pytest
import numpy as np
from numpy.testing import assert_equal, assert_almost_equal

from MDAnalysis.lib.pkdtree import PeriodicKDTree
from MDAnalysis.lib.mdamath import triclinic_vectors, triclinic_box
from MDAnalysis.lib.distances import (_box_check, transform_RtoS,
                                      transform_StoR, apply_PBC)

#
# Testing initialization with different boxes.
#
boxes_1 = (np.array([1, 2, 3, 90, 90, 90], dtype=np.float32),  # ortho
         np.array([1, 2, 3, 30, 45, 60], dtype=np.float32),  # tri_box
         triclinic_vectors(  # tri_vecs
             np.array([1, 2, 3, 90, 90, 45], dtype=np.float32)),
         np.array([[0.5, 0.9, 1.9],  # tri_vecs_bad
                   [2.0, 0.4, 0.1],
                   [0.0, 0.6, 0.5]], dtype=np.float32)
         )
rec_m = (np.array([[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, 1]], dtype=np.float32),
         np.array([[0.67044002, -0.38707867, -0.6329931],
                   [0, 0.54741204, -0.83686334],
                   [0, 0, 1]], dtype=np.float32),
         np.array([[0.707106829, -0.707106829, 0],
                   [0, 1, 0],
                   [0, 0, 1]], dtype=np.float32),
         np.array([[0.42783618, -0.16049984, -0.88949198],
                   [-0, 0.95654362, 0.29158944],
                   [0, 0, 1]], dtype=np.float32),
         )


@pytest.mark.parametrize('box, rm', zip(boxes_1, rec_m))
def test_initialize_bm(box, rm):
    """
    Assert the construction of the recripocal box matrix.
    """
    assert_almost_equal(PeriodicKDTree(box)._rm, rm, decimal=7)


def test_set_coords():
    match = 'coords must be a sequence of 3 dimensional coordinates'
    with pytest.raises(ValueError, match=match):
        xy = np.array([[2, 2], [5, 5], [1.1, 1.1]], dtype=np.float32)
        tree = PeriodicKDTree(boxes_1[0])
        tree.set_coords(xy)


#
# Testing for correct images generation for a given query point.
#

# fractional coordinates for data points
f_dataset = np.array([[0.2, 0.2, 0.2],  # center of the box
                      [0.5, 0.5, 0.5],
                      [0.11, 0.11, 0.11],
                      [1.1, -1.1, 1.1],  # wrapped to [1, 9, 1]
                      [2.1, 2.1, 0.3]],  # wrapped to [1, 1, 3]
                     dtype=np.float32)

radius = 1.5

boxes_2 = ([10, 10, 10, 90, 90, 90],  # ortho
          [10, 10, 10, 45, 60, 90])  # tri_box

# Find images for several query points, here in fractional coordinates.
queries = ([0.1, 0.5, 0.5],  # box face
           [0.5, 0.5, 0.5],  # box center
           [0.5, -0.1, 0.5],  # box face
           [0.1, 0.1, 0.5],  # box edge
           [0.5, -0.1, 1.1],  # box edge
           [0.1, 0.1, 0.1],  # box vertex
           [0.1, -0.1, 1.1],  # box vertex
           [2.1, -3.1, 0.1]  # box vertex
           )

# Images for the previous query vectors, here in fractional coordinates.
images = (([1.1, 0.5, 0.5],),
          (),
          ([0.5, -0.1, 0.5],),
          ([1.1, 0.1, 0.5], [0.1, 1.1, 0.5], [1.1, 1.1, 0.5]),
          ([0.5, -0.1, 0.1], [0.5, 0.9, 1.1], [0.5, -0.1, 1.1]),
          ([1.1, 0.1, 0.1], [0.1, 1.1, 0.1], [0.1, 0.1, 1.1],
           [0.1, 1.1, 1.1], [1.1, 1.1, 0.1], [1.1, 0.1, 1.1],
           [1.1, 1.1, 1.1]),
          ([1.1, 0.9, 0.1], [0.1, -0.1, 0.1], [0.1, 0.9, 1.1],
           [0.1, -0.1, 1.1], [1.1, -0.1, 0.1], [1.1, 0.9, 1.1],
           [1.1, -0.1, 1.1]),
          ([1.1, 0.9, 0.1], [0.1, -0.1, 0.1], [0.1, 0.9, 1.1],
           [0.1, -0.1, 1.1], [1.1, -0.1, 0.1], [1.1, 0.9, 1.1],
           [1.1, -0.1, 1.1]))


@pytest.mark.parametrize('b, qcs', product(boxes_2, zip(queries, images)))
def test_find_images(b, qcs):
    """
    Test the generation of images for a given query vector and type of box.

    Parameters
    ----------
    b : list
        MDAnalysis dimensions like list
    qns : tuple
        a query point and a list of expected neighbors.
    """
    b = np.array(b, dtype=np.float32)
    q = transform_StoR(np.array(qcs[0], dtype=np.float32), b)
    tree = PeriodicKDTree(b)
    coords = transform_StoR(f_dataset, b)
    tree.set_coords(coords)  # Input real space coordinates
    cs = np.sort(transform_StoR(np.array(qcs[1], dtype=np.float32), b), axis=0)
    q_wrapped = apply_PBC(q.reshape((1, 3)), b)
    found_images = np.sort(tree.find_images(q_wrapped, radius), axis=0)
    assert_almost_equal(found_images, cs, decimal=6)


#
# Testing for neighbor finding
#

# Find neighbors for a given query vector, here in fractional coordinates.
queries_2 = ([0.5, 0.5, 0.5],  # case box center
             [0, 1, 0.07],  # box face
             [0.1, 0.1, 0.5],  # box edge
             [0.1, 0.1, 0.1],  # box vertex
             ([-1.9, 4.2, 0.2], [2.1, -3.1, 0.1])  # multiple queries
            )
# Expected neighbors of queries_2 in the orthogonal box case, here in
# fractional coordinates.
n_ortho = (([0.5, 0.5, 0.5],),
           ([1.1, -1.1, 1.1], ),
           (),
           ([0.11, 0.11, 0.11], ),
           ([0.2, 0.2, 0.2], [0.11, 0.11, 0.11], [2.1, 2.1, 0.3],
               [1.1, -1.1, 1.1])
           )
# Expected neighbors for previous queries in the trigonal box case
n_tric = (([0.5, 0.5, 0.5],),
          ([1.1, -1.1, 1.1], ),
          (),
          ([0.11, 0.11, 0.11], ),
          ([0.2, 0.2, 0.2], [2.1, 2.1, 0.3], [1.1, -1.1, 1.1])
          )

# Combinations of boxes, query points, and expected neighbors.
doublets = list()
for b, n in zip(boxes_2, (n_ortho, n_tric)):
    doublets.extend(list(product([b], zip(queries_2, n))))


@pytest.mark.parametrize('b, qns', doublets)
def test_search(b, qns):
    """
    Test finding neighbors for a given query vector and type of box.

    Parameters
    ----------
    b : list
        MDAnalysis dimensions like list
    qns : tuple
        a query point and a list of expected neighbors.
    """
    b = np.array(b, dtype=np.float32)
    q = transform_StoR(np.array(qns[0], dtype=np.float32), b)
    # Setting up the periodic tree
    tree = PeriodicKDTree(b)
    coords = transform_StoR(f_dataset, b)
    tree.set_coords(coords)  # Input real space coordinates
    # Carry out the search and retrieve results
    tree.search(q, radius)
    indices = tree.get_indices()
    if indices:
        found_neighbors = np.sort(coords[indices], axis=0)
    else:
        found_neighbors = list()
    if qns[1]:
        expected_neighbors = transform_StoR(np.array(qns[1],dtype=np.float32), b)
        expected_neighbors = np.sort(expected_neighbors, axis=0)
    else:
        expected_neighbors = list()
    assert_equal(found_neighbors, expected_neighbors)
