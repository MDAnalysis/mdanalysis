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

from __future__ import print_function, absolute_import
from six.moves import zip

import pytest
import numpy as np
from numpy.testing import assert_equal, assert_almost_equal

from MDAnalysis.lib.pkdtree import PeriodicKDTree
from MDAnalysis.lib.mdamath import triclinic_vectors, triclinic_box


boxes = (np.array([1, 2, 3, 90, 90, 90], dtype=np.float32),  # ortho
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


@pytest.mark.parametrize('box, rm', zip(boxes, rec_m))
def test_initialize_bm(box, rm):
    assert_almost_equal(PeriodicKDTree(box)._rm, rm, decimal=7)


@pytest.fixture
def ptree_ortho():
    b = np.array([10, 10, 10, 90, 90, 90], dtype=np.float32)
    coords = np.array([[2, 2, 2],
                       [5, 5, 5],
                       [1.1, 1.1, 1.1],
                       [11, -11, 11],  # wrapped to [1, 9, 1]
                       [21, 21, 3]],  # wrapped to [1, 1, 3]
                      dtype=np.float32)
    t = PeriodicKDTree(b)
    t.set_coords(coords)
    return {'coords': coords, 'tree': t, 'box': b, 'radius': 1.5}


def test_set_coords(ptree_ortho):
    with pytest.raises(ValueError) as excinfo:
        xy = np.array([[2, 2], [5, 5], [1.1, 1.1]], dtype=np.float32)
        tree = PeriodicKDTree(ptree_ortho['box'])
        tree.set_coords(xy)
    assert_equal(str(excinfo.value),
                 'coords must be a sequence of 3 dimensional coordinates')


queries = ([5, 5, 5],  # case box center
           [1, 5, 5],  # box face
           [5, -1, 5],  # box face
           [1, 1, 5],  # box edge
           [5, -1, 11],  # box edge
           [1, 1, 1],  # box vertex
           [1, -1, 11],  # box vertex
           [21, -31, 1]  # box vertex
           )
centers = (([5, 5, 5], ),
           ([1, 5, 5], [11, 5, 5]),  # centers for first case box face
           ([5, 9, 5], [5, -1, 5]),
           ([1, 1, 5], [11, 1, 5], [1, 11, 5], [11, 11, 5]),
           ([5, 9, 1], [5, -1, 1], [5, 9, 11], [5, -1, 11]),
           ([1, 1, 1], [11, 1, 1], [1, 11, 1], [1, 1, 11],
            [1, 11, 11], [11, 11, 1], [11, 1, 11], [11, 11, 11]),
           ([1, 9, 1], [11, 9, 1], [1, -1, 1], [1, 9, 11],
            [1, -1, 11], [11, -1, 1], [11, 9, 11], [11, -1, 11]),
           ([1, 9, 1], [11, 9, 1], [1, -1, 1], [1, 9, 11],
            [1, -1, 11], [11, -1, 1], [11, 9, 11], [11, -1, 11])
           )


@pytest.mark.parametrize('q, cs', zip(queries, centers))
def test_find_centers(ptree_ortho, q, cs):
    q = np.array(q, dtype=np.float32)
    cs = [np.array(c, dtype=np.float32) for c in cs]
    assert_equal(ptree_ortho['tree'].find_centers(q,
                                                  ptree_ortho['radius']), cs)


queries = ([5, 5, 5],  # case box center
           [-8.5, 11.5, 2.2],  # wrapped to [1.5, 1.5, 2.2]
           [0, 100, 0.7],  # box face
           [1, 1, 5],  # box edge
           [1, 1, 1],  # box vertex
           [-19, 42, 2],  # box vertex
           [21, -31, 1]  # box vertex
           )
neighbors = (([5, 5, 5], ),
             ([2, 2, 2], [1.1, 1.1, 1.1], [21, 21, 3]),
             ([11, -11, 11], ),
             (),
             ([1.1, 1.1, 1.1], ),
             ([2, 2, 2], [1.1, 1.1, 1.1], [21, 21, 3]),
             ([11, -11, 11], )
             )


@pytest.mark.parametrize('q, ns', zip(queries, neighbors))
def test_search(ptree_ortho, q, ns):
    ptree_ortho['tree'].search(np.array(q, dtype=np.float32),
                               ptree_ortho['radius'])
    indices = ptree_ortho['tree'].get_indices()
    found_neighbors = list() if indices is None \
        else [ptree_ortho['coords'][i] for i in indices]
    ns = [np.array(n, dtype=np.float32) for n in ns]
    assert_equal(found_neighbors, ns)
