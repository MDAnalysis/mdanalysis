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
import pytest
import numpy as np
from numpy.testing import assert_equal

import itertools

import MDAnalysis as mda

from MDAnalysis.lib.mdamath import triclinic_vectors

def test_transform_StoR_pass():
    box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
    s = np.array([[0.5, -0.1, 0.5]], dtype=np.float32)

    original_r = np.array([[ 5.75,  0.36066014, 0.75000012]], dtype=np.float32)

    test_r = mda.lib.distances.transform_StoR(s, box)

    assert_equal(original_r, test_r)


def test_transform_StoR_fail():
    box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
    s = np.array([[0.5, -0.1, 0.5]])

    with pytest.raises(TypeError, match='S must be of type float32'):
        r = mda.lib.distances.transform_StoR(s, box)


def test_capped_distance_noresults():
    point1 = np.array([0.1, 0.1, 0.1], dtype=np.float32)
    point2 = np.array([0.95, 0.1, 0.1], dtype=np.float32)

    pairs, distances = mda.lib.distances.capped_distance(point1,
                                                        point2,
                                                        max_cutoff=0.2)

    assert_equal(len(pairs), 0)


boxes_1 = (np.array([1, 2, 3, 90, 90, 90], dtype=np.float32),  # ortho
           np.array([1, 2, 3, 30, 45, 60], dtype=np.float32),  # tri_box
           triclinic_vectors(  # tri_vecs
           np.array([1, 2, 3, 90, 90, 45], dtype=np.float32)),
           np.array([[0.5, 0.9, 1.9],  # tri_vecs_bad
                     [2.0, 0.4, 0.1],
                     [0.0, 0.6, 0.5]], dtype=np.float32))

query_1 = (np.array([0.1, 0.1, 0.1], dtype=np.float32),
           np.array([[0.1, 0.1, 0.1],
                     [0.2, 0.1, 0.1]], dtype=np.float32))

method_1 = ('bruteforce', 'pkdtree')

np.random.seed(90003)
points = (np.random.uniform(low=0, high=1.0,
                        size=(100, 3))*(boxes_1[0][:3])).astype(np.float32)


@pytest.mark.parametrize('box, query , method',
                         itertools.product(boxes_1, query_1, method_1))
def test_capped_distance_checkbrute(box, query, method):
    max_cutoff = 0.3
    # capped distance should be able to handle array of vectors
    # as well as single vectors.
    pairs, dist = mda.lib.distances.capped_distance(query,
                                                    points,
                                                    max_cutoff,
                                                    box=box,
                                                    method=method)
    if(query.shape[0] == 3):
        distances = mda.lib.distances.distance_array(query[None, :],
                                                     points, box=box)
    else:
        distances = mda.lib.distances.distance_array(query,
                                                     points, box=box)
    indices = np.where(distances < max_cutoff)

    assert_equal(pairs[:, 1], indices[1])
