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
from numpy.testing import assert_equal, assert_almost_equal

import MDAnalysis as mda

@pytest.mark.parametrize('coord_dtype', (np.float32, np.float64))
def test_transform_StoR_pass(coord_dtype):
    box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
    s = np.array([[0.5, -0.1, 0.5]], dtype=coord_dtype)

    original_r = np.array([[ 5.75,  0.36066014, 0.75000012]], dtype=np.float32)

    test_r = mda.lib.distances.transform_StoR(s, box)

    assert_equal(original_r, test_r)


def test_capped_distance_noresults():
    point1 = np.array([0.1, 0.1, 0.1], dtype=np.float32)
    point2 = np.array([0.95, 0.1, 0.1], dtype=np.float32)

    pairs, distances = mda.lib.distances.capped_distance(point1,
                                                        point2,
                                                        max_cutoff=0.2)

    assert_equal(len(pairs), 0)


npoints_1 = (1, 100)

boxes_1 = (np.array([1, 2, 3, 90, 90, 90], dtype=np.float32),  # ortho
           np.array([1, 2, 3, 30, 45, 60], dtype=np.float32),  # tri_box
           None,  # Non Periodic
           )


query_1 = (np.array([0.1, 0.1, 0.1], dtype=np.float32),
           np.array([[0.1, 0.1, 0.1],
                     [0.2, 0.1, 0.1]], dtype=np.float32))

method_1 = ('bruteforce', 'pkdtree', 'nsgrid')

min_cutoff_1 = (None, 0.1)


@pytest.mark.parametrize('npoints', npoints_1)
@pytest.mark.parametrize('box', boxes_1)
@pytest.mark.parametrize('query', query_1)
@pytest.mark.parametrize('method', method_1)
@pytest.mark.parametrize('min_cutoff', min_cutoff_1)
def test_capped_distance_checkbrute(npoints, box, query, method, min_cutoff):
    np.random.seed(90003)
    points = (np.random.uniform(low=0, high=1.0,
                        size=(npoints, 3))*(boxes_1[0][:3])).astype(np.float32)
    max_cutoff = 0.3
    # capped distance should be able to handle array of vectors
    # as well as single vectors.
    pairs, dist = mda.lib.distances.capped_distance(query,
                                                    points,
                                                    max_cutoff,
                                                    min_cutoff=min_cutoff,
                                                    box=box,
                                                    method=method)

    if pairs.shape != (0, ):
        found_pairs = pairs[:, 1]
    else:
        found_pairs = list()

    if(query.shape[0] == 3):
        query = query.reshape((1, 3))

    distances = mda.lib.distances.distance_array(query,
                                                 points, box=box)

    if min_cutoff is None:
        min_cutoff = 0.
    indices = np.where((distances < max_cutoff) & (distances > min_cutoff))

    assert_equal(np.sort(found_pairs, axis=0), np.sort(indices[1], axis=0))

@pytest.mark.parametrize('npoints', npoints_1)
@pytest.mark.parametrize('box', boxes_1)
@pytest.mark.parametrize('method', method_1)
@pytest.mark.parametrize('min_cutoff', min_cutoff_1)
def test_self_capped_distance(npoints, box, method, min_cutoff):
    np.random.seed(90003)
    points = (np.random.uniform(low=0, high=1.0,
                         size=(npoints, 3))*(boxes_1[0][:3])).astype(np.float32)
    max_cutoff = 0.2
    pairs, distance = mda.lib.distances.self_capped_distance(points,
                                                             max_cutoff,
                                                             min_cutoff=min_cutoff,
                                                             box=box,
                                                             method=method)
    found_pairs, found_distance = [], []
    for i, coord in enumerate(points):
        dist = mda.lib.distances.distance_array(coord[None, :],
                                                points[i+1:],
                                                box=box)
        if min_cutoff is not None:
            idx = np.where((dist < max_cutoff) & (dist > min_cutoff))[1]
        else:
            idx = np.where((dist < max_cutoff))[1]
        for other_idx in idx:
            j = other_idx + 1 + i
            found_pairs.append((i, j))
            found_distance.append(dist[0, other_idx])
    assert_equal(len(pairs), len(found_pairs))


@pytest.mark.parametrize('box', (None,
                                 np.array([1, 1, 1,  90, 90, 90], dtype=np.float32),
                                 np.array([1, 1, 1, 60, 75, 80], dtype=np.float32)))
@pytest.mark.parametrize('npoints,cutoff,meth',
                         [(1, 0.02, '_bruteforce_capped_self'),
                          (1, 0.2, '_bruteforce_capped_self'),
                          (6000, 0.02, '_pkdtree_capped_self'),
                          (6000, 0.2, '_pkdtree_capped_self'),
                          (200000, 0.02, '_pkdtree_capped_self'),
                          (200000, 0.2, '_bruteforce_capped_self')])
def test_method_selfselection(box, npoints, cutoff, meth):
    np.random.seed(90003)
    points = (np.random.uniform(low=0, high=1.0,
                        size=(npoints, 3))).astype(np.float32)
    method = mda.lib.distances._determine_method_self(points, cutoff, box=box)
    assert_equal(method.__name__, meth)


@pytest.mark.parametrize('box', (None,
                                 np.array([1, 1, 1,  90, 90, 90], dtype=np.float32),
                                 np.array([1, 1, 1, 60, 75, 80], dtype=np.float32)))
@pytest.mark.parametrize('npoints,cutoff,meth',
                         [(1, 0.02, '_bruteforce_capped'),
                          (1, 0.2, '_bruteforce_capped'),
                          (6000, 0.02, '_pkdtree_capped'),
                          (6000, 0.2, '_pkdtree_capped'),
                          (200000, 0.02, '_pkdtree_capped'),
                          (200000, 0.2, '_bruteforce_capped')])
def test_method_selection(box, npoints, cutoff, meth):
    np.random.seed(90003)
    points = (np.random.uniform(low=0, high=1.0,
                        size=(npoints, 3)).astype(np.float32))
    method = mda.lib.distances._determine_method(points, points, cutoff, box=box)
    assert_equal(method.__name__, meth)


# different boxlengths to shift a coordinate
shifts = [
    ((0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)),  # no shifting
    ((1, 0, 0), (0, 1, 1), (0, 0, 1), (1, 1, 0)),  # single box lengths
    ((-1, 0, 1), (0, -1, 0), (1, 0, 1), (-1, -1, -1)),  # negative single
    ((4, 3, -2), (-2, 2, 2), (-5, 2, 2), (0, 2, 2)),  # multiple boxlengths
]


@pytest.mark.parametrize('shift', shifts)
@pytest.mark.parametrize('periodic', [True, False])
def test_calc_distance(shift, periodic):
    box = np.array([10, 20, 30, 90., 90., 90.], dtype=np.float32)

    coords = np.array([[1, 1, 1], [3, 1, 1]], dtype=np.float32)

    shift1, shift2, _, _ = shift

    coords[0] += shift1 * box[:3]
    coords[1] += shift2 * box[:3]

    box = box if periodic else None
    result = mda.lib.distances.calc_distance(coords[0], coords[1], box)

    reference = 2.0 if periodic else np.linalg.norm(coords[0] - coords[1])

    assert_almost_equal(result, reference, decimal=3)


@pytest.mark.parametrize('case', [
    # 90 degree angle
    (np.array([[1, 1, 1], [1, 2, 1], [2, 2, 1]], dtype=np.float32), 90.),
    # straight line / 180.
    (np.array([[1, 1, 1], [1, 2, 1], [1, 3, 1]], dtype=np.float32), 180.),
    # 45
    (np.array([[1, 1, 1], [1, 2, 1], [2, 1, 1]], dtype=np.float32), 45.),
])
@pytest.mark.parametrize('shift', shifts)
@pytest.mark.parametrize('periodic', [True, False])
def test_calc_angle(case, shift, periodic):
    def manual_angle(x, y, z):
        return np.rad2deg(mda.lib.mdamath.angle(y - x, y - z))

    box = np.array([10, 20, 30, 90., 90., 90.], dtype=np.float32)
    (a, b, c), ref = case

    shift1, shift2, shift3, _ = shift

    a += shift1 * box[:3]
    b += shift2 * box[:3]
    c += shift3 * box[:3]

    box = box if periodic else None
    result = mda.lib.distances.calc_angle(a, b, c, box)

    reference = ref if periodic else manual_angle(a, b, c)

    assert_almost_equal(result, reference, decimal=3)


@pytest.mark.parametrize('case', [
    # 0 degree angle (cis)
    (np.array([[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 2, 1]], dtype=np.float32), 0.),
    # 180 degree (trans)
    (np.array([[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 0, 1]], dtype=np.float32), 180.),
    # 90 degree
    (np.array([[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 1, 2]], dtype=np.float32), 90.),
    # other 90 degree
    (np.array([[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 1, 0]], dtype=np.float32), 90.),
    # 45 degree
    (np.array([[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 2, 2]], dtype=np.float32), 45.),
    # 135
    (np.array([[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 0, 2]], dtype=np.float32), 135.),
])
@pytest.mark.parametrize('shift', shifts)
@pytest.mark.parametrize('periodic', [True, False])
def test_calc_dihedral(case, shift, periodic):
    def manual_dihedral(a, b, c, d):
        return np.rad2deg(mda.lib.mdamath.dihedral(b - a, c - b, d - c))

    box = np.array([10., 10., 10., 90., 90., 90.], dtype=np.float32)

    (a, b, c, d), ref = case

    shift1, shift2, shift3, shift4 = shift

    a += shift1 * box[:3]
    b += shift2 * box[:3]
    c += shift3 * box[:3]
    d += shift4 * box[:3]

    box = box if periodic else None
    result = mda.lib.distances.calc_dihedral(a, b, c, d, box)

    reference = ref if periodic else manual_dihedral(a, b, c, d)

    assert_almost_equal(abs(result), abs(reference), decimal=3)
