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


class TestCheckBox(object):

    prec = 6
    ref_ortho = np.ones(3, dtype=np.float32)
    ref_tri_vecs = np.array([[1, 0, 0], [0, 1, 0], [0, 2 ** 0.5, 2 ** 0.5]],
                            dtype=np.float32)

    @pytest.mark.parametrize('box',
        ([1, 1, 1, 90, 90, 90],
         (1, 1, 1, 90, 90, 90),
         ['1', '1', 1, 90, '90', '90'],
         ('1', '1', 1, 90, '90', '90'),
         np.array(['1', '1', 1, 90, '90', '90']),
         np.array([1, 1, 1, 90, 90, 90], dtype=np.float32),
         np.array([1, 1, 1, 90, 90, 90], dtype=np.float64),
         np.array([1, 1, 1, 1, 1, 1, 90, 90, 90, 90, 90, 90],
                  dtype=np.float32)[::2]))
    def test_ckeck_box_ortho(self, box):
        boxtype, checked_box = mda.lib.distances._check_box(box)
        assert boxtype == 'ortho'
        assert_equal(checked_box, self.ref_ortho)
        assert checked_box.dtype == np.float32
        assert checked_box.flags['C_CONTIGUOUS']

    @pytest.mark.parametrize('box',
         ([1, 1, 2, 45, 90, 90],
          (1, 1, 2, 45, 90, 90),
          ['1', '1', 2, 45, '90', '90'],
          ('1', '1', 2, 45, '90', '90'),
          np.array(['1', '1', 2, 45, '90', '90']),
          np.array([1, 1, 2, 45, 90, 90], dtype=np.float32),
          np.array([1, 1, 2, 45, 90, 90], dtype=np.float64),
          np.array([1, 1, 1, 1, 2, 2, 45, 45, 90, 90, 90, 90],
                   dtype=np.float32)[::2]))
    def test_check_box_tri_vecs(self, box):
        boxtype, checked_box = mda.lib.distances._check_box(box)
        assert boxtype == 'tri_vecs'
        assert_almost_equal(checked_box, self.ref_tri_vecs, self.prec)
        assert checked_box.dtype == np.float32
        assert checked_box.flags['C_CONTIGUOUS']

    def test_check_box_wrong_data(self):
        with pytest.raises(ValueError):
            wrongbox = ['invalid', 1, 1, 90, 90, 90]
            boxtype, checked_box = mda.lib.distances._check_box(wrongbox)

    def test_check_box_wrong_shape(self):
        with pytest.raises(ValueError):
            wrongbox = np.ones((3, 3), dtype=np.float32)
            boxtype, checked_box = mda.lib.distances._check_box(wrongbox)


class TestCheckResultArray(object):

    ref = np.zeros(1, dtype=np.float64)

    def test_check_result_array_pass(self):
        # Assert input array is returned if it has correct shape and dtype:
        res = mda.lib.distances._check_result_array(self.ref, self.ref.shape)
        assert res is self.ref
        # Assert correct array is returned if input is None:
        res = mda.lib.distances._check_result_array(None, self.ref.shape)
        assert_equal(res, self.ref)
        assert res.dtype == np.float64

    def test_check_result_array_wrong_shape(self):
        wrong_shape = (1,) + self.ref.shape
        with pytest.raises(ValueError) as err:
            res = mda.lib.distances._check_result_array(self.ref, wrong_shape)
            assert err.msg == ("Result array has incorrect shape, should be "
                               "{0}, got {1}.".format(self.ref.shape,
                                                      wrong_shape))

    def test_check_result_array_wrong_dtype(self):
        wrong_dtype = np.int64
        ref_wrong_dtype = self.ref.astype(wrong_dtype)
        with pytest.raises(TypeError) as err:
            res = mda.lib.distances._check_result_array(ref_wrong_dtype,
                                                        self.ref.shape)
            assert err.msg == ("Result array must be of type numpy.float64, "
                               "got {}.".format(wrong_dtype))


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

# for coverage
@pytest.mark.parametrize('npoints', npoints_1)
@pytest.mark.parametrize('box', boxes_1)
@pytest.mark.parametrize('query', query_1)
@pytest.mark.parametrize('method', method_1)
@pytest.mark.parametrize('min_cutoff', min_cutoff_1)
def test_capped_distance_return(npoints, box, query, method, min_cutoff):
    np.random.seed(90003)
    points = (np.random.uniform(low=0, high=1.0,
                        size=(npoints, 3))*(boxes_1[0][:3])).astype(np.float32)
    max_cutoff = 0.3
    # capped distance should be able to handle array of vectors
    # as well as single vectors.
    pairs = mda.lib.distances.capped_distance(query,
                                                    points,
                                                    max_cutoff,
                                                    min_cutoff=min_cutoff,
                                                    box=box,
                                                    method=method,
                                                    return_distances=False)

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
                          (600, 0.02, '_pkdtree_capped_self'),
                          (600, 0.2, '_nsgrid_capped_self')])
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
                          (200, 0.02, '_nsgrid_capped'),
                          (200, 0.35, '_bruteforce_capped'),
                          (10000, 0.35, '_nsgrid_capped')])
def test_method_selection(box, npoints, cutoff, meth):
    np.random.seed(90003)
    points = (np.random.uniform(low=0, high=1.0,
                        size=(npoints, 3)).astype(np.float32))
    method = mda.lib.distances._determine_method(points, points, cutoff, box=box)
    assert_equal(method.__name__, meth)
