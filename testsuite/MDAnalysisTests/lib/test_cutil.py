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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

import pytest
import numpy as np
from numpy.testing import assert_equal

from MDAnalysis.lib._cutil import (coords_add_vec, unique_int_1d,
                                   unique_masks_int_1d, iscontiguous_int_1d,
                                   argwhere_int_1d, indices_to_slice_1d,
                                   find_fragments)


@pytest.mark.parametrize('positions', (
    np.zeros((0, 3), dtype=np.float32),
    np.array([[1.0000001, -1.0000001, 0.0000001]], dtype=np.float32),
    np.zeros((10, 3), dtype=np.float32),
    np.ones((10, 3), dtype=np.float32),
    np.reshape(np.arange(-30, 30, dtype=np.float32), (-1, 3))
))
@pytest.mark.parametrize('vector', (
    np.array([-1, 2, 3], dtype=np.float32),
    np.array([-1, 2, 3], dtype=np.float64),
    np.array([-1, 2, 3], dtype=int),
    np.array([-1, 0, 3], dtype=bool),
    np.zeros(3, dtype=np.float32),
    np.zeros(3, dtype=np.float64),
    np.zeros(3, dtype=int),
    np.zeros(3, dtype=bool),
    np.array([1.000000051, 0.000000051, -1.000000051], dtype=np.float32),
    np.array([1.000000051, 0.000000051, -1.000000051], dtype=np.float64),
    np.array([np.nan, np.nan, np.nan], dtype=np.float32),
    np.array([np.nan, np.nan, np.nan], dtype=np.float64),
    np.array([np.inf, -np.inf, np.nan], dtype=np.float32),
    np.array([np.inf, -np.inf, np.nan], dtype=np.float64),
    np.array([1, 2, 3, 4, 5, 6], dtype=np.float32)[::2],
    np.array([1, 2, 3, 4, 5, 6], dtype=np.float64)[::2],
    np.array([1, 2, 3, 4, 5, 6], dtype=np.float32)[::-2],
    np.array([1, 2, 3, 4, 5, 6], dtype=np.float64)[::-2],
))
def test_coords_add_vec(positions, vector):
    ref = positions.copy()
    ref += vector
    coords_add_vec(positions, vector)
    assert_equal(positions, ref)


def test_coords_add_vec_wrong_types():
    # wrong pos dtypes:
    with pytest.raises(ValueError):
        coords_add_vec(np.zeros((1, 3), dtype=np.float64),
                       np.zeros(3, dtype=np.float32))
    with pytest.raises(ValueError):
        coords_add_vec(np.zeros((1, 3), dtype=np.int32),
                       np.zeros(3, dtype=np.float64))
    # wrong pos type:
    with pytest.raises(TypeError):
        coords_add_vec([[1, 2, 3]], np.zeros(3, dtype=np.float64))
    # wrong vec dtype:
    with pytest.raises(ValueError):
        coords_add_vec(np.zeros((1, 3), dtype=np.float32),
                       np.array(["x", "y", "z"], dtype='S8'))
    # wrong vec type:
    with pytest.raises(TypeError):
        coords_add_vec(np.zeros((1, 3), dtype=np.float32), [0, 0, 0])


@pytest.mark.parametrize('positions', (
    np.array([[]], dtype=np.float32),
    np.zeros(3, dtype=np.float32),
    np.zeros((2, 2), dtype=np.float32),
    np.zeros((1, 1, 3), dtype=np.float32),
    np.zeros((1, 4), dtype=np.float32)
))
def test_coords_add_vec_wrong_pos_shape(positions):
    with pytest.raises(ValueError):
        coords_add_vec(positions, np.zeros(3, dtype=np.float32))


@pytest.mark.parametrize('vector', (
    np.array([], dtype=np.float32),
    np.zeros((1, 3), dtype=np.float32),
    np.zeros((1, 2, 3), dtype=np.float32),
    np.zeros(1, dtype=np.float32),
    np.zeros(4, dtype=np.float32)
))
def test_coords_add_vec_wrong_vec_shape(vector):
    positions = np.zeros((1, 3), dtype=np.float32)
    with pytest.raises(ValueError):
        coords_add_vec(positions, vector)
    with pytest.raises(ValueError):
        coords_add_vec(positions, vector.astype(np.float64))


@pytest.mark.parametrize('values', (
    [],                  # empty array
    [-999],              # single value array
    [1, 1, 1, 1],        # all identical
    [2, 3, 5, 7],        # all different, monotonic
    [5, 2, 7, 3],        # all different, non-monotonic
    [1, 2, 2, 4, 4, 6],  # duplicates, monotonic
    [1, 2, 2, 6, 4, 4],  # duplicates, non-monotonic
    [4, 2, 6, 1, 4, 2]   # duplicates, scrambled
))
@pytest.mark.parametrize('counts', (True, False))
@pytest.mark.parametrize('masks', (True, False))
def test_unique_int_1d(values, counts, masks):
    array = np.array(values, dtype=np.intp)
    ref = np.unique(array)
    res = unique_int_1d(array, return_counts=counts, return_masks=masks)
    if counts or masks:
        res = res[0]
    assert_equal(res, ref)
    assert type(res) == type(ref)
    assert res.dtype == ref.dtype


@pytest.mark.parametrize('values', (
    [],                  # empty array
    [-999],              # single value array
    [1, 1, 1, 1],        # all identical
    [2, 3, 5, 7],        # all different, monotonic
    [5, 2, 7, 3],        # all different, non-monotonic
    [1, 2, 2, 4, 4, 6],  # duplicates, monotonic
    [1, 2, 2, 6, 4, 4],  # duplicates, non-monotonic
    [4, 2, 6, 1, 4, 2]   # duplicates, scrambled
))
@pytest.mark.parametrize('masks', (True, False))
def test_unique_int_1d_return_counts(values, masks):
    array = np.array(values, dtype=np.intp)
    _, ref_counts = np.unique(array, return_counts=True)
    res = unique_int_1d(array, return_counts=True, return_masks=masks)
    counts = res[1]
    assert_equal(counts, ref_counts)
    assert type(counts) == type(ref_counts)
    assert counts.dtype == ref_counts.dtype


@pytest.mark.parametrize('values', (
    [],                  # empty array
    [-999],              # single value array
    [1, 1, 1, 1],        # all identical
    [2, 3, 5, 7],        # all different, monotonic
    [5, 2, 7, 3],        # all different, non-monotonic
    [1, 2, 2, 4, 4, 6],  # duplicates, monotonic
    [1, 2, 2, 6, 4, 4],  # duplicates, non-monotonic
    [4, 2, 6, 1, 4, 2]   # duplicates, scrambled
))
@pytest.mark.parametrize('counts', (True, False))
def test_unique_int_1d_return_masks(values, counts):
    array = np.array(values, dtype=np.intp)
    ref_masks = unique_masks_int_1d(array)
    res = unique_int_1d(array, return_counts=counts, return_masks=True)
    if counts:
        masks = res[2]
    else:
        masks = res[1]
    assert type(masks) == type(ref_masks)
    assert masks.dtype == ref_masks.dtype
    assert len(masks) == len(ref_masks)
    for i in range(len(masks)):
        assert isinstance(masks[i], (tuple, slice))
        assert_equal(masks[i], ref_masks[i])


@pytest.mark.parametrize('values', (
    [],                  # empty array
    [-999],              # single value array
    [1, 1, 1, 1],        # all identical
    [2, 3, 5, 7],        # all different, monotonic
    [5, 2, 7, 3],        # all different, non-monotonic
    [1, 2, 2, 4, 4, 6],  # duplicates, monotonic
    [1, 2, 2, 6, 4, 4],  # duplicates, non-monotonic
    [4, 2, 6, 1, 4, 2]   # duplicates, scrambled
))
def test_unique_masks_int_1d(values):
    array = np.array(values, dtype=np.intp)
    ismonotonic = len(array) > 0
    if len(array) > 1:
        ismonotonic = np.all((array[1:] - array[:-1]) >= 0)
    ref_unique = np.unique(values)
    ref_masks = np.empty(len(ref_unique), dtype=object)
    for i in range(len(ref_unique)):
        ref_masks[i] = (array == ref_unique[i]).nonzero()
        if ismonotonic:
            ref_masks[i] = slice(ref_masks[i][0][0], ref_masks[i][0][-1] + 1, 1)
    masks = unique_masks_int_1d(array)
    assert type(masks) == type(ref_masks)
    assert masks.dtype == ref_masks.dtype
    assert len(masks) == len(ref_masks)
    for i in range(len(masks)):
        if ismonotonic:
            assert isinstance(masks[i], slice)
        else:
            assert isinstance(masks[i], tuple)
            assert len(masks[i]) == 1
            assert isinstance(masks[i][0], np.ndarray)
            assert masks[i][0].dtype == np.intp
        assert_equal(masks[i], ref_masks[i])

@pytest.mark.parametrize('values, ref', (
    ([-3], True),                   # length-1 array, negative
    ([3], True),                    # length-1 array, positive
    ([0], True),                    # length-1 array, zero
    ([-5, -4], True),               # length-2 contiguous range
    ([1, 2, 3, 4, 5], True),        # contiguous range, positive
    ([-5, -4, -3, -2, -1], True),   # contiguous range, negative
    ([-2, -1, 0, 1, 2], True),      # contiguous range, neg to pos
    ([], False),                    # empty array
    ([0, 0], False),                # length-2 array, zeros
    ([-3, -4], False),              # length-2 inverted range
    ([5, 4, 3, 2, 1], False),       # inverted range, positive
    ([-1, -2, -3, -4, -5], False),  # inverted range, negative
    ([2, 1, 0, -1, -2], False),     # inverted range, pos to neg
    ([1, 3, 5, 7, 9], False),       # strided range, positive
    ([-9, -7, -5, -3, -1], False),  # strided range, negative
    ([-5, -3, -1, 1, 3], False),    # strided range, neg to pos
    ([3, 1, -1, -3, -5], False),    # inverted strided range, pos to neg
    ([1, 1, 1, 1], False),          # all identical
    ([2, 3, 5, 7], False),          # monotonic
    ([5, 2, 7, 3], False),          # non-monotonic
    ([1, 2, 2, 3, 4], False),       # range with middle duplicates
    ([1, 2, 3, 4, 4], False),       # range with end duplicates
    ([1, 1, 2, 3, 4], False),       # range with start duplicates
    ([-1, 2, 2, 4, 3], False)       # duplicates, non-monotonic
))
def test_iscontiguous_int_1d(values, ref):
    array = np.array(values, dtype=np.intp)
    res = iscontiguous_int_1d(array)
    assert_equal(res, ref)
    assert type(res) == bool


@pytest.mark.parametrize('arr', (
    [],                       # empty array
    [1],                      # single value array
    [1, 1, 1, 1],             # all identical
    [0, 3, 5, 7],             # all different, monotonic
    [5, 2, 7, 3],             # all different, non-monotonic
    [-1, -1, 2, 2, 4, 4, 6],  # duplicates, monotonic
    [1, 2, 2, 6, 4, 4, -1],   # duplicates, non-monotonic
    [4, 2, 6, 1, 4, 2, -1]    # duplicates, scrambled
))
@pytest.mark.parametrize('value', (-1, 0, 1, 2, 3, 4, 5, 6, 7))
def test_argwhere_int_1d(arr, value):
    arr = np.array(arr, dtype=np.intp)
    ref = np.argwhere(arr == value).ravel()
    res = argwhere_int_1d(arr, value)
    assert_equal(res, ref)
    assert type(res) == type(ref)
    assert res.dtype == ref.dtype


@pytest.mark.parametrize('indices, ref', (
    ([5], slice(5, 6, 1)),               # single value array
    ([0, 1], slice(0, 2, 1)),            # two value array, stride 1
    ([2, 12], slice(2, 13, 10)),         # two value array, stride 10
    ([1, 0], slice(1, None, -1)),        # two value array, stride -1
    ([12, 1], slice(12, 0, -11)),        # two value array, stride -11
    ([0, 1, 2, 3], slice(0, 4, 1)),      # monotonic, stride 1
    ([0, 2, 4, 6], slice(0, 7, 2)),      # monotonic, stride 2
    ([2, 3, 4, 5], slice(2, 6, 1)),      # monotonic with offset, stride 1
    ([3, 5, 7, 9], slice(3, 10, 2)),     # monotonic with offset, stride 2
    ([3, 2, 1, 0], slice(3, None, -1)),  # monotonic, stride -1
    ([9, 6, 3, 0], slice(9, None, -3)),  # monotonic, stride -3
    ([7, 6, 5, 4], slice(7, 3, -1)),     # monotonic with offset, stride -1
    ([14, 10, 6], slice(14, 5, -4))      # monotonic with offset, stride -4
))
def test_indices_to_slice_1d_slice(indices, ref):
    ix = np.array(indices, dtype=np.intp)
    res = indices_to_slice_1d(ix)
    assert isinstance(res, slice)
    assert res == ref
    values = np.arange(60, dtype=np.float32).reshape((-1, 3))
    fancy_indexed = values[ix]
    view = values[res]
    assert not view.flags['OWNDATA']  # the purpose of indices_to_slice_1d
    assert_equal(view, fancy_indexed)


@pytest.mark.parametrize('indices', (
    [],                       # empty array
    [1, 1, 1, 1],             # all identical
    [0, 3, 5, 7],             # all different, monotonic, non-uniform stride
    [5, 2, 7, 3],             # all different, non-monotonic
    [0, 0, 2, 2, 4, 4, 6],    # duplicates, monotonic
    [1, 2, 2, 6, 4, 4, 1],    # duplicates, non-monotonic
    [4, 2, 6, 1, 4, 2, 1],    # duplicates, scrambled
    [-1, 0, 1, 2, 3, 4, 5],   # monotonic, stride 1, negative values
    [1, 0, -1, -2, -3, -4],   # monotonic, stride -1, negative values
    [-3, -1, 1, 3, 5, 7],     # monotonic, stride 2, negative values
    [3, 0, -3, -6, -9, -12]   # monotonic, stride -3, negative values
))
def test_indices_to_slice_1d_noslice(indices):
    ix = np.array(indices, dtype=np.intp)
    assert indices_to_slice_1d(ix) is ix


@pytest.mark.parametrize('edges,ref', [
    ([[0, 1], [1, 2], [2, 3], [3, 4]],
     [[0, 1, 2, 3, 4]]),  # linear chain
    ([[0, 1], [1, 2], [2, 3], [3, 4], [4, 10]],
     [[0, 1, 2, 3, 4]]),  # unused edge (4, 10)
    ([[0, 1], [1, 2], [2, 3]],
     [[0, 1, 2, 3], [4]]),  # lone atom
    ([[0, 1], [1, 2], [2, 0], [3, 4], [4, 3]],
     [[0, 1, 2], [3, 4]]),  # circular
])
def test_find_fragments(edges, ref):
    atoms = np.arange(5)

    fragments = find_fragments(atoms, edges)

    assert len(fragments) == len(ref)
    for frag, r in zip(fragments, ref):
        assert_equal(frag, r)
