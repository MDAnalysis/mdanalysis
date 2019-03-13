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

from MDAnalysis.lib._cutil import (unique_int_1d, unique_masks_int_1d,
                                   iscontiguous_int_1d, argwhere_int_1d,
                                   find_fragments)


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
