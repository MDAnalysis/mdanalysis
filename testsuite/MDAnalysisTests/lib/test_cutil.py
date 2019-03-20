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
from numpy.testing import assert_equal, assert_almost_equal

from MDAnalysis.lib._cutil import (coords_add_vector, _coords_add_vectors,
                                   coords_center, unique_int_1d,
                                   unique_masks_int_1d, iscontiguous_int_1d,
                                   argwhere_int_1d, indices_to_slice_1d,
                                   find_fragments)



class TestCoordFunctions(object):

    precision = 14
    seed = 1337
    rng = None

    @classmethod
    def setup_class(self):
        """Setup class-local PRNG"""
        self.rng = np.random.RandomState(seed=self.seed)

    def gen_masks(self, n_coords, duplicates=False, random=False):
        ix = np.arange(n_coords, dtype=np.intp)
        if n_coords > 1 and duplicates:
            ix[:-1:2] = ix[1::2]
        if random:
            self.rng.shuffle(ix)
        print(ix)
        masks = unique_masks_int_1d(ix)
        return masks


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
    def test_coords_add_vector(self, positions, vector):
        ref = positions.copy()
        ref += vector
        coords_add_vector(positions, vector)
        assert_equal(positions, ref)


    def test_coords_add_vector_wrong_types(self):
        # wrong pos dtypes:
        with pytest.raises(ValueError):
            coords_add_vector(np.zeros((1, 3), dtype=np.float64),
                              np.zeros(3, dtype=np.float32))
        with pytest.raises(ValueError):
            coords_add_vector(np.zeros((1, 3), dtype=np.int32),
                              np.zeros(3, dtype=np.float64))
        # wrong pos type:
        with pytest.raises(TypeError):
            coords_add_vector([[1, 2, 3]], np.zeros(3, dtype=np.float64))
        # wrong vec dtype:
        with pytest.raises(ValueError):
            coords_add_vector(np.zeros((1, 3), dtype=np.float32),
                              np.array(["x", "y", "z"], dtype='S8'))
        # wrong vec type:
        with pytest.raises(TypeError):
            coords_add_vector(np.zeros((1, 3), dtype=np.float32), [0, 0, 0])


    @pytest.mark.parametrize('positions', (
        np.array([[]], dtype=np.float32),
        np.zeros(3, dtype=np.float32),
        np.zeros((2, 2), dtype=np.float32),
        np.zeros((1, 1, 3), dtype=np.float32),
        np.zeros((1, 4), dtype=np.float32)
    ))
    def test_coords_add_vector_wrong_pos_shape(self, positions):
        with pytest.raises(ValueError):
            coords_add_vector(positions, np.zeros(3, dtype=np.float32))


    @pytest.mark.parametrize('vector', (
        np.array([], dtype=np.float32),
        np.zeros((1, 3), dtype=np.float32),
        np.zeros((1, 2, 3), dtype=np.float32),
        np.zeros(1, dtype=np.float32),
        np.zeros(4, dtype=np.float32)
    ))
    def test_coords_add_vector_wrong_vec_shape(self, vector):
        positions = np.zeros((1, 3), dtype=np.float32)
        with pytest.raises(ValueError):
            coords_add_vector(positions, vector)
        with pytest.raises(ValueError):
            coords_add_vector(positions, vector.astype(np.float64))


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
    def test_coords_add_vectors_nomasks(self, positions, vector):
        vectors = np.repeat(vector[None, :], len(positions), axis=0)
        ref = positions.copy()
        ref += vectors
        _coords_add_vectors(positions, vectors, compound_masks=None)
        assert_equal(positions, ref)


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
    @pytest.mark.parametrize('dup', (False, True))
    @pytest.mark.parametrize('rnd', (False, True))
    def test_coords_add_vectors_masks(self, positions, vector, dup, rnd):
        comp_masks = self.gen_masks(len(positions), duplicates=dup, random=rnd)
        vectors = np.repeat(vector[None, :], len(comp_masks), axis=0)
        ref = positions.copy()
        for i, mask in enumerate(comp_masks):
            if isinstance(mask, slice):
                ref[mask] += vectors[i]
            else:
                print(np.asarray(mask))
                for j in mask:
                    ref[j] += vectors[i]
        _coords_add_vectors(positions, vectors, compound_masks=comp_masks)
        assert_equal(positions, ref)


    @pytest.mark.parametrize('positions', (
        np.array([[]], dtype=np.float32),
        np.zeros(3, dtype=np.float32),
        np.zeros((2, 2), dtype=np.float32),
        np.zeros((1, 1, 3), dtype=np.float32),
        np.zeros((1, 4), dtype=np.float32)
    ))
    @pytest.mark.parametrize('masks', (None, np.empty(2, dtype=object)))
    def test_coords_add_vectors_wrong_pos_shape(self, positions, masks):
        if masks is not None:
            vectors = np.zeros((len(masks), 3), dtype=np.float32)
        else:
            vectors = np.zeros((len(positions), 3), dtype=np.float32)
        with pytest.raises(ValueError):
            _coords_add_vectors(positions, vectors, compound_masks=masks)


    @pytest.mark.parametrize('vectors', (
        np.array([[]], dtype=np.float32),
        np.zeros((3, 1), dtype=np.float32),
        np.zeros((1, 2, 3), dtype=np.float32),
        np.zeros((0, 5), dtype=np.float32),
        np.zeros(4, dtype=np.float32)
    ))
    def test_coords_add_vectors_wrong_vec_shape(self, vectors):
        positions = np.zeros((len(vectors), 3), dtype=np.float32)
        with pytest.raises(ValueError):
            _coords_add_vectors(positions, vectors, compound_masks=None)
        masks = np.empty(len(vectors), dtype=object)
        with pytest.raises(ValueError):
            _coords_add_vectors(positions, vectors, compound_masks=masks)


    def test_coords_add_vectors_wrong_masks_shape(self):
        positions = np.zeros((2, 3), dtype=np.float32)
        vectors = np.zeros((2, 3), dtype=np.float32)
        masks = np.empty(1, dtype=object)
        with pytest.raises(ValueError):
            _coords_add_vectors(positions, vectors, compound_masks=masks)


    @pytest.mark.parametrize('coords', (
        np.empty((0, 3), dtype=np.float32),
        np.array([[1, 2, 3]], dtype=np.float32),
        np.reshape(np.arange(-30, 30, dtype=np.float32), (-1, 3))
    ))
    @pytest.mark.parametrize('check_weights', (False, True))
    @pytest.mark.parametrize('ret_masks', (False, True))
    def test_coords_center_noweights_nocomp(self, coords, check_weights,
                                            ret_masks):
        if len(coords) < 2:
            ref = coords.astype(np.float64)
        else:
            ref = coords.mean(axis=0, dtype=np.float64).reshape(-1, 3)
        res = coords_center(coords, weights=None, compound_indices=None,
                            check_weights=check_weights,
                            return_compound_masks=ret_masks)
        assert res.dtype == np.float64
        assert res.shape == (int(len(coords) != 0), 3)
        assert_equal(res, ref)


    @pytest.mark.parametrize('coords', (
        np.empty((0, 3), dtype=np.float32),
        np.array([[1, 2, 3]], dtype=np.float32),
        np.reshape(np.arange(-30, 30, dtype=np.float32), (-1, 3))
    ))
    @pytest.mark.parametrize('check_weights', (False, True))
    @pytest.mark.parametrize('ret_masks', (False, True))
    def test_coords_center_weights_nocomp(self, coords, check_weights,
                                          ret_masks):
        weights = np.arange(1, len(coords) + 1, dtype=np.float64)
        if len(coords) == 0:
            ref = coords.astype(np.float64)
        else:
            ref = np.average(coords, axis=0, weights=weights).reshape((1, 3))
        res = coords_center(coords, weights=weights, compound_indices=None,
                            check_weights=check_weights,
                            return_compound_masks=ret_masks)
        assert res.dtype == np.float64
        assert res.shape == (int(len(coords) != 0), 3)
        assert_almost_equal(res, ref, decimal=self.precision)


    @pytest.mark.parametrize('coords, comp_ix', (
        (np.empty((0, 3), dtype=np.float32),
         np.empty(0, dtype=np.intp)),
        (np.array([[1, 2, 3]], dtype=np.float32),
         np.array([-999], dtype=np.intp)),
        (np.reshape(np.arange(-10, 20, dtype=np.float32), (-1, 3)),
         np.array([0, 0, 0, 1, 2, 2, 2, 3, 3, 4], dtype=np.intp)),
        (np.reshape(np.arange(-20, 10, dtype=np.float32), (-1, 3)),
         np.array([4, -2, 0, 1, -2, 0, 3, 3, 0, -2], dtype=np.intp)),
    ))
    @pytest.mark.parametrize('check_weights', (False, True))
    @pytest.mark.parametrize('ret_masks', (False, True))
    def test_coords_center_noweights_comp(self, coords, comp_ix, check_weights,
                                          ret_masks):
        unique_comp_ix = np.unique(comp_ix)
        if len(coords) < 2:
            ref = coords.astype(np.float64)
        else:
            ref = np.empty((len(unique_comp_ix), 3), dtype=np.float64)
            for i, uix in enumerate(unique_comp_ix):
                mask = (comp_ix == uix).nonzero()
                ref[i] = coords[mask].mean(axis=0, dtype=np.float64)
        res = coords_center(coords, weights=None, compound_indices=comp_ix,
                            check_weights=check_weights,
                            return_compound_masks=ret_masks)
        if ret_masks:
            ref_masks = unique_masks_int_1d(comp_ix, assume_unsorted=True)
            res, res_masks = res
            assert res_masks.dtype == object
            assert len(res_masks) == len(ref_masks)
            for i in range(len(res_masks)):
                assert_equal(np.asarray(res_masks[i]), np.asarray(ref_masks[i]))
        assert res.dtype == np.float64
        assert res.shape == (len(unique_comp_ix), 3)
        assert_almost_equal(res, ref, decimal=self.precision)


    @pytest.mark.parametrize('coords, weights, comp_ix', (
        (np.empty((0, 3), dtype=np.float32),
         np.empty(0, dtype=np.float64),
         np.empty(0, dtype=np.intp)),
        (np.array([[1, 2, 3]], dtype=np.float32),
         np.array([2], dtype=np.float64),
         np.array([-999], dtype=np.intp)),
        (np.reshape(np.arange(-10, 20, dtype=np.float32), (-1, 3)),
         np.array([1, 1, 1, 2, 3, 3, 3, 4, 4, 5], dtype=np.float64),
         np.array([0, 0, 0, 1, 2, 2, 2, 3, 3, 4], dtype=np.intp)),
        (np.reshape(np.arange(-20, 10, dtype=np.float32), (-1, 3)),
         np.array([4, -1, 1, 1, -2, 1, 2, 3, 5, -2], dtype=np.float64),
         np.array([4, -2, 0, 1, -2, 0, 3, 3, 0, -2], dtype=np.intp)),
    ))
    @pytest.mark.parametrize('check_weights', (False, True))
    @pytest.mark.parametrize('ret_masks', (False, True))
    def test_coords_center_weights_comp(self, coords, weights, comp_ix,
                                        check_weights, ret_masks):
        unique_comp_ix = np.unique(comp_ix)
        if len(coords) == 0:
            ref = coords.astype(np.float64)
        else:
            ref = np.empty((len(unique_comp_ix), 3), dtype=np.float64)
            for i, uix in enumerate(unique_comp_ix):
                mask = (comp_ix == uix).nonzero()
                ref[i] = np.average(coords[mask], axis=0, weights=weights[mask])
        res = coords_center(coords, weights=weights, compound_indices=comp_ix,
                            check_weights=check_weights,
                            return_compound_masks=ret_masks)
        if ret_masks:
            ref_masks = unique_masks_int_1d(comp_ix, assume_unsorted=True)
            res, res_masks = res
            assert res_masks.dtype == object
            assert len(res_masks) == len(ref_masks)
            for i in range(len(res_masks)):
                assert_equal(np.asarray(res_masks[i]), np.asarray(ref_masks[i]))
        assert res.dtype == np.float64
        assert res.shape == (len(unique_comp_ix), 3)
        assert_almost_equal(res, ref, decimal=self.precision)


    @pytest.mark.parametrize('coords', (
        np.empty((0, 3), dtype=np.float32),
        np.array([[1, 2, 3]], dtype=np.float32),
        np.reshape(np.arange(-30, 30, dtype=np.float32), (-1, 3))
    ))
    @pytest.mark.parametrize('check_weights', (False, True))
    @pytest.mark.parametrize('ret_masks', (False, True))
    def test_coords_center_zero_weights_nocomp(self, coords, check_weights,
                                               ret_masks):
        weights = np.zeros(len(coords), dtype=np.float64)
        if check_weights and len(coords) > 0:
            with pytest.raises(ValueError):
                res = coords_center(coords, weights=weights,
                                    compound_indices=None,
                                    check_weights=check_weights,
                                    return_compound_masks=ret_masks)
        else:
            if len(coords) == 0:
                ref = coords.astype(np.float64)
            else:
                ref = np.full((1, 3), np.nan, dtype=np.float64)
            res = coords_center(coords, weights=weights, compound_indices=None,
                                check_weights=check_weights,
                                return_compound_masks=ret_masks)
            assert res.dtype == np.float64
            assert res.shape == (int(len(coords) != 0), 3)
            assert_equal(res, ref)


    @pytest.mark.parametrize('coords, weights, comp_ix', (
        (np.empty((0, 3), dtype=np.float32),
         np.empty(0, dtype=np.float64),
         np.empty(0, dtype=np.intp)),
        (np.array([[1, 2, 3]], dtype=np.float32),
         np.array([0], dtype=np.float64),
         np.array([-999], dtype=np.intp)),
        (np.reshape(np.arange(-10, 20, dtype=np.float32), (-1, 3)),
         np.array([1, 1, 1, 2, -1, 2, -1, 1, 1, 5], dtype=np.float64),
         np.array([0, 0, 0, 1, 2, 2, 2, 3, 3, 4], dtype=np.intp)),
        (np.reshape(np.arange(-20, 10, dtype=np.float32), (-1, 3)),
         np.array([4, -1, 1, 1, -2, 1, 2, -2, 5, -2], dtype=np.float64),
         np.array([4, -2, 0, 1, -2, 0, 3, 3, 0, -2], dtype=np.intp)),
    ))
    @pytest.mark.parametrize('check_weights', (False, True))
    @pytest.mark.parametrize('ret_masks', (False, True))
    def test_coords_center_zero_weights_comp(self, coords, weights, comp_ix,
                                             check_weights, ret_masks):
        if check_weights and len(coords) > 0:
            with pytest.raises(ValueError):
                res = coords_center(coords, weights=weights,
                                    compound_indices=comp_ix,
                                    check_weights=check_weights,
                                    return_compound_masks=ret_masks)
        else:
            unique_comp_ix = np.unique(comp_ix)
            if len(coords) == 0:
                ref = coords.astype(np.float64)
            else:
                ref = np.empty((len(unique_comp_ix), 3), dtype=np.float64)
                for i, uix in enumerate(unique_comp_ix):
                    m = (comp_ix == uix).nonzero()
                    inv_sum_weights = 1.0 / weights[m].sum()
                    if not np.isfinite(inv_sum_weights):
                        ref[i] = np.nan
                    else:
                        ref[i] = (coords[m] * weights[m][:, None]).sum(axis=0)
                        ref[i] *= inv_sum_weights
                res = coords_center(coords, weights=weights,
                                    compound_indices=comp_ix,
                                    check_weights=check_weights,
                                    return_compound_masks=ret_masks)
                if ret_masks:
                    ref_masks = unique_masks_int_1d(comp_ix,
                                                    assume_unsorted=True)
                    res, res_masks = res
                    assert res_masks.dtype == object
                    assert len(res_masks) == len(ref_masks)
                    for i in range(len(res_masks)):
                        assert_equal(np.asarray(res_masks[i]),
                                     np.asarray(ref_masks[i]))
                assert res.dtype == np.float64
                assert res.shape == (len(unique_comp_ix), 3)
                assert_almost_equal(res, ref, decimal=self.precision)


    @pytest.mark.parametrize('coords, weights, comp_ix', (
        (np.empty((0, 2), dtype=np.float32),
         np.empty(0, dtype=np.float64),
         np.empty(0, dtype=np.intp)),
        (np.empty((0, 2), dtype=np.float32),
         None,
         np.empty(0, dtype=np.intp)),
        (np.empty((0, 2), dtype=np.float32),
         None,
         None),
        (np.array([[1, 2, 3]], dtype=np.float32),
         np.array([1, 1], dtype=np.float64),
         np.array([0], dtype=np.intp)),
        (np.array([[1, 2, 3]], dtype=np.float32),
         None,
         np.array([0, 0], dtype=np.intp)),
        (np.array([[1, 2, 3]], dtype=np.float32),
         np.array([1, 1], dtype=np.float64),
         None),
        (np.array([[1, 2, 3]], dtype=np.float32),
         np.array([1, 1], dtype=np.float64),
         np.array([0, 0], dtype=np.intp)),
    ))
    @pytest.mark.parametrize('check_weights', (False, True))
    @pytest.mark.parametrize('ret_masks', (False, True))
    def test_coords_center_wrong_shapes(self, coords, weights, comp_ix,
                                        check_weights, ret_masks):
        with pytest.raises(ValueError):
            res = coords_center(coords, weights=weights,
                                compound_indices=comp_ix,
                                check_weights=check_weights,
                                return_compound_masks=ret_masks)



class TestInt1dFunctions(object):

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
    @pytest.mark.parametrize('unsorted', (True, False))
    def test_unique_int_1d(self, values, counts, masks, unsorted):
        array = np.array(values, dtype=np.intp)
        ref = np.unique(array)
        res = unique_int_1d(array, return_counts=counts, return_masks=masks,
                            assume_unsorted=unsorted)
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
    @pytest.mark.parametrize('unsorted', (True, False))
    def test_unique_int_1d_return_counts(self, values, masks, unsorted):
        array = np.array(values, dtype=np.intp)
        _, ref_counts = np.unique(array, return_counts=True)
        res = unique_int_1d(array, return_counts=True, return_masks=masks,
                            assume_unsorted=unsorted)
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
    @pytest.mark.parametrize('unsorted', (True, False))
    def test_unique_int_1d_return_masks(self, values, counts, unsorted):
        array = np.array(values, dtype=np.intp)
        ref_masks = unique_masks_int_1d(array, assume_unsorted=unsorted)
        res = unique_int_1d(array, return_counts=counts, return_masks=True,
                            assume_unsorted=unsorted)
        if counts:
            masks = res[2]
        else:
            masks = res[1]
        assert type(masks) == type(ref_masks)
        assert masks.dtype == ref_masks.dtype
        assert len(masks) == len(ref_masks)
        for i in range(len(masks)):
            assert isinstance(masks[i], slice) or \
                masks[i].__class__.__name__ == '_memoryviewslice'
            if isinstance(masks[i], slice):
                assert_equal(masks[i], ref_masks[i])
            else:
                assert_equal(np.asarray(masks[i]), np.asarray(ref_masks[i]))


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
    @pytest.mark.parametrize('unsorted', (True, False))
    def test_unique_masks_int_1d(self, values, unsorted):
        array = np.array(values, dtype=np.intp)
        ismonotonic = len(array) > 0
        if len(array) > 1:
            ismonotonic = np.all((array[1:] - array[:-1]) >= 0)
        ref_unique = np.unique(values)
        ref_masks = np.empty(len(ref_unique), dtype=object)
        for i in range(len(ref_unique)):
            ref_masks[i] = (array == ref_unique[i]).nonzero()[0]
            if ismonotonic and (not unsorted or len(array) == 1):
                ref_masks[i] = slice(ref_masks[i][0], ref_masks[i][-1] + 1, 1)
        masks = unique_masks_int_1d(array, assume_unsorted=unsorted)
        assert type(masks) == type(ref_masks)
        assert masks.dtype == ref_masks.dtype
        assert len(masks) == len(ref_masks)
        for i in range(len(masks)):
            if ismonotonic and (not unsorted or len(array) == 1):
                assert isinstance(masks[i], slice)
                assert_equal(masks[i], ref_masks[i])
            else:
                assert masks[i].__class__.__name__ == '_memoryviewslice'
                assert len(masks[i]) == len(ref_masks[i])
                assert np.asarray(masks[i]).dtype == np.intp
                assert_equal(np.asarray(masks[i]), ref_masks[i])

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
    def test_iscontiguous_int_1d(self, values, ref):
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
    def test_argwhere_int_1d(self, arr, value):
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
    def test_indices_to_slice_1d_slice(self, indices, ref):
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
    def test_indices_to_slice_1d_noslice(self, indices):
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
