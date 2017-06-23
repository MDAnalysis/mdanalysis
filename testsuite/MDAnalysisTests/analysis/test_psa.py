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

import MDAnalysis as mda
import MDAnalysis.analysis.psa as PSA

from numpy.testing import (TestCase, dec, assert_array_less,
                           assert_array_almost_equal, assert_,
                           assert_almost_equal, assert_equal)
import numpy as np
import scipy
import scipy.spatial

from MDAnalysisTests.datafiles import PSF, DCD, DCD2
from MDAnalysisTests import parser_not_found, tempdir, module_not_found


class TestPSAnalysis(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.iu1 = np.triu_indices(3, k=1)
        self.universe1 = mda.Universe(PSF, DCD)
        self.universe2 = mda.Universe(PSF, DCD2)
        self.universe_rev = mda.Universe(PSF, DCD)
        self.universes = [self.universe1, self.universe2, self.universe_rev]
        self.psa = PSA.PSAnalysis(self.universes,
                                  path_select='name CA',
                                  targetdir=self.tmpdir.name)

        self.psa.generate_paths(align=True)
        self.psa.paths[-1] = self.psa.paths[-1][::-1,:,:] # reverse third path
        self._run()
        self._plot()

    def _run(self):
        self.psa.run(metric='hausdorff')
        self.hausd_matrix = self.psa.get_pairwise_distances()
        self.psa.run(metric='discrete_frechet')
        self.frech_matrix = self.psa.get_pairwise_distances()
        self.hausd_dists = self.hausd_matrix[self.iu1]
        self.frech_dists = self.frech_matrix[self.iu1]

    def _plot(self):
        self.plot_data = self.psa.plot()

    def tearDown(self):
        del self.universe1
        del self.universe2
        del self.universe_rev
        del self.psa
        del self.tmpdir

    def test_hausdorff_bound(self):
        err_msg = "Some Frechet distances are smaller than corresponding "      \
                + "Hausdorff distances"
        assert_array_less(self.hausd_dists, self.frech_dists, err_msg)

    def test_reversal_hausdorff(self):
        err_msg = "Hausdorff distances changed after path reversal"
        assert_array_almost_equal(self.hausd_matrix[1,2],
                                  self.hausd_matrix[0,1],
                                  decimal=3, err_msg=err_msg)

    def test_reversal_frechet(self):
        err_msg = "Frechet distances did not increase after path reversal"
        assert_(self.frech_matrix[1,2] >= self.frech_matrix[0,1], err_msg)

    def test_dendrogram_produced(self):
        err_msg = "Dendrogram dictionary object was not produced"
        assert_(type(self.plot_data[1]) is dict, err_msg)

    def test_dist_mat_to_vec_i_less_j(self):
        """Test the index of corresponding distance vector is correct if i < j"""
        err_msg = "dist_mat_to_vec function returning wrong values"
        assert_equal(PSA.dist_mat_to_vec(5, 3, 4), 9, err_msg)

    def test_dist_mat_to_vec_i_greater_j(self):
        """Test the index of corresponding distance vector is correct if i > j"""
        err_msg = "dist_mat_to_vec function returning wrong values"
        assert_equal(PSA.dist_mat_to_vec(5, 4, 3), 9, err_msg)

    def test_dist_mat_to_vec_input_numpy_integer_32(self):
        """Test whether inputs are supported as numpy integers rather than normal Integers"""
        err_msg = "dist_mat_to_vec function returning wrong values"
        assert_equal(PSA.dist_mat_to_vec(np.int32(5), np.int32(3), np.int32(4)), np.int32(9), err_msg)

    def test_dist_mat_to_vec_input_numpy_integer_16(self):
        """Test whether inputs are supported as numpy integers rather than normal Integers"""
        err_msg = "dist_mat_to_vec function returning wrong values"
        assert_equal(PSA.dist_mat_to_vec(np.int16(5), np.int16(3), np.int16(4)), np.int16(9), err_msg)

class TestPSAExceptions(TestCase):
    '''Tests for exceptions that should be raised
    or caught by code in the psa module.'''

    def test_get_path_metric_func_bad_key(self):
        '''Test that KeyError is caught by
        get_path_metric_func().'''

        try:
           PSA.get_path_metric_func('123456')
        except KeyError:
            self.fail('KeyError should be caught')

    def test_get_coord_axes_bad_dims(self):
        """Test that ValueError is raised when
        numpy array with incorrect dimensions
        is fed to get_coord_axes()."""

        with self.assertRaises(ValueError):
            PSA.get_coord_axes(np.zeros((5,5,5,5)))

    def test_dist_mat_to_vec_func_out_of_bounds(self):
        """Test that ValueError is raised when i or j or both are
        out of bounds of N"""

        # Check if i is out of bounds of N
        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(5, 6, 4)

        # Check if j is out of bounds of N
        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(5, 4, 6)

        # Check if both i and j are out of bounds of N
        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(5, 6, 7)

        # Check if i is negative
        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(5, -1, 2)

        # Check if j is negative
        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(5, 1, -2)

        # Check if N is less than 2
        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(1, 0, 0)

    def test_dist_mat_to_vec_func_i_equals_j(self):
        """Test that ValueError is raised when i == j or i,j == N"""

        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(5, 4, 4)

        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(4, 6, 4)

    def test_dist_mat_to_vec_func_bad_integers(self):
        """Test that ValueError is raised when i or j are
        not Integers"""

        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(5, '6', '7')

        with self.assertRaises(ValueError):
            PSA.dist_mat_to_vec(5, float(6), 7)


class _BaseHausdorffDistance(TestCase):
    '''Base Class setup and unit tests
    for various Hausdorff distance
    calculation properties.'''

    def setUp(self):
        self.random_angles = np.random.random((100,)) * np.pi * 2
        self.random_columns = np.column_stack((self.random_angles,
                                               self.random_angles,
                                               np.zeros((100,))))
        self.random_columns[...,0] = np.cos(self.random_columns[...,0])
        self.random_columns[...,1] = np.sin(self.random_columns[...,1])
        self.random_columns_2 = np.column_stack((self.random_angles,
                                                 self.random_angles,
                                                 np.zeros((100,))))
        self.random_columns_2[1:,0] = np.cos(self.random_columns_2[1:,0]) * 2.0
        self.random_columns_2[1:,1] = np.sin(self.random_columns_2[1:,1]) * 2.0
        # move one point farther out so we don't have two perfect circles
        self.random_columns_2[0,0] = np.cos(self.random_columns_2[0,0]) * 3.3
        self.random_columns_2[0,1] = np.sin(self.random_columns_2[0,1]) * 3.3
        self.path_1 = self.random_columns
        self.path_2 = self.random_columns_2

    def tearDown(self):
        del self.random_angles
        del self.random_columns
        del self.random_columns_2
        del self.path_1
        del self.path_2

    def test_symmetry(self):
        '''Ensure that the undirected (symmetric)
        Hausdorff distance is actually symmetric
        for a given Hausdorff metric, h.'''
        forward = self.h(self.path_1, self.path_2)
        reverse = self.h(self.path_2, self.path_1)
        # lower precision on 32bit
        assert_almost_equal(forward, reverse, decimal=15)

    def test_hausdorff_value(self):
        '''Test that the undirected Hausdorff
        distance matches expected value for
        the simple case here.'''
        actual = self.h(self.path_1, self.path_2)
        # unless I pin down the random generator
        # seems unstable to use decimal > 2
        assert_almost_equal(actual, self.expected,
                            decimal=2)

class TestHausdorffSymmetric(_BaseHausdorffDistance):
    '''Tests for conventional and symmetric (undirected)
    Hausdorff distance between point sets in 3D.'''

    def setUp(self):
        super(TestHausdorffSymmetric, self).setUp()
        self.h = PSA.hausdorff
        # radii differ by ~ 2.3 for outlier
        self.expected = 2.3

class TestWeightedAvgHausdorffSymmetric(_BaseHausdorffDistance):
    '''Tests for weighted average and symmetric (undirected)
    Hausdorff distance between point sets in 3D.'''

    def setUp(self):
        super(TestWeightedAvgHausdorffSymmetric, self).setUp()
        self.h = PSA.hausdorff_wavg
        self.distance_matrix = scipy.spatial.distance.cdist(self.path_1,
                                                            self.path_2)
        self.expected = (np.mean(np.amin(self.distance_matrix, axis=0)) +
                         np.mean(np.amin(self.distance_matrix, axis = 1))) / 2.

    def test_asymmetric_weight(self):
        '''Test to ensure that increasing N points in one of the paths
        does NOT increase the weight of its contributions.'''
        inflated_path_1 = np.concatenate((self.path_1, self.path_1))
        inflated_path_2 = np.concatenate((self.path_2, self.path_2))
        d_inner_inflation = self.h(inflated_path_1, self.path_2)
        d_outer_inflation = self.h(self.path_1, inflated_path_2)
        assert_almost_equal(d_inner_inflation,
                            d_outer_inflation)

class TestAvgHausdorffSymmetric(_BaseHausdorffDistance):
    '''Tests for unweighted average and symmetric (undirected)
    Hausdorff distance between point sets in 3D.'''

    def setUp(self):
        super(TestAvgHausdorffSymmetric, self).setUp()
        self.h = PSA.hausdorff_avg
        self.distance_matrix = scipy.spatial.distance.cdist(self.path_1,
                                                            self.path_2)
        self.expected = np.mean(np.append(np.amin(self.distance_matrix, axis=0),
                                np.amin(self.distance_matrix, axis = 1)))

    def test_asymmetric_weight(self):
        '''Test to ensure that increasing N points in one of the paths
        increases the weight of its contributions.'''
        inflated_path_1 = np.concatenate((self.path_1, self.path_1))
        inflated_path_2 = np.concatenate((self.path_2, self.path_2))
        d_inner_inflation = self.h(inflated_path_1, self.path_2)
        d_outer_inflation = self.h(self.path_1, inflated_path_2)
        assert_array_less(d_inner_inflation,
                          d_outer_inflation)

class DiscreteFrechetDistance(TestCase):
    # unit tests for the discrete Frechet distance

    def setUp(self):
        np.random.seed(50)
        random_angles = np.random.random((100,)) * np.pi * 2
        random_columns = np.column_stack((random_angles, random_angles,
                                          np.zeros((100,))))
        random_columns[...,0] = np.cos(random_columns[...,0])
        random_columns[...,1] = np.sin(random_columns[...,1])
        random_columns_2 = np.column_stack((random_angles, random_angles,
                                            np.zeros((100,))))
        random_columns_2[...,0] = np.cos(random_columns_2[...,0]) * 5.5
        random_columns_2[...,1] = np.sin(random_columns_2[...,1]) * 5.5
        self.path_1 = random_columns
        self.path_2 = random_columns_2

    def tearDown(self):
        del self.path_1
        del self.path_2

    def test_discrete_Frechet_concentric_circles(self):
        # test for the simple case of the discrete Frechet distance
        # between concentric circular paths, which for a sufficiently
        # high random discrete point density around each circle
        # should be the absolute difference between their respective
        # radii

        expected = 4.5
        actual = PSA.discrete_frechet(self.path_1,
                                      self.path_2)
        assert_almost_equal(actual, expected)


