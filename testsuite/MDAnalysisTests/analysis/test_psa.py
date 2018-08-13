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
from __future__ import division, absolute_import, print_function

import pytest
from six.moves import range

import MDAnalysis as mda
import MDAnalysis.analysis.psa as PSA

from numpy.testing import (assert_array_less,
                           assert_almost_equal, assert_equal)
import numpy as np
import scipy
import scipy.spatial
import matplotlib

from MDAnalysisTests.datafiles import PSF, DCD, DCD2


class TestPSAnalysis(object):
    iu1 = np.triu_indices(3, k=1)

    @pytest.fixture()
    def psa(self, tmpdir):
        universe1 = mda.Universe(PSF, DCD)
        universe2 = mda.Universe(PSF, DCD2)
        universe_rev = mda.Universe(PSF, DCD)

        psa = PSA.PSAnalysis([universe1, universe2, universe_rev],
                             path_select='name CA',
                             targetdir=str(tmpdir))

        psa.generate_paths(align=True)
        psa.paths[-1] = psa.paths[-1][::-1, :, :]  # reverse third path
        return psa

    @pytest.fixture()
    def hausd_matrix(self, psa):
        psa.run(metric='hausdorff')
        return psa.get_pairwise_distances()

    @pytest.fixture()
    def hausd_dists(self, hausd_matrix):
        return hausd_matrix[self.iu1]

    @pytest.fixture()
    def frech_matrix(self, psa):
        psa.run(metric='discrete_frechet')
        return psa.get_pairwise_distances()

    @pytest.fixture()
    def frech_dists(self, frech_matrix):
        return frech_matrix[self.iu1]

    @pytest.fixture()
    def plot_data(self, psa, tmpdir):
        psa.run(metric='hausdorff')
        psa.run(metric='discrete_frechet')
        with tmpdir.as_cwd():
            results = psa.plot(filename="distmat.png")
        return results

    @pytest.fixture()
    def plot_annotated_heatmap(self, psa, tmpdir):
        pytest.importorskip('seaborn')
        psa.run(metric='hausdorff')
        with tmpdir.as_cwd():
            results = psa.plot_annotated_heatmap(filename="annotated.png")
        return results

    @pytest.fixture()
    def plot_nearest_neighbors(self, psa, tmpdir):
        pytest.importorskip('seaborn')
        psa.run(metric='hausdorff')
        psa.run_pairs_analysis(neighbors=True)
        with tmpdir.as_cwd():
            results = psa.plot_nearest_neighbors(filename="nn.png")
        return results

    @pytest.fixture()
    def hausd_pairs_dists(self, psa):
        psa.run_pairs_analysis(neighbors=True, hausdorff_pairs=True)
        hausd_pairs_matrix = psa.get_pairwise_distances()

        return hausd_pairs_matrix[self.iu1]

    def test_hausdorff_bound(self, hausd_dists, frech_dists):
        """Test whether Frechet distances are smaller than corresponding
        Hausdorff distances"""
        err_msg = ("Some Frechet distances are smaller than corresponding "
                   "Hausdorff distances")
        assert_array_less(hausd_dists, frech_dists, err_msg)

    def test_explicit_metric(self, psa, hausd_dists):
        """Test whether explicitly specifying Hausdorff metric gives same result
        as specifying Hausdorff metric with string name"""
        psa.run(metric=PSA.hausdorff)
        hausd_matrix_explicit = psa.get_pairwise_distances()
        hausd_explicit_dists = hausd_matrix_explicit[self.iu1]

        err_msg = ("Specifying Python function for Hausdoff gives different "
                   "distances than specifying Hausdorff with string name")
        assert_equal(hausd_dists, hausd_explicit_dists, err_msg)

    def test_reversal_hausdorff(self, hausd_matrix):
        """Test whether Hausdorff distances are invariant to path reversal"""
        err_msg = "Hausdorff distances changed after path reversal"
        assert_almost_equal(hausd_matrix[1,2],
                            hausd_matrix[0,1],
                            decimal=3, err_msg=err_msg)

    def test_reversal_frechet(self, frech_matrix):
        """Test whether Frechet distances are same/larger after path reversal"""
        err_msg = "Frechet distances did not increase after path reversal"
        assert frech_matrix[1,2] >= frech_matrix[0,1], err_msg

    def test_get_num_paths(self, psa):
        assert psa.get_num_paths() == 3

    def test_get_paths(self, psa):
        paths = psa.get_paths()
        assert len(paths) == 3
        assert isinstance(paths, list)

    def test_psa_pairs_ValueError(self, psa):
        with pytest.raises(ValueError):
            psa.psa_pairs

    def test_psa_pairs(self, psa):
        psa.run_pairs_analysis()
        assert len(psa.psa_pairs) == 3

    def test_hausdorff_pairs_ValueError(self, psa):
        with pytest.raises(ValueError):
            psa.hausdorff_pairs

    def test_hausdorff_pairs(self, psa):
        psa.run_pairs_analysis(hausdorff_pairs=True)
        assert len(psa.hausdorff_pairs) == 3

    def test_nearest_neighbors_ValueError(self, psa):
        with pytest.raises(ValueError):
            psa.nearest_neighbors

    def test_nearest_neighbors(self, psa):
        psa.run_pairs_analysis(neighbors=True)
        assert len(psa.nearest_neighbors) == 3

    @pytest.mark.xfail
    def test_load(self, psa):
        """Test that the automatically saved files can be loaded"""
        expected_path_names = psa.path_names[:]
        expected_paths = [p.copy() for p in psa.paths]
        psa.save_paths()
        psa.load()
        assert psa.path_names == expected_path_names
        # manually compare paths because
        #         assert_almost_equal(psa.paths, expected_paths, decimal=6)
        # raises a ValueError in the assertion code itself
        assert len(psa.paths) == len(expected_paths)
        for ipath, (observed, expected) in enumerate(zip(psa.paths, expected_paths)):
            assert_almost_equal(observed, expected, decimal=6,
                                err_msg="loaded path {} does not agree with input".format(ipath))

    def test_dendrogram_produced(self, plot_data):
        """Test whether Dendrogram dictionary object was produced"""
        err_msg = "Dendrogram dictionary object was not produced"
        assert isinstance(plot_data[1], dict), err_msg

    def test_dendrogram_produced_annotated(self, plot_annotated_heatmap):
        """Test whether Dendrogram dictionary object was produced"""
        err_msg = "Dendrogram dictionary object was not produced"
        assert isinstance(plot_annotated_heatmap[1], dict), err_msg

    def test_plot_nearest_neighbors(self, plot_nearest_neighbors):
        assert isinstance(plot_nearest_neighbors, matplotlib.axes.Axes)

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
        assert_equal(
            PSA.dist_mat_to_vec(np.int32(5), np.int32(3), np.int32(4)), np.int32(9),
            err_msg)

    def test_dist_mat_to_vec_input_numpy_integer_16(self):
        """Test whether inputs are supported as numpy integers rather than normal Integers"""
        err_msg = "dist_mat_to_vec function returning wrong values"
        assert_equal(
            PSA.dist_mat_to_vec(np.int16(5), np.int16(3), np.int16(4)), np.int16(9),
            err_msg)

    def test_hausdorff_pairs_distances(self, hausd_dists, hausd_pairs_dists):
        """Test whether Hausdorff pairs analysis distances are
        identical to those from standard Hausdorff metric"""
        err_msg = ("Some Hausdorff distances from pairs analysis vary "
                   "significantly from usual Hausdorff calculation")
        assert_almost_equal(hausd_dists, hausd_pairs_dists,
                            decimal=6, err_msg=err_msg)

    def test_distances_from_hausdorff_pairs_frames(self, psa):
        """Test whether actual distances between frames of Hausdorff
        pairs of a path give the expected Hausdorff distance"""
        psa.run_pairs_analysis(neighbors=True, hausdorff_pairs=True)
        hausd_pairs = psa.hausdorff_pairs
        npairs = int(psa.npaths * (psa.npaths - 1) / 2)
        hausd_pairs_dists2 = np.array([hausd_pairs[i]['distance']
                                       for i in range(npairs)])


        err_msg = ("A Hausdorff pair analysis distance when accessed "
                   "by frame varies from expected Hausdorff distance")
        dists = np.zeros((psa.npaths, psa.npaths))
        for i in range(0, psa.npaths-1):
            for j in range(i+1, psa.npaths):
                pairidx = PSA.dist_mat_to_vec(psa.npaths, i, j)
                p, q = hausd_pairs[pairidx]['frames']
                dists[i,j] = (PSA.sqnorm(psa.paths[i][p,:,:] -
                                         psa.paths[j][q,:,:]) /
                                         psa.natoms)**0.5
        assert_almost_equal(hausd_pairs_dists2,
                            dists[self.iu1],
                            decimal=6, err_msg=err_msg)

class TestPSAExceptions(object):
    '''Tests for exceptions that should be raised
    or caught by code in the psa module.'''

    def test_get_path_metric_func_bad_key(self):
        '''Test that KeyError is caught by
        get_path_metric_func().'''
        with pytest.raises(KeyError):
            PSA.get_path_metric_func('123456')

    def test_get_coord_axes_bad_dims(self):
        """Test that ValueError is raised when
        numpy array with incorrect dimensions
        is fed to get_coord_axes()."""
        with pytest.raises(ValueError):
            PSA.get_coord_axes(np.zeros((5,5,5,5)))

    @pytest.mark.parametrize('N, i, j', (
        (5, 6, 4),
        (5, 4, 6),
        (5, 6, 7),
        (5, -1, 2),
        (5, 1, -2),
        (1, 0, 0)

    ))
    def test_dist_mat_to_vec_func_out_of_bounds(self, N, i, j):
        """Test that ValueError is raised when i or j or both are
        out of bounds of N"""

        # Check if i is out of bounds of N
        with pytest.raises(ValueError):
            PSA.dist_mat_to_vec(N, i, j)

    @pytest.mark.parametrize('N, i, j', (
        (5, 4, 4),
        (4, 6, 4)
    ))
    def test_dist_mat_to_vec_func_i_equals_j(self, N, i, j):
        """Test that ValueError is raised when i == j or i,j == N"""

        with pytest.raises(ValueError):
            PSA.dist_mat_to_vec(N, i, j)

    def test_dist_mat_to_vec_func_bad_integers(self):
        """Test that ValueError is raised when i or j are
        not Integers"""

        with pytest.raises(ValueError) as err:
            PSA.dist_mat_to_vec(5, '6', '7')
        assert 'all must be of type int' in str(err.value)

        with pytest.raises(ValueError):
            PSA.dist_mat_to_vec(5, float(6), 7)


class _BaseHausdorffDistance(object):
    '''Base Class setup and unit tests
    for various Hausdorff distance
    calculation properties.'''

    @pytest.fixture()
    def random_angles(self):
        return np.random.random((100,)) * np.pi * 2

    @staticmethod
    @pytest.fixture()
    def path_1(random_angles):
        random_columns = np.column_stack((random_angles,
                                          random_angles,
                                          np.zeros((100,))))
        random_columns[..., 0] = np.cos(random_columns[..., 0])
        random_columns[..., 1] = np.sin(random_columns[..., 1])

        return random_columns

    @staticmethod
    @pytest.fixture()
    def path_2(random_angles):
        random_columns_2 = np.column_stack((random_angles,
                                            random_angles,
                                            np.zeros((100,))))
        random_columns_2[1:, 0] = np.cos(random_columns_2[1:, 0]) * 2.0
        random_columns_2[1:, 1] = np.sin(random_columns_2[1:, 1]) * 2.0
        # move one point farther out so we don't have two perfect circles
        random_columns_2[0, 0] = np.cos(random_columns_2[0, 0]) * 3.3
        random_columns_2[0, 1] = np.sin(random_columns_2[0, 1]) * 3.3
        return random_columns_2

    def test_symmetry(self, path_1, path_2, h):
        '''Ensure that the undirected (symmetric)
        Hausdorff distance is actually symmetric
        for a given Hausdorff metric, h.'''
        forward = h(path_1, path_2)
        reverse = h(path_2, path_1)
        # lower precision on 32bit
        assert_almost_equal(forward, reverse, decimal=15)

    def test_hausdorff_value(self, path_1, path_2, h, expected):
        '''Test that the undirected Hausdorff
        distance matches expected value for
        the simple case here.'''
        actual = h(path_1, path_2)
        # unless I pin down the random generator
        # seems unstable to use decimal > 2
        assert_almost_equal(actual, expected, decimal=2)


class TestHausdorffSymmetric(_BaseHausdorffDistance):
    '''Tests for conventional and symmetric (undirected)
    Hausdorff distance between point sets in 3D.'''

    # expected = 2.3

    @pytest.fixture()
    def expected(self):
        return 2.3

    @pytest.fixture()
    def h(self):
        return PSA.hausdorff


class TestWeightedAvgHausdorffSymmetric(_BaseHausdorffDistance):
    '''Tests for weighted average and symmetric (undirected)
    Hausdorff distance between point sets in 3D.'''

    @pytest.fixture()
    def expected(self, path_1, path_2):
        distance_matrix = scipy.spatial.distance.cdist(path_1, path_2)
        return (np.mean(np.amin(distance_matrix, axis=0)) +
                np.mean(np.amin(distance_matrix, axis=1))) / 2.0

    # params instead of 2 (I think it sends self as a parameter too.)
    @pytest.fixture()
    def h(self):
        return PSA.hausdorff_wavg

    def test_asymmetric_weight(self, path_1, path_2, h):
        '''Test for WAvg Hausdorff to ensure that increasing N points in one
        of the paths does NOT increase the weight of its contributions.'''
        inflated_path_1 = np.concatenate((path_1, path_1))
        inflated_path_2 = np.concatenate((path_2, path_2))
        d_inner_inflation = h(inflated_path_1, path_2)
        d_outer_inflation = h(path_1, inflated_path_2)
        assert_almost_equal(d_inner_inflation, d_outer_inflation)


class TestAvgHausdorffSymmetric(_BaseHausdorffDistance):
    '''Tests for unweighted average and symmetric (undirected)
    Hausdorff distance between point sets in 3D.'''

    @pytest.fixture()
    def expected(self, path_1, path_2):
        distance_matrix = scipy.spatial.distance.cdist(path_1, path_2)
        return np.mean(np.append(np.amin(distance_matrix, axis=0),
                                 np.amin(distance_matrix, axis=1)))

    # params instead of 2 (I think it sends self as a parameter too.)
    @pytest.fixture()
    def h(self):
        return PSA.hausdorff_avg

    def test_asymmetric_weight(self, path_1, path_2, h):
        '''Test to ensure that increasing N points in one of the paths
        increases the weight of its contributions.'''
        inflated_path_1 = np.concatenate((path_1, path_1))
        inflated_path_2 = np.concatenate((path_2, path_2))
        d_inner_inflation = h(inflated_path_1, path_2)
        d_outer_inflation = h(path_1, inflated_path_2)

        assert_array_less(d_inner_inflation, d_outer_inflation)


class DiscreteFrechetDistance(object):
    @staticmethod
    @pytest.fixture()
    def random_angles():
        return np.random.random((100,)) * np.pi * 2

    @staticmethod
    @pytest.fixture()
    def path_1(random_angles):
        random_columns = np.column_stack((random_angles, random_angles,
                                          np.zeros((100,))))
        random_columns[..., 0] = np.cos(random_columns[..., 0])
        random_columns[..., 1] = np.sin(random_columns[..., 1])
        return random_columns

    @staticmethod
    @pytest.fixture()
    def path_2(random_angles):
        random_columns_2 = np.column_stack((random_angles, random_angles,
                                            np.zeros((100,))))
        random_columns_2[..., 0] = np.cos(random_columns_2[..., 0]) * 5.5
        random_columns_2[..., 1] = np.sin(random_columns_2[..., 1]) * 5.5
        return random_columns_2

    def test_discrete_Frechet_concentric_circles(self, path_1, path_2):
        # test for the simple case of the discrete Frechet distance
        # between concentric circular paths, which for a sufficiently
        # high random discrete point density around each circle
        # should be the absolute difference between their respective
        # radii

        expected = 4.5
        actual = PSA.discrete_frechet(path_1, path_2)
        assert_almost_equal(actual, expected)
