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
from __future__ import print_function, absolute_import

import MDAnalysis as mda
import MDAnalysis.analysis.encore as encore

import importlib
import tempfile
import numpy as np
import sys
import os
import warnings

import pytest
from numpy.testing import assert_equal, assert_almost_equal

from MDAnalysisTests.datafiles import DCD, DCD2, PSF, TPR, XTC
from MDAnalysisTests import block_import

import MDAnalysis.analysis.rms as rms
import MDAnalysis.analysis.align as align
import MDAnalysis.analysis.encore.confdistmatrix as confdistmatrix


def function(x):
    return x**2

class TestEncore(object):
    @pytest.fixture(scope='class')
    def ens1_template(self):
        template = mda.Universe(PSF, DCD)
        template.transfer_to_memory(step=5)
        return template

    @pytest.fixture(scope='class')
    def ens2_template(self):
        template = mda.Universe(PSF, DCD2)
        template.transfer_to_memory(step=5)
        return template

    @pytest.fixture()
    def ens1(self, ens1_template):
        return mda.Universe(
            ens1_template.filename,
            ens1_template.trajectory.timeseries(order='fac'),
            format=mda.coordinates.memory.MemoryReader)

    @pytest.fixture()
    def ens2(self, ens2_template):
        return mda.Universe(
            ens2_template.filename,
            ens2_template.trajectory.timeseries(order='fac'),
            format=mda.coordinates.memory.MemoryReader)

    def test_triangular_matrix(self):
        scalar = 2
        size = 3
        expected_value = 1.984
        filename = tempfile.mktemp()+".npz"

        triangular_matrix = encore.utils.TriangularMatrix(size = size)

        triangular_matrix[0,1] = expected_value

        assert_equal(triangular_matrix[0,1], expected_value,
                     err_msg="Data error in TriangularMatrix: read/write are not consistent")

        assert_equal(triangular_matrix[0,1], triangular_matrix[1,0],
                        err_msg="Data error in TriangularMatrix: matrix non symmetrical")
        triangular_matrix.savez(filename)

        triangular_matrix_2 = encore.utils.TriangularMatrix(size = size, loadfile = filename)
        assert_equal(triangular_matrix_2[0,1], expected_value,
                        err_msg="Data error in TriangularMatrix: loaded matrix non symmetrical")

        triangular_matrix_3 = encore.utils.TriangularMatrix(size = size)
        triangular_matrix_3.loadz(filename)
        assert_equal(triangular_matrix_3[0,1], expected_value,
                        err_msg="Data error in TriangularMatrix: loaded matrix non symmetrical")

        incremented_triangular_matrix = triangular_matrix + scalar
        assert_equal(incremented_triangular_matrix[0,1], expected_value + scalar,
                     err_msg="Error in TriangularMatrix: addition of scalar gave"
                     "inconsistent results")

        triangular_matrix += scalar
        assert_equal(triangular_matrix[0,1], expected_value + scalar,
                     err_msg="Error in TriangularMatrix: addition of scalar gave"
                     "inconsistent results")

        multiplied_triangular_matrix_2 = triangular_matrix_2 * scalar
        assert_equal(multiplied_triangular_matrix_2[0,1], expected_value * scalar,
                     err_msg="Error in TriangularMatrix: multiplication by scalar gave"
                     "inconsistent results")

        triangular_matrix_2 *= scalar
        assert_equal(triangular_matrix_2[0,1], expected_value * scalar,
                        err_msg="Error in TriangularMatrix: multiplication by scalar gave\
inconsistent results")

    @pytest.mark.xfail(os.name == 'nt',
                       reason="Not yet supported on Windows.")
    def test_parallel_calculation(self):

        arguments = [tuple([i]) for i in np.arange(0,100)]

        parallel_calculation = encore.utils.ParallelCalculation(function=function,
                                                                n_jobs=4,
                                                                args=arguments)
        results = parallel_calculation.run()

        for i,r in enumerate(results):
            assert_equal(r[1], arguments[i][0]**2,
                err_msg="Unexpeted results from ParallelCalculation")

    def test_rmsd_matrix_with_superimposition(self, ens1):
        conf_dist_matrix = encore.confdistmatrix.conformational_distance_matrix(
            ens1,
            encore.confdistmatrix.set_rmsd_matrix_elements,
            select="name CA",
            pairwise_align=True,
            weights='mass',
            n_jobs=1)

        reference = rms.RMSD(ens1, select="name CA")
        reference.run()

        for i,rmsd in enumerate(reference.rmsd):
            assert_almost_equal(conf_dist_matrix[0,i], rmsd[2], decimal=3,
                                err_msg = "calculated RMSD values differ from the reference implementation")

    def test_rmsd_matrix_with_superimposition_custom_weights(self, ens1):
        conf_dist_matrix = encore.confdistmatrix.conformational_distance_matrix(
            ens1,
            encore.confdistmatrix.set_rmsd_matrix_elements,
            select="name CA",
            pairwise_align=True,
            weights='mass',
            n_jobs=1)

        conf_dist_matrix_custom = encore.confdistmatrix.conformational_distance_matrix(
            ens1,
            encore.confdistmatrix.set_rmsd_matrix_elements,
            select="name CA",
            pairwise_align=True,
            weights=(ens1.select_atoms('name CA').masses, ens1.select_atoms('name CA').masses),
            n_jobs=1)

        for i in range(conf_dist_matrix_custom.size):
            assert_almost_equal(conf_dist_matrix_custom[0, i], conf_dist_matrix[0, i])

    def test_rmsd_matrix_without_superimposition(self, ens1):
        selection_string = "name CA"
        selection = ens1.select_atoms(selection_string)
        reference_rmsd = []
        coordinates = ens1.trajectory.timeseries(selection, order='fac')
        for coord in coordinates:
            reference_rmsd.append(rms.rmsd(coordinates[0], coord, superposition=False))

        confdist_matrix = encore.confdistmatrix.conformational_distance_matrix(
            ens1,
            encore.confdistmatrix.set_rmsd_matrix_elements,
            select=selection_string,
            pairwise_align=False,
            weights='mass',
            n_jobs=1)

        print (repr(confdist_matrix.as_array()[0,:]))
        assert_almost_equal(confdist_matrix.as_array()[0,:], reference_rmsd, decimal=3,
                            err_msg="calculated RMSD values differ from reference")

    def test_ensemble_superimposition(self):
        aligned_ensemble1 = mda.Universe(PSF, DCD)
        align.AlignTraj(aligned_ensemble1, aligned_ensemble1,
                        select="name CA",
                        in_memory=True).run()
        aligned_ensemble2 = mda.Universe(PSF, DCD)
        align.AlignTraj(aligned_ensemble2, aligned_ensemble2,
                        select="name *",
                        in_memory=True).run()

        rmsfs1 = rms.RMSF(aligned_ensemble1.select_atoms('name *'))
        rmsfs1.run()

        rmsfs2 = rms.RMSF(aligned_ensemble2.select_atoms('name *'))
        rmsfs2.run()

        assert sum(rmsfs1.rmsf) > sum(rmsfs2.rmsf),"Ensemble aligned on all " \
                                                   "atoms should have lower full-atom RMSF than ensemble aligned on only CAs."

    def test_ensemble_superimposition_to_reference_non_weighted(self):
        aligned_ensemble1 = mda.Universe(PSF, DCD)
        align.AlignTraj(aligned_ensemble1, aligned_ensemble1,
                        select="name CA",
                        in_memory=True).run()
        aligned_ensemble2 = mda.Universe(PSF, DCD)
        align.AlignTraj(aligned_ensemble2, aligned_ensemble2,
                        select="name *",
                        in_memory=True).run()

        rmsfs1 = rms.RMSF(aligned_ensemble1.select_atoms('name *'))
        rmsfs1.run()

        rmsfs2 = rms.RMSF(aligned_ensemble2.select_atoms('name *'))
        rmsfs2.run()

        assert sum(rmsfs1.rmsf) > sum(rmsfs2.rmsf), "Ensemble aligned on all " \
                                                    "atoms should have lower full-atom RMSF than ensemble aligned on only CAs."

    def test_hes_to_self(self, ens1):
        results, details = encore.hes([ens1, ens1])
        result_value = results[0, 1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value,
                            err_msg="Harmonic Ensemble Similarity to itself not zero: {0:f}".format(result_value))

    def test_hes(self, ens1, ens2):
        results, details = encore.hes([ens1, ens2], weights='mass')
        result_value = results[0, 1]
        min_bound = 1E5
        assert result_value > min_bound, "Unexpected value for Harmonic " \
                                          "Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, min_bound)

    def test_hes_custom_weights(self, ens1, ens2):
        results, details = encore.hes([ens1, ens2], weights='mass')
        results_custom, details_custom = encore.hes([ens1, ens2],
                                                    weights=(ens1.select_atoms('name CA').masses,
                                                             ens2.select_atoms('name CA').masses))
        result_value = results[0, 1]
        result_value_custom = results_custom[0, 1]
        assert_almost_equal(result_value, result_value_custom)

    def test_hes_align(self, ens1, ens2):
        # This test is massively sensitive!
        # Get 5260 when masses were float32?
        results, details = encore.hes([ens1, ens2], align=True)
        result_value = results[0,1]
        expected_value = 2047.05
        assert_almost_equal(result_value, expected_value, decimal=-3,
                            err_msg="Unexpected value for Harmonic Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))

    def test_ces_to_self(self, ens1):
        results, details = \
            encore.ces([ens1, ens1],
            clustering_method=encore.AffinityPropagationNative(preference = -3.0))
        result_value = results[0,1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value,
                            err_msg="ClusteringEnsemble Similarity to itself not zero: {0:f}".format(result_value))

    def test_ces(self, ens1, ens2):
        results, details = encore.ces([ens1, ens2])
        result_value = results[0,1]
        expected_value = 0.51
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Unexpected value for Cluster Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))

    def test_dres_to_self(self, ens1):
        results, details = encore.dres([ens1, ens1])
        result_value = results[0,1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Dim. Reduction Ensemble Similarity to itself not zero: {0:f}".format(result_value))

    def test_dres(self, ens1, ens2):
        results, details = encore.dres([ens1, ens2], select="name CA and resnum 1-10")
        result_value = results[0,1]
        upper_bound = 0.6
        assert result_value < upper_bound, "Unexpected value for Dim. " \
                                            "reduction Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, upper_bound)

    @pytest.mark.xfail  # sporadically fails, see Issue #2158
    def test_dres_without_superimposition(self, ens1, ens2):
        distance_matrix = encore.get_distance_matrix(
            encore.merge_universes([ens1, ens2]),
            superimpose=False)
        results, details = encore.dres([ens1, ens2],
                                       distance_matrix = distance_matrix)
        result_value = results[0,1]
        expected_value = 0.68
        assert_almost_equal(result_value, expected_value, decimal=1,
                            err_msg="Unexpected value for Dim. reduction Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))

    def test_ces_convergence(self, ens1):
        expected_values = [0.3443593, 0.1941854, 0.06857104,  0.]
        results = encore.ces_convergence(ens1, 5)
        for i,ev in enumerate(expected_values):
            assert_almost_equal(ev, results[i], decimal=2,
                                err_msg="Unexpected value for Clustering Ensemble similarity in convergence estimation")

    def test_dres_convergence(self, ens1):
        # Due to encore.dres_convergence() involving random numbers, the
        # following assertion is allowed to fail once. This significantly
        # reduces the probability of a random test failure.
        expected_values = [0.3, 0.]
        results = encore.dres_convergence(ens1, 10)
        try:
            assert_almost_equal(results[:,0], expected_values, decimal=1)
        except AssertionError:
            # Random test failure is very rare, but repeating the failed test
            # just once would only assert that the test passes with 50%
            # probability. To be a little safer, we raise a warning and repeat
            # the test 10 times:
            warnings.warn(message="Test 'test_dres_convergence' failed, "
                                  "repeating test 10 times.",
                          category=RuntimeWarning)
            for i in range(10):
                results = encore.dres_convergence(ens1, 10)
                assert_almost_equal(results[:,0], expected_values, decimal=1,
                                    err_msg="Unexpected value for Dim. "
                                            "reduction Ensemble similarity in "
                                            "convergence estimation")

    @pytest.mark.xfail  # sporadically fails, see Issue #2158
    def test_hes_error_estimation(self, ens1):
        expected_average = 10
        expected_stdev = 12
        averages, stdevs = encore.hes([ens1, ens1], estimate_error = True, bootstrapping_samples=10, select="name CA and resnum 1-10")
        average = averages[0,1]
        stdev = stdevs[0,1]

        assert_almost_equal(average, expected_average, decimal=-2,
                            err_msg="Unexpected average value for bootstrapped samples in Harmonic Ensemble imilarity")
        assert_almost_equal(stdev, expected_stdev, decimal=-2,
                            err_msg="Unexpected standard daviation  for bootstrapped samples in Harmonic Ensemble imilarity")

    def test_ces_error_estimation(self, ens1):
        expected_average = 0.03
        expected_stdev = 0.31
        averages, stdevs = encore.ces([ens1, ens1],
                                      estimate_error = True,
                                      bootstrapping_samples=10,
                                      clustering_method=encore.AffinityPropagationNative(preference=-2.0),
                                      select="name CA and resnum 1-10")
        average = averages[0,1]
        stdev = stdevs[0,1]

        assert_almost_equal(average, expected_average, decimal=1,
                            err_msg="Unexpected average value for bootstrapped samples in Clustering Ensemble similarity")
        assert_almost_equal(stdev, expected_stdev, decimal=0,
                            err_msg="Unexpected standard daviation  for bootstrapped samples in Clustering Ensemble similarity")

    def test_ces_error_estimation_ensemble_bootstrap(self, ens1):
        # Error estimation using a method that does not take a distance
        # matrix as input, and therefore relies on bootstrapping the ensembles
        # instead

        pytest.importorskip('sklearn')

        expected_average = 0.03
        expected_stdev = 0.02
        averages, stdevs = encore.ces([ens1, ens1],
                                      estimate_error = True,
                                      bootstrapping_samples=10,
                                      clustering_method=encore.KMeans(n_clusters=2),
                                      select="name CA and resnum 1-10")
        average = averages[0,1]
        stdev = stdevs[0,1]

        assert_almost_equal(average, expected_average, decimal=1,
                            err_msg="Unexpected average value for bootstrapped samples in Clustering Ensemble similarity")
        assert_almost_equal(stdev, expected_stdev, decimal=1,
                            err_msg="Unexpected standard daviation  for bootstrapped samples in Clustering Ensemble similarity")

    def test_dres_error_estimation(self, ens1):
        average_upper_bound = 0.3
        stdev_upper_bound = 0.2
        averages, stdevs = encore.dres([ens1, ens1], estimate_error = True,
                                       bootstrapping_samples=10,
                                       select="name CA and resnum 1-10")
        average = averages[0,1]
        stdev = stdevs[0,1]

        assert average < average_upper_bound, "Unexpected average value for " \
                                               "bootstrapped samples in Dim. reduction Ensemble similarity"
        assert stdev < stdev_upper_bound, "Unexpected standard deviation for" \
                                           " bootstrapped samples in Dim. reduction Ensemble imilarity"

class TestEncoreClustering(object):
    @pytest.fixture(scope='class')
    def ens1_template(self):
        template = mda.Universe(PSF, DCD)
        template.transfer_to_memory(step=5)
        return template

    @pytest.fixture(scope='class')
    def ens2_template(self):
        template = mda.Universe(PSF, DCD2)
        template.transfer_to_memory(step=5)
        return template

    @pytest.fixture(scope='class')
    def cc(self):
        return encore.ClusterCollection([1, 1, 1, 3, 3, 5, 5, 5])

    @pytest.fixture(scope='class')
    def cluster(self):
        return encore.Cluster(elem_list=np.array([0, 1, 2]), centroid=1)

    @pytest.fixture()
    def ens1(self, ens1_template):
        return mda.Universe(
                ens1_template.filename,
                ens1_template.trajectory.timeseries(order='fac'),
                format=mda.coordinates.memory.MemoryReader)

    @pytest.fixture()
    def ens2(self, ens2_template):
        return mda.Universe(
            ens2_template.filename,
            ens2_template.trajectory.timeseries(order='fac'),
            format=mda.coordinates.memory.MemoryReader)
    
    def test_clustering_one_ensemble(self, ens1):
        cluster_collection = encore.cluster(ens1)
        expected_value = 7
        assert len(cluster_collection) == expected_value, "Unexpected " \
                                                          "results: {0}".format(cluster_collection)

    def test_clustering_two_ensembles(self, ens1, ens2):
        cluster_collection = encore.cluster([ens1, ens2])
        expected_value = 14
        assert len(cluster_collection) == expected_value, "Unexpected " \
                                                          "results: {0}".format(cluster_collection)

    def test_clustering_three_ensembles_two_identical(self, ens1, ens2):
        cluster_collection = encore.cluster([ens1, ens2, ens1])
        expected_value = 40
        assert len(cluster_collection) == expected_value, "Unexpected result:" \
                                                          " {0}".format(cluster_collection)

    def test_clustering_two_methods(self, ens1):
        cluster_collection = encore.cluster(
            [ens1],
            method=[encore.AffinityPropagationNative(),
                    encore.AffinityPropagationNative()])
        assert len(cluster_collection[0]) == len(cluster_collection[1]), \
                     "Unexpected result: {0}".format(cluster_collection)

    def test_clustering_AffinityPropagationNative_direct(self, ens1):
        method = encore.AffinityPropagationNative()
        distance_matrix = encore.get_distance_matrix(ens1)
        cluster_assignment = method(distance_matrix)
        expected_value = 7
        assert len(set(cluster_assignment)) == expected_value, \
                     "Unexpected result: {0}".format(cluster_assignment)

    def test_clustering_AffinityPropagation_direct(self, ens1):
        pytest.importorskip('sklearn')
        method = encore.AffinityPropagation()
        distance_matrix = encore.get_distance_matrix(ens1)
        cluster_assignment = method(distance_matrix)
        expected_value = 7
        assert len(set(cluster_assignment)) == expected_value, \
                     "Unexpected result: {0}".format(cluster_assignment)

    def test_clustering_KMeans_direct(self, ens1):
        pytest.importorskip('sklearn')
        clusters = 10
        method = encore.KMeans(clusters)
        coordinates = ens1.trajectory.timeseries(order='fac')
        coordinates = np.reshape(coordinates,
                                 (coordinates.shape[0], -1))
        cluster_assignment = method(coordinates)
        assert len(set(cluster_assignment)) == clusters, \
                     "Unexpected result: {0}".format(cluster_assignment)

    def test_clustering_DBSCAN_direct(self, ens1):
        pytest.importorskip('sklearn')
        method = encore.DBSCAN(eps=0.5, min_samples=2)
        distance_matrix = encore.get_distance_matrix(ens1)
        cluster_assignment = method(distance_matrix)
        expected_value = 2
        assert len(set(cluster_assignment)) == expected_value, \
                     "Unexpected result: {0}".format(cluster_assignment)

    def test_clustering_two_different_methods(self, ens1):
        pytest.importorskip('sklearn')
        cluster_collection = encore.cluster(
            [ens1],
            method=[encore.AffinityPropagation(preference=-7.5),
                    encore.DBSCAN(min_samples=2)])
        assert len(cluster_collection[0]) == len(cluster_collection[1]), \
                     "Unexpected result: {0}".format(cluster_collection)

    def test_clustering_method_w_no_distance_matrix(self, ens1):
        pytest.importorskip('sklearn')
        cluster_collection = encore.cluster(
            [ens1],
            method=encore.KMeans(10))
        assert len(cluster_collection) == 10, \
                     "Unexpected result: {0}".format(cluster_collection)

    def test_clustering_two_methods_one_w_no_distance_matrix(self, ens1):
        pytest.importorskip('sklearn')
        cluster_collection = encore.cluster(
            [ens1],
            method=[encore.KMeans(17),
                    encore.AffinityPropagationNative()])
        assert len(cluster_collection[0]) == len(cluster_collection[0]), \
                     "Unexpected result: {0}".format(cluster_collection)

    def test_sklearn_affinity_propagation(self, ens1):
        pytest.importorskip('sklearn')
        cc1 = encore.cluster([ens1])
        cc2 = encore.cluster([ens1],
                             method=encore.AffinityPropagation())
        assert len(cc1) == len(cc2), \
                     "Native and sklearn implementations of affinity "\
                              "propagation don't agree: mismatch in number of "\
                              "clusters: {0} {1}".format(len(cc1), len(cc2))

    def test_ClusterCollection_init(self, cc):
        assert np.all(cc.clusters[0].elements == [0, 1, 2]) and \
               np.all(cc.clusters[1].elements == [3, 4   ]) and \
               np.all(cc.clusters[2].elements == [5, 6, 7]) and \
               cc.clusters[0].centroid == 1 and \
               cc.clusters[1].centroid == 3 and \
               cc.clusters[2].centroid == 5, \
                      "ClusterCollection was not constructed correctly"

    def test_Cluster_init(self, cluster):
        assert np.all(cluster.elements == [0, 1, 2]) and \
               cluster.centroid == 1, \
                      "Cluster was not constructed correctly"

    def test_ClusterCollection_get_ids(self, cc):
        assert cc.get_ids() == [0, 1, 2], \
                     "ClusterCollection ids aren't as expected"

    def test_ClusterCollection_get_centroids(self, cc):
        assert cc.get_centroids() == [1, 3, 5], \
                     "ClusterCollection centroids aren't as expected"

    def test_Cluster_add_metadata(self, cluster):
        metadata = cluster.elements*10
        cluster.add_metadata('test', metadata)
        assert np.all(cluster.metadata['test'] == metadata), \
                     "Cluster metadata isn't as expected"

class TestEncoreClusteringSklearn(object):
    """The tests in this class were duplicated from the affinity propagation
    tests in scikit-learn"""

    n_clusters = 3

    @pytest.fixture()
    def distance_matrix(self):
        X = np.array([[8.73101582, 8.85617874],
                      [11.61311169, 11.58774351],
                      [10.86083514, 11.06253959],
                      [9.45576027, 8.50606967],
                      [11.30441509, 11.04867001],
                      [8.63708065, 9.02077816],
                      [8.34792066, 9.1851129],
                      [11.06197897, 11.15126501],
                      [11.24563175, 9.36888267],
                      [10.83455241, 8.70101808],
                      [11.49211627, 11.48095194],
                      [10.6448857, 10.20768141],
                      [10.491806, 9.38775868],
                      [11.08330999, 9.39065561],
                      [10.83872922, 9.48897803],
                      [11.37890079, 8.93799596],
                      [11.70562094, 11.16006288],
                      [10.95871246, 11.1642394],
                      [11.59763163, 10.91793669],
                      [11.05761743, 11.5817094],
                      [8.35444086, 8.91490389],
                      [8.79613913, 8.82477028],
                      [11.00420001, 9.7143482],
                      [11.90790185, 10.41825373],
                      [11.39149519, 11.89635728],
                      [8.31749192, 9.78031016],
                      [11.59530088, 9.75835567],
                      [11.17754529, 11.13346973],
                      [11.01830341, 10.92512646],
                      [11.75326028, 8.46089638],
                      [11.74702358, 9.36241786],
                      [10.53075064, 9.77744847],
                      [8.67474149, 8.30948696],
                      [11.05076484, 9.16079575],
                      [8.79567794, 8.52774713],
                      [11.18626498, 8.38550253],
                      [10.57169895, 9.42178069],
                      [8.65168114, 8.76846013],
                      [11.12522708, 10.6583617],
                      [8.87537899, 9.02246614],
                      [9.29163622, 9.05159316],
                      [11.38003537, 10.93945712],
                      [8.74627116, 8.85490353],
                      [10.65550973, 9.76402598],
                      [8.49888186, 9.31099614],
                      [8.64181338, 9.154761],
                      [10.84506927, 10.8790789],
                      [8.98872711, 9.17133275],
                      [11.7470232, 10.60908885],
                      [10.89279865, 9.32098256],
                      [11.14254656, 9.28262927],
                      [9.02660689, 9.12098876],
                      [9.16093666, 8.72607596],
                      [11.47151183, 8.92803007],
                      [11.76917681, 9.59220592],
                      [9.97880407, 11.26144744],
                      [8.58057881, 8.43199283],
                      [10.53394006, 9.36033059],
                      [11.34577448, 10.70313399],
                      [9.07097046, 8.83928763]])

        XX = np.einsum('ij,ij->i', X, X)[:, np.newaxis]
        YY = XX.T
        distances = np.dot(X, X.T)
        distances *= -2
        distances += XX
        distances += YY
        np.maximum(distances, 0, out=distances)
        distances.flat[::distances.shape[0] + 1] = 0.0
        dimension = len(distances)

        distance_matrix = encore.utils.TriangularMatrix(len(distances))
        for i in range(dimension):
            for j in range(i, dimension):
                distance_matrix[i, j] = distances[i, j]
        return distance_matrix

    def test_one(self, distance_matrix):
        preference = -float(np.median(distance_matrix.as_array()) * 10.)
        clustering_method = encore.AffinityPropagationNative(preference=preference)
        ccs = encore.cluster(None,
                             distance_matrix=distance_matrix,
                             method=clustering_method)
        assert self.n_clusters == len(ccs), \
                     "Basic clustering test failed to give the right"\
                    "number of clusters: {0} vs {1}".format(self.n_clusters, len(ccs))


class TestEncoreDimensionalityReduction(object):
    @pytest.fixture(scope='class')
    def ens1_template(self):
        template = mda.Universe(PSF, DCD)
        template.transfer_to_memory(step=5)
        return template

    @pytest.fixture(scope='class')
    def ens2_template(self):
        template = mda.Universe(PSF, DCD2)
        template.transfer_to_memory(step=5)
        return template

    @pytest.fixture()
    def ens1(self, ens1_template):
        return mda.Universe(
            ens1_template.filename,
            ens1_template.trajectory.timeseries(order='fac'),
            format=mda.coordinates.memory.MemoryReader)

    @pytest.fixture()
    def ens2(self, ens2_template):
        return mda.Universe(
            ens2_template.filename,
            ens2_template.trajectory.timeseries(order='fac'),
            format=mda.coordinates.memory.MemoryReader)

    def test_dimensionality_reduction_one_ensemble(self, ens1):
        dimension = 2
        coordinates, details = encore.reduce_dimensionality(ens1)
        assert_equal(coordinates.shape[0], dimension,
                     err_msg="Unexpected result in dimensionality reduction: {0}".format(coordinates))


    def test_dimensionality_reduction_two_ensembles(self, ens1, ens2):
        dimension = 2
        coordinates, details = \
            encore.reduce_dimensionality([ens1, ens2])
        assert_equal(coordinates.shape[0], dimension,
                     err_msg="Unexpected result in dimensionality reduction: {0}".format(coordinates))


    def test_dimensionality_reduction_three_ensembles_two_identical(self,
                                                                    ens1, ens2):
        coordinates, details = \
            encore.reduce_dimensionality([ens1, ens2, ens1])
        coordinates_ens1 = coordinates[:,np.where(details["ensemble_membership"]==1)]
        coordinates_ens3 = coordinates[:,np.where(details["ensemble_membership"]==3)]
        assert_almost_equal(coordinates_ens1, coordinates_ens3, decimal=0,
                     err_msg="Unexpected result in dimensionality reduction: {0}".format(coordinates))


    def test_dimensionality_reduction_specified_dimension(self, ens1, ens2):
        dimension = 3
        coordinates, details = encore.reduce_dimensionality(
            [ens1, ens2],
            method=encore.StochasticProximityEmbeddingNative(dimension=dimension))
        assert_equal(coordinates.shape[0], dimension,
                     err_msg="Unexpected result in dimensionality reduction: {0}".format(coordinates))


    def test_dimensionality_reduction_SPENative_direct(self, ens1):
        dimension = 2
        method = encore.StochasticProximityEmbeddingNative(dimension=dimension)
        distance_matrix = encore.get_distance_matrix(ens1)
        coordinates, details = method(distance_matrix)
        assert_equal(coordinates.shape[0], dimension,
                     err_msg="Unexpected result in dimensionality reduction: {0}".format(
                     coordinates))

    def test_dimensionality_reduction_PCA_direct(self, ens1):
        pytest.importorskip('sklearn')
        dimension = 2
        method = encore.PrincipalComponentAnalysis(dimension=dimension)
        coordinates = ens1.trajectory.timeseries(order='fac')
        coordinates = np.reshape(coordinates,
                                 (coordinates.shape[0], -1))
        coordinates, details = method(coordinates)
        assert_equal(coordinates.shape[0], dimension,
                     err_msg="Unexpected result in dimensionality reduction: {0}".format(
                     coordinates))


    def test_dimensionality_reduction_different_method(self, ens1, ens2):
        pytest.importorskip('sklearn')
        dimension = 3
        coordinates, details = \
            encore.reduce_dimensionality(
                [ens1, ens2],
                method=encore.PrincipalComponentAnalysis(dimension=dimension))
        assert_equal(coordinates.shape[0], dimension,
                     err_msg="Unexpected result in dimensionality reduction: {0}".format(coordinates))


    def test_dimensionality_reduction_two_methods(self, ens1, ens2):
        dims = [2,3]
        coordinates, details = \
            encore.reduce_dimensionality(
                [ens1, ens2],
                method=[encore.StochasticProximityEmbeddingNative(dims[0]),
                        encore.StochasticProximityEmbeddingNative(dims[1])])
        assert_equal(coordinates[1].shape[0], dims[1])

    def test_dimensionality_reduction_two_different_methods(self, ens1, ens2):
        pytest.importorskip('sklearn')
        dims = [2,3]
        coordinates, details = \
            encore.reduce_dimensionality(
                [ens1, ens2],
                method=[encore.StochasticProximityEmbeddingNative(dims[0]),
                        encore.PrincipalComponentAnalysis(dims[1])])
        assert_equal(coordinates[1].shape[0], dims[1])


class TestEncoreConfDistMatrix(object):
    def test_get_distance_matrix(self):
        # Issue #1324
        u = mda.Universe(TPR,XTC)
        dm = confdistmatrix.get_distance_matrix(u)

class TestEncoreImportWarnings(object):
    @block_import('sklearn')
    def _check_sklearn_import_warns(self, package, recwarn):
        for mod in list(sys.modules):  # list as we're changing as we iterate
            if 'encore' in mod:
                sys.modules.pop(mod, None)
        warnings.simplefilter('always')
        # assert_warns(ImportWarning, importlib.import_module, package)
        importlib.import_module(package)
        assert recwarn.pop(ImportWarning)

    def test_import_warnings(self, recwarn):
        for mod in list(sys.modules):  # list as we're changing as we iterate
            if 'encore' in mod:
                sys.modules.pop(mod, None)
        for pkg in (
                'MDAnalysis.analysis.encore.dimensionality_reduction.DimensionalityReductionMethod',
                'MDAnalysis.analysis.encore.clustering.ClusteringMethod',
        ):
            self._check_sklearn_import_warns(pkg, recwarn)
            # This is a quickfix! Convert this to a parametrize call in future.
