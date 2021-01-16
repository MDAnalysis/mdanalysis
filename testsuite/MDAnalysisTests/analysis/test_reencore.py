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
import MDAnalysis as mda
from MDAnalysis.analysis import encore, reencore, align, pca
from MDAnalysis.analysis.diffusionmap import DistanceMatrix

from MDAnalysis.analysis.clustering import Clusters, methods
from MDAnalysis.analysis.encore.covariance import (ml_covariance_estimator,
                                                   shrinkage_covariance_estimator)
from MDAnalysis.analysis.encore.dimensionality_reduction import DimensionalityReductionMethod as drm

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

np.random.seed(0)

@pytest.fixture()
def data():
    return np.arange(24, dtype=np.float64).reshape((4, 6))


def test_max_likelihood_estimator(data):
    old = ml_covariance_estimator(data, 0)
    new = reencore.utils.max_likelihood_covariance(data, False)
    assert_almost_equal(new, old)

def test_accurate_default(data):
    old = shrinkage_covariance_estimator(data)
    new = reencore.utils.shrinkage_covariance(data)
    assert_almost_equal(new, old)


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

    @pytest.fixture()
    def ensemble(self, ens1, ens2):
        return mda.Ensemble([ens1, ens2]).select_atoms("name CA")
    
    @pytest.fixture()
    def ensemble1(self, ens1):
        return mda.Ensemble([ens1, ens1]).select_atoms("name CA")

    @pytest.fixture()
    def ensemble_aligned(self):
        u1 = mda.Universe(PSF, DCD)
        u1.transfer_to_memory(step=5)
        u2 = mda.Universe(PSF, DCD2)
        u2.transfer_to_memory(step=5)
        ens = mda.Ensemble([u1, u2]).select_atoms("name CA")
        a = align.AlignTraj(ens, ens, in_memory=True).run()
        return ens
    
    @pytest.fixture()
    def ensemble1_aligned(self):
        u1 = mda.Universe(PSF, DCD)
        u1.transfer_to_memory(step=5)
        u2 = mda.Universe(PSF, DCD)
        u2.transfer_to_memory(step=5)
        ens = mda.Ensemble([u1, u2]).select_atoms("name CA")
        a = align.AlignTraj(ens, ens, in_memory=True).run()
        return ens

    @pytest.fixture()
    def ens1_aligned(self):
        u = mda.Universe(PSF, DCD)
        u.transfer_to_memory(step=5)
        align.AlignTraj(u, u, select="name CA", in_memory=True).run()
        return u

    @pytest.fixture()
    def dist_mat(self, ensemble_aligned):
        return encore.get_distance_matrix(ensemble_aligned)

    @pytest.fixture()
    def dist_mat1(self, ensemble1_aligned):
        return encore.get_distance_matrix(ensemble1_aligned)

    def test_affinity_propagation(self, ens1):
        dist_mat = encore.get_distance_matrix(ens1).as_array()
        clusters = Clusters(methods.AffinityPropagation)
        clusters.run(-dist_mat)
        assert len(clusters.cluster_indices) == 7

    def test_hes_to_self(self, ensemble1):
        result_value = reencore.hes(ensemble1)[0, 1]
        assert_almost_equal(result_value, 0,
                            err_msg=f"hes() to itself not zero: {result_value}")
    
    def test_hes(self, ensemble):
        result = reencore.hes(ensemble, weights='mass')[0, 1]
        min_bound = 1E5
        assert result > min_bound, "Unexpected value for Harmonic Ensemble Similarity"

    def test_hes_align(self, ens1, ens2):
        ensemble = mda.Ensemble([ens1, ens2]).select_atoms("name CA")
        result_value = reencore.hes(ensemble, align=True)[0, 1]
        old, _ = encore.hes(ensemble.universes)
        assert_almost_equal(result_value, old[0, 1], decimal=-2)
    
    def test_hes_estimate_error(self, ensemble1):
        omean, ostd = encore.hes(ensemble1.universes, estimate_error=True,
                                 bootstrapping_samples=10,
                                 select="name CA and resnum 1-10")
        mean, std = reencore.hes(ensemble1, estimate_error=True,
                                 select="name CA and resnum 1-10",
                                 n_bootstrap_samples=10)
        assert_almost_equal(mean[0, 1], mean[0, 1], decimal=1)
        assert_almost_equal(std[0, 1], std[0, 1], decimal=1)


    def test_ces_to_self(self, ensemble1_aligned, dist_mat1):
        clusters = Clusters(methods.AffinityPropagation(preference=-3.0))
        clusters.run(-dist_mat1.as_array())
        result_value = reencore.ces(ensemble1_aligned, clusters)[0, 1]
        assert_almost_equal(result_value, 0,
                            err_msg=f"ces() to itself not zero: {result_value}")

    def test_ces_rmsd_enc(self, ensemble_aligned, dist_mat):
        rmsd_mat_enc = dist_mat.as_array()
        clusters = Clusters(methods.AffinityPropagation())
        clusters.run(-rmsd_mat_enc)
        result_value = reencore.ces(ensemble_aligned, clusters)[0, 1]
        assert_almost_equal(result_value, 0.51, decimal=2,
                            err_msg=f"unexpected value")
    
    def test_ces(self, ensemble_aligned):
        dm = DistanceMatrix(ensemble_aligned, select="name CA").run()
        rmsd_mat = dm.dist_matrix
        clusters = Clusters(methods.AffinityPropagation())
        clusters.run(-rmsd_mat)
        result_value = reencore.ces(ensemble_aligned, clusters)[0, 1]
        assert_almost_equal(result_value, 0.51, decimal=2,
                            err_msg=f"unexpected value")

    def test_ces_estimate_error(self, ensemble1_aligned, dist_mat1):
        omean, ostd = encore.ces(ensemble1_aligned.universes,
                                 estimate_error=True,
                                 bootstrapping_samples=10,
                                 clustering_method=encore.AffinityPropagationNative(preference=-2.0),
                                 select="name CA and resnum 1-10")
        rmsd_mat_enc = dist_mat1.as_array()
        clusters = Clusters(methods.AffinityPropagation(preference=-2.0))
        clusters.run(-rmsd_mat_enc)
        mean, std = reencore.ces(ensemble1_aligned, clusters, estimate_error=True,
                                 select="name CA and resnum 1-10",
                                 n_bootstrap_samples=10)
        assert_almost_equal(mean[0, 1], mean[0, 1])
        assert_almost_equal(std[0, 1], std[0, 1])

    
    def test_dres_to_self(self, ensemble1_aligned):
        result = reencore.dres(ensemble1_aligned, pca.PCA, n_components=3)[0, 1]
        assert_almost_equal(result, 0)
    
    def test_dres(self, ensemble_aligned, dist_mat):
        spe = drm.StochasticProximityEmbeddingNative(dimension=3,
                                                     distance_cutoff=1.5,
                                                     min_lam=0.1,
                                                     max_lam=2.0,
                                                     ncycle=100,
                                                     nstep=10000)
        dimred, _ = spe(dist_mat)
        old, _ = encore.dres(ensemble_aligned.universes)
        new = reencore.dres(ensemble_aligned, dimred.T, seed=0)
        assert_almost_equal(old[0, 1], new[0, 1], decimal=1)

    def test_dres_estimate_error(self, ensemble1_aligned):
        omean, ostd = encore.dres(ensemble1_aligned.universes,
                                 estimate_error=True,
                                 bootstrapping_samples=10,
                                 select="name CA and resnum 1-10")
        mean, std = reencore.dres(ensemble1_aligned, estimate_error=True,
                                  select="name CA and resnum 1-10",
                                  n_bootstrap_samples=10)
        assert_almost_equal(mean[0, 1], mean[0, 1])
        assert_almost_equal(std[0, 1], std[0, 1])
    
    def test_ces_convergence(self, ens1_aligned):
        # clusters
        dm = DistanceMatrix(ens1_aligned, select="name CA").run()
        rmsd_mat = dm.dist_matrix
        clusters = Clusters(methods.AffinityPropagation())
        clusters.run(-rmsd_mat)

        results = reencore.ces_convergence(ens1_aligned, clusters,
                                           window_size=5)
        exp = np.array([0.34, 0.19, 0.07,  0.])
        assert_almost_equal(results, exp, decimal=2)

    def test_dres_convergence(self, ens1_aligned):
        results = reencore.dres_convergence(ens1_aligned, pca.PCA,
                                            window_size=10,
                                            select="name CA",
                                            n_components=3,
                                            seed=0)
        exp = np.array([0.5, 0.])
        assert_almost_equal(results, exp, decimal=1)

    def test_dres_convergence_spe(self, ens1_aligned):
        dist_mat = encore.get_distance_matrix(ens1_aligned)
        spe = drm.StochasticProximityEmbeddingNative(dimension=3,
                                                     nstep=10000)
        dimred, _ = spe(dist_mat)
        results = reencore.dres_convergence(ens1_aligned, dimred.T,
                                            window_size=10)
        exp = np.array([0.3, 0.])
        assert_almost_equal(results, exp, decimal=1)





