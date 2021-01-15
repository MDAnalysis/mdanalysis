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
        a = align.AlignTraj(ens, ens).run()
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

    def test_ces_to_self(self, ensemble1_aligned):
        dm = DistanceMatrix(ensemble1_aligned, select="name CA").run()
        rmsd_mat1 = dm.dist_matrix
        dm = DistanceMatrix(ensemble1_aligned).run()
        clusters = Clusters(methods.AffinityPropagation(preference=-3.0))
        clusters.run(-rmsd_mat1)
        result_value = reencore.ces(ensemble1_aligned, clusters)[0, 1]
        assert_almost_equal(result_value, 0,
                            err_msg=f"ces() to itself not zero: {result_value}")

    def test_ces_rmsd_enc(self, ensemble_aligned):
        rmsd_mat_enc = encore.get_distance_matrix(ensemble_aligned).as_array()
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
        assert_almost_equal(result_value, 0.69, decimal=2,
                            err_msg=f"unexpected value")
    
    def test_dres_to_self(self, ensemble1_aligned):
        result = reencore.dres(ensemble1_aligned)[0, 1]
        assert_almost_equal(result, 0)
        


