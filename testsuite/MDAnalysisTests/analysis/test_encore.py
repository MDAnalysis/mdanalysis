# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function

import MDAnalysis as mda
import MDAnalysis.analysis.encore as encore

from numpy.testing import (TestCase, dec, assert_equal, assert_almost_equal)

from MDAnalysisTests.datafiles import DCD, DCD2, PDB_small, PDB,XTC
from MDAnalysisTests import parser_not_found

import MDAnalysis.analysis.rms as rms

class TestEnsemble(TestCase):

    def test_from_reader_w_timeseries(self):
        ensemble = encore.Ensemble(topology=PDB_small, trajectory=DCD)
        assert_equal(len(ensemble.atoms.coordinates()), 3341,
                     err_msg="Unexpected number of atoms in trajectory")

    def test_from_reader_wo_timeseries(self):
        ensemble = encore.Ensemble(topology=PDB, trajectory=XTC)
        assert_equal(len(ensemble.atoms.coordinates()), 47681,
                     err_msg="Unexpected number of atoms in trajectory")

    def test_trajectories_list(self):
        ensemble = encore.Ensemble(topology=PDB_small, trajectory=[DCD])
        assert_equal(len(ensemble.atoms.coordinates()), 3341,
                     err_msg="Unexpected number of atoms in trajectory")

class TestEncore(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.ens1 = encore.Ensemble(topology=PDB_small, trajectory=DCD)
        self.ens2 = encore.Ensemble(topology=PDB_small, trajectory=DCD2)

    def tearDown(self):
        del self.ens1
        del self.ens2

    def test_ensemble_frame_filtering(self):
        total_frames = len(self.ens1.get_coordinates("", format='fac'))
        interval = 10
        filtered_ensemble = encore.Ensemble(topology=PDB_small, trajectory=DCD,
                                            frame_interval=interval)
        filtered_frames = len(filtered_ensemble.get_coordinates("", format='fac'))
        assert_equal(filtered_frames, total_frames//interval,
                     err_msg="Incorrect frame number in Ensemble filtering: {0:f} out of {1:f}"
                     .format(filtered_frames, total_frames//interval))

    def test_ensemble_atom_selection_default(self):
        coordinates_per_frame_default = len(self.ens1.atoms.coordinates())
        expected_value = 3341
        assert_equal(coordinates_per_frame_default, expected_value,
                     err_msg="Unexpected atom number in default selection: {0:f}. "
                             "Expected {1:f}.".format(coordinates_per_frame_default, expected_value))

    def test_ensemble_superimposition(self):
        aligned_ensemble1 = encore.Ensemble(topology=PDB_small, trajectory=DCD)
        aligned_ensemble1.align(selection="name CA")
        aligned_ensemble2 = encore.Ensemble(topology=PDB_small, trajectory=DCD)
        aligned_ensemble2.align(selection="name *")

        rmsfs1 = rms.RMSF(aligned_ensemble1.select_atoms('name *'))
        rmsfs1.run()

        rmsfs2 = rms.RMSF(aligned_ensemble2.select_atoms('name *'))
        rmsfs2.run()

        assert_equal(sum(rmsfs1.rmsf)>sum(rmsfs2.rmsf), True,
                     err_msg="Ensemble aligned on all atoms should have lower full-atom RMSF "
                             "than ensemble aligned on only CAs.")

    def test_ensemble_superimposition_to_reference_non_weighted(self):
        aligned_ensemble1 = encore.Ensemble(topology=PDB_small, trajectory=DCD)
        aligned_ensemble1.align(selection="name CA", weighted=False,
                                reference=mda.Universe(PDB_small))
        aligned_ensemble2 = encore.Ensemble(topology=PDB_small, trajectory=DCD)
        aligned_ensemble2.align(selection="name *", weighted=False,
                                reference=mda.Universe(PDB_small))

        rmsfs1 = rms.RMSF(aligned_ensemble1.select_atoms('name *'))
        rmsfs1.run()

        rmsfs2 = rms.RMSF(aligned_ensemble2.select_atoms('name *'))
        rmsfs2.run()

        assert_equal(sum(rmsfs1.rmsf)>sum(rmsfs2.rmsf), True,
                     err_msg="Ensemble aligned on all atoms should have lower full-atom RMSF "
                             "than ensemble aligned on only CAs.")

    @dec.slow
    def test_hes_to_self(self):
        results, details = encore.hes([self.ens1, self.ens1])
        result_value = results[0,1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value,
                            err_msg="Harmonic Ensemble Similarity to itself not zero: {0:f}".format(result_value))

    @dec.slow
    def test_hes(self):
        results, details = encore.hes([self.ens1, self.ens2])
        result_value = results[0,1]
        expected_value = 13946090.576
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Unexpected value for Harmonic Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))

    @dec.slow
    def atest_ces_to_self(self):
        results, details = encore.ces([self.ens1, self.ens1])
        result_value = results[0,1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value,
                            err_msg="ClusteringEnsemble Similarity to itself not zero: {0:f}".format(result_value))

    @dec.slow
    def test_ces(self):
        results, details = encore.ces([self.ens1, self.ens2])
        result_value = results[0,1]
        expected_value = 0.55392
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Unexpected value for Cluster Ensemble Similarity: {}. Expected {}.".format(result_value, expected_value))

    @dec.slow
    def test_dres_to_self(self):
        results, details = encore.dres([self.ens1, self.ens1])
        result_value = results[0,1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Dim. Reduction Ensemble Similarity to itself not zero: {0:f}".format(result_value))

    @dec.slow
    def test_dres(self):
        results, details = encore.dres([self.ens1, self.ens2])
        result_value = results[0,1]
        expected_value = 0.68
        assert_almost_equal(result_value, expected_value, decimal=1,
                            err_msg="Unexpected value for Dim. reduction Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))

    @dec.slow
    def test_ces_convergence(self):
        expected_values = [0.51124, 0.38618, 0.28370, 0.26927, 0.19035, 0.12918, 0.08996, 0.06434, 0.00000]
        results = encore.ces_convergence(self.ens1, 10)
        for i,ev in enumerate(expected_values):
            assert_almost_equal(ev, results[i], decimal=2, 
                                err_msg="Unexpected value for Clustering Ensemble similarity in convergence estimation")
    @dec.slow
    def test_dres_convergence(self):
        expected_values = [0.62387, 0.55965, 0.48308, 0.39526, 0.29047, 0.18011, 0.12844, 0.06337, 0.00000]
        #import numpy
        results = encore.dres_convergence(self.ens1, 10)
        for i,ev in enumerate(expected_values):
            assert_almost_equal(ev, results[i], decimal=1, 
                                err_msg="Unexpected value for Dim. reduction Ensemble similarity in convergence estimation")
    @dec.slow
    def test_hes_error_estimation(self):
        expected_average = 0.086
        expected_stdev = 0.009
        averages, stdevs = encore.hes([self.ens1, self.ens1], estimate_error = True, bootstrapping_samples=10)
        average = averages[0,1]
        stdev = stdevs[0,1]

        assert_almost_equal(expected_average, average, decimal=1, 
                            err_msg="Unexpected average value for bootstrapped samples in Harmonic Ensemble imilarity")
        assert_almost_equal(expected_average, average, decimal=1, 
                            err_msg="Unexpected standard daviation  for bootstrapped samples in Harmonic Ensemble imilarity")
    @dec.slow
    def test_ces_error_estimation(self):
        expected_average = 0.02
        expected_stdev = 0.008
        averages, stdevs = encore.ces([self.ens1, self.ens1], estimate_error = True, bootstrapping_samples=10)
        average = averages[0,1]
        stdev = stdevs[0,1]

        assert_almost_equal(expected_average, average, decimal=1, 
                            err_msg="Unexpected average value for bootstrapped samples in Clustering Ensemble similarity")
        assert_almost_equal(expected_average, average, decimal=1, 
                            err_msg="Unexpected standard daviation  for bootstrapped samples in Clustering Ensemble similarity")        
    @dec.slow
    def test_dres_error_estimation(self):
        expected_average = 0.02
        expected_stdev = 0.01
        averages, stdevs = encore.dres([self.ens1, self.ens1], estimate_error = True, bootstrapping_samples=10)
        average = averages[0,1]
        stdev = stdevs[0,1]

        assert_almost_equal(expected_average, average, decimal=1, 
                            err_msg="Unexpected average value for bootstrapped samples in Dim. reduction Ensemble similarity")
        assert_almost_equal(expected_average, average, decimal=1, 
                            err_msg="Unexpected standard daviation for bootstrapped samples in Dim. reduction Ensemble imilarity")        
