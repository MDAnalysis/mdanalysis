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

from __future__ import print_function

import MDAnalysis as mda
import MDAnalysis.analysis.encore as encore

import tempfile
import numpy

from numpy.testing import (TestCase, dec, assert_equal, assert_almost_equal)

from MDAnalysisTests.datafiles import DCD, DCD2, PDB_small
from MDAnalysisTests import parser_not_found, module_not_found

import MDAnalysis.analysis.rms as rms
import MDAnalysis.analysis.align as align


class FakePBarCounter(object):
    def __init__(self):
        self.value = 0


class TestEncore(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.ens1 = mda.Universe(PDB_small, DCD)
        self.ens2 = mda.Universe(PDB_small, DCD2)

    def tearDown(self):
        del self.ens1
        del self.ens2

    @staticmethod
    def test_triangular_matrix():
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

    @staticmethod
    def test_parallel_calculation():

        def function(x):
            return x**2

        arguments = [tuple([i]) for i in numpy.arange(0,100)]

        parallel_calculation = encore.utils.ParallelCalculation(function = function,
                                                                ncores = 4,
                                                                args = arguments)
        results = parallel_calculation.run()

        for i,r in enumerate(results):
            assert_equal(r[1], arguments[i][0]**2,
                err_msg="Unexpeted results from ParallelCalculation")



    def test_rmsd_matrix_with_superimposition(self):
        
        generator = encore.confdistmatrix.RMSDMatrixGenerator()
        confdist_matrix = generator(self.ens1,
                                    selection = "name CA",
                                    pairwise_align = True,
                                    mass_weighted = True,
                                    ncores = 1)

        reference = rms.RMSD(self.ens1, select = "name CA")
        reference.run()


        tasks = ((0, 0), (1, 0))
        n_tasks = len(list(encore.utils.trm_indeces(tasks[0],tasks[1])))
        distmatrix = numpy.zeros(n_tasks)
        coordinates = self.ens1.trajectory.timeseries(
            self.ens1.select_atoms("name CA"), format = 'fac')
        masses = numpy.ones(coordinates.shape[1])
        pbar_counter = FakePBarCounter()

        generator._fitter_worker(tasks,
                                 coordinates,
                                 coordinates,
                                 masses,
                                 masses,
                                 distmatrix,
                                 pbar_counter)

        for i in range(n_tasks):
            assert_almost_equal(reference.rmsd[i,2], distmatrix[i], decimal = 3,
                                err_msg = "calculated RMSD values differ from the reference implementation")
                                       
        for i,rmsd in enumerate(reference.rmsd):
            assert_almost_equal(rmsd[2], confdist_matrix[0,i], decimal=3,
                                err_msg = "calculated RMSD values differ from the reference implementation")

    def test_minus_rmsd_matrix_with_superimposition(self):
        
        generator = encore.confdistmatrix.MinusRMSDMatrixGenerator()
        confdist_matrix = generator(self.ens1,
                                    selection = "name CA",
                                    pairwise_align = True,
                                    mass_weighted = True,
                                    ncores = 1)

        reference = rms.RMSD(self.ens1, select = "name CA")
        reference.run()

        for i,rmsd in enumerate(reference.rmsd):
            assert_almost_equal(-rmsd[2], confdist_matrix[0,i], decimal=3,
                                err_msg = "calculated RMSD values differ from the reference implementation")            

    def test_rmsd_matrix_without_superimposition(self):
        
        # calculated with gromacs - gmx rms -fit none
        reference_rmsd =[0.0000001,
                         0.0425684,
                         0.0595158,
                         0.0738680,
                         0.0835519,
                         0.0924640,
                         0.1010487,
                         0.1131771,
                         0.1227527,
                         0.1343707,
                         0.1433841,
                         0.1545489,
                         0.1638420,
                         0.1720007,
                         0.1818408,
                         0.1897694,
                         0.1979185,
                         0.2050228,
                         0.2190710,
                         0.2282337,
                         0.2392368,
                         0.2467754,
                         0.2559295,
                         0.2634292,
                         0.2758299,
                         0.2815295,
                         0.2889598,
                         0.2988116,
                         0.3075704,
                         0.3168339,
                         0.3252532,
                         0.3335701,
                         0.3421980,
                         0.3499905,
                         0.3576347,
                         0.3648850,
                         0.3746280,
                         0.3787407,
                         0.3876532,
                         0.3949267,
                         0.4022163,
                         0.4123725,
                         0.4171653,
                         0.4270313,
                         0.4339235,
                         0.4441433,
                         0.4535998,
                         0.4629753,
                         0.4738565,
                         0.4778692,
                         0.4846473,
                         0.4921997,
                         0.5025109,
                         0.5078515,
                         0.5176530,
                         0.5236758,
                         0.5279259,
                         0.5359889,
                         0.5479882,
                         0.5513062,
                         0.5550882,
                         0.5616842,
                         0.5691664,
                         0.5797819,
                         0.5860255,
                         0.5929349,
                         0.6031308,
                         0.6075997,
                         0.6206015,
                         0.6300921,
                         0.6396201,
                         0.6409384,
                         0.6439900,
                         0.6467734,
                         0.6527478,
                         0.6543783,
                         0.6585453,
                         0.6659292,
                         0.6674148,
                         0.6699741,
                         0.6713669,
                         0.6696672,
                         0.6695362,
                         0.6699672,
                         0.6765218,
                         0.6806746,
                         0.6801361,
                         0.6786651,
                         0.6828524,
                         0.6851146,
                         0.6872993,
                         0.6837722,
                         0.6852713,
                         0.6838173,
                         0.6822636,
                         0.6829022,
                         0.6846855,
                         0.6843332 ]

        selection_string = "name CA"
        generator = encore.confdistmatrix.RMSDMatrixGenerator()
        confdist_matrix = generator(self.ens1,
                                    selection = selection_string,
                                    pairwise_align = False,
                                    mass_weighted = True,
                                    ncores = 1)
        
        for i,rmsd in enumerate(reference_rmsd):
            assert_almost_equal(confdist_matrix[0,i]/10.0, rmsd, decimal=3,
                                err_msg = "calculated RMSD values differ from the reference implementation")

    def test_minus_rmsd_matrix_without_superimposition(self):
        
        # calculated with gromacs - gmx rms -fit none
        reference_rmsd =[0.0000001,
                         0.0425684,
                         0.0595158,
                         0.0738680,
                         0.0835519,
                         0.0924640,
                         0.1010487,
                         0.1131771,
                         0.1227527,
                         0.1343707,
                         0.1433841,
                         0.1545489,
                         0.1638420,
                         0.1720007,
                         0.1818408,
                         0.1897694,
                         0.1979185,
                         0.2050228,
                         0.2190710,
                         0.2282337,
                         0.2392368,
                         0.2467754,
                         0.2559295,
                         0.2634292,
                         0.2758299,
                         0.2815295,
                         0.2889598,
                         0.2988116,
                         0.3075704,
                         0.3168339,
                         0.3252532,
                         0.3335701,
                         0.3421980,
                         0.3499905,
                         0.3576347,
                         0.3648850,
                         0.3746280,
                         0.3787407,
                         0.3876532,
                         0.3949267,
                         0.4022163,
                         0.4123725,
                         0.4171653,
                         0.4270313,
                         0.4339235,
                         0.4441433,
                         0.4535998,
                         0.4629753,
                         0.4738565,
                         0.4778692,
                         0.4846473,
                         0.4921997,
                         0.5025109,
                         0.5078515,
                         0.5176530,
                         0.5236758,
                         0.5279259,
                         0.5359889,
                         0.5479882,
                         0.5513062,
                         0.5550882,
                         0.5616842,
                         0.5691664,
                         0.5797819,
                         0.5860255,
                         0.5929349,
                         0.6031308,
                         0.6075997,
                         0.6206015,
                         0.6300921,
                         0.6396201,
                         0.6409384,
                         0.6439900,
                         0.6467734,
                         0.6527478,
                         0.6543783,
                         0.6585453,
                         0.6659292,
                         0.6674148,
                         0.6699741,
                         0.6713669,
                         0.6696672,
                         0.6695362,
                         0.6699672,
                         0.6765218,
                         0.6806746,
                         0.6801361,
                         0.6786651,
                         0.6828524,
                         0.6851146,
                         0.6872993,
                         0.6837722,
                         0.6852713,
                         0.6838173,
                         0.6822636,
                         0.6829022,
                         0.6846855,
                         0.6843332 ]

        selection_string = "name CA"
        generator = encore.confdistmatrix.MinusRMSDMatrixGenerator()
        confdist_matrix = generator(self.ens1,
                                    selection = selection_string,
                                    pairwise_align = False,
                                    mass_weighted = True,
                                    ncores = 1)
        
        for i,rmsd in enumerate(reference_rmsd):
            assert_almost_equal(-confdist_matrix[0,i]/10.0, rmsd, decimal=3,
                                err_msg = "calculated RMSD values differ from the reference implementation")

        
    def test_ensemble_frame_filtering(self):
        total_frames = len(self.ens1.trajectory.timeseries(format='fac'))
        interval = 10
        filtered_ensemble = mda.Universe(PDB_small, DCD,
                                         in_memory=True,
                                         in_memory_frame_interval=interval)
        filtered_frames = len(filtered_ensemble.trajectory.timeseries(format='fac'))
        assert_equal(filtered_frames, total_frames//interval,
                     err_msg="Incorrect frame number in Ensemble filtering: {0:f} out of {1:f}"
                     .format(filtered_frames, total_frames//interval))

    def test_ensemble_atom_selection_default(self):
        coordinates_per_frame_default = len(self.ens1.atoms.coordinates())
        expected_value = 3341
        assert_equal(coordinates_per_frame_default, expected_value,
                     err_msg="Unexpected atom number in default selection: {0:f}. "
                             "Expected {1:f}.".format(coordinates_per_frame_default, expected_value))

    @staticmethod
    def test_ensemble_superimposition():
        aligned_ensemble1 = mda.Universe(PDB_small, DCD)
        align.rms_fit_trj(aligned_ensemble1, aligned_ensemble1,
                          select="name CA",
                          in_memory=True)
        aligned_ensemble2 = mda.Universe(PDB_small, DCD)
        align.rms_fit_trj(aligned_ensemble2, aligned_ensemble2,
                          select="name *",
                          in_memory=True)

        rmsfs1 = rms.RMSF(aligned_ensemble1.select_atoms('name *'))
        rmsfs1.run()

        rmsfs2 = rms.RMSF(aligned_ensemble2.select_atoms('name *'))
        rmsfs2.run()

        assert_equal(sum(rmsfs1.rmsf)>sum(rmsfs2.rmsf), True,
                     err_msg="Ensemble aligned on all atoms should have lower full-atom RMSF "
                             "than ensemble aligned on only CAs.")

    @staticmethod
    def test_ensemble_superimposition_to_reference_non_weighted():
        ensemble0 = mda.Universe(PDB_small, DCD)
        filename = align.rms_fit_trj(ensemble0, ensemble0,
                                     select="name CA", mass_weighted=False)
        aligned_ensemble0 = mda.Universe(PDB_small, filename)
        aligned_ensemble1 = mda.Universe(PDB_small, DCD)
        align.rms_fit_trj(aligned_ensemble1, aligned_ensemble1,
                          select="name CA", mass_weighted=False,
                          in_memory=True)
        aligned_ensemble2 = mda.Universe(PDB_small, DCD)
        align.rms_fit_trj(aligned_ensemble2, aligned_ensemble2,
                          select="name *", mass_weighted=False,
                          in_memory=True)

        rmsfs0 = rms.RMSF(aligned_ensemble0.select_atoms('name *'))
        rmsfs0.run()

        rmsfs1 = rms.RMSF(aligned_ensemble1.select_atoms('name *'))
        rmsfs1.run()

        rmsfs2 = rms.RMSF(aligned_ensemble2.select_atoms('name *'))
        rmsfs2.run()

        import logging
        logging.info("{0} {1} {2}".format(sum(rmsfs1.rmsf), sum(rmsfs2.rmsf), sum(rmsfs0.rmsf)))

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
    def test_hes_align(self):
        results, details = encore.hes([self.ens1, self.ens2], align=True)
        result_value = results[0,1]
        expected_value = 6868.28
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Unexpected value for Harmonic Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))
    @dec.slow
    def test_hes_ml_cov(self):
        results, details = encore.hes([self.ens1, self.ens2], cov_estimator="ml")
        result_value = results[0,1]
        expected_value = 50687.12
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Unexpected value for Harmonic Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))

    @dec.slow
    def test_ces_to_self(self):
        results, details = encore.ces([self.ens1, self.ens1], preference_values = -3.0)
        result_value = results[0,1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value,
                            err_msg="ClusteringEnsemble Similarity to itself not zero: {0:f}".format(result_value))

    @dec.slow
    def test_ces(self):
        results, details = encore.ces([self.ens1, self.ens2])
        result_value = results[0,1]
        expected_value = 0.68070
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Unexpected value for Cluster Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))
        
    @dec.slow
    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def test_dres_to_self(self):
        results, details = encore.dres([self.ens1, self.ens1])
        result_value = results[0,1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Dim. Reduction Ensemble Similarity to itself not zero: {0:f}".format(result_value))

    @dec.slow
    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def test_dres(self):
        results, details = encore.dres([self.ens1, self.ens2])
        result_value = results[0,1]
        expected_value = 0.68
        assert_almost_equal(result_value, expected_value, decimal=1,
                            err_msg="Unexpected value for Dim. reduction Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))

    @dec.slow
    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def test_dres_without_superimposition(self):
        results, details = encore.dres([self.ens1, self.ens2], superimpose=False)
        result_value = results[0,1]
        expected_value = 0.68
        assert_almost_equal(result_value, expected_value, decimal=1,
                            err_msg="Unexpected value for Dim. reduction Ensemble Similarity: {0:f}. Expected {1:f}.".format(result_value, expected_value))
        
    @dec.slow
    def test_ces_convergence(self):
        expected_values = [ 0.48194205,  0.40284672,  0.31699026,  0.25220447,  0.19829817,
         0.14642725,  0.09911411,  0.05667391,  0.        ]
        results = encore.ces_convergence(self.ens1, 10)
        for i,ev in enumerate(expected_values):
            assert_almost_equal(ev, results[i], decimal=2, 
                                err_msg="Unexpected value for Clustering Ensemble similarity in convergence estimation")
    @dec.slow
    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def test_dres_convergence(self):
        expected_values = [ 0.53998088,  0.40466411,  0.30709079,  0.26811765,  0.19571984,
         0.11489109,  0.06484937,  0.02803273,  0.        ]
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
        averages, stdevs = encore.ces([self.ens1, self.ens1], estimate_error = True, bootstrapping_samples=10, preference_values=-2.0)
        average = averages[0,1]
        stdev = stdevs[0,1]

        assert_almost_equal(expected_average, average, decimal=1, 
                            err_msg="Unexpected average value for bootstrapped samples in Clustering Ensemble similarity")
        assert_almost_equal(expected_average, average, decimal=1, 
                            err_msg="Unexpected standard daviation  for bootstrapped samples in Clustering Ensemble similarity")        

    @dec.slow
    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
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
