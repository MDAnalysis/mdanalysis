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

import MDAnalysis.analysis.encore as encore

from numpy.testing import (TestCase, dec, assert_equal, assert_almost_equal)

from MDAnalysisTests.datafiles import DCD, DCD2, PDB_small
from MDAnalysisTests import parser_not_found

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
        total_frames = len(self.ens1.get_coordinates())
        filtered_ensemble = encore.Ensemble(topology=PDB_small, trajectory=DCD,
                                            frame_interval=10)
        filtered_frames = len(filtered_ensemble.get_coordinates())
        assert_equal(filtered_frames, total_frames//10,
                     err_msg="Incorrect frame number in Ensemble filtering: {0:f} out of {1:f}"
                     .format(filtered_frames, total_frames))

    def test_ensemble_atom_selection_default(self):
        coordinates_per_frame_default = len(self.ens1.get_coordinates()[0])
        expected_value = 214
        assert_equal(coordinates_per_frame_default, expected_value,
                     err_msg="Unexpected value for Harmonic Ensemble Similarity: {0:f}. "
                             "Expected {1:f}.".format(coordinates_per_frame_default, expected_value))

    def test_ensemble_atom_selection_full(self):
        ensemble_full = encore.Ensemble(topology=PDB_small, trajectory=DCD, atom_selection_string="name *")
        coordinates_per_frame_full = len(ensemble_full.get_coordinates()[0])
        expected_value = 3341
        assert_equal(coordinates_per_frame_full, expected_value,
                     err_msg="Unexpected value for Harmonic Ensemble Similarity: {0:f}. "
                             "Expected {1:f}.".format(coordinates_per_frame_full, expected_value))

    @dec.slow
    def test_hes_to_self(self):
        results, details = encore.hes([self.ens1, self.ens1])
        result_value = results[0,1]
        expected_value = 0.
        assert_almost_equal(results[0, 1], expected_value,
                            err_msg="Harmonic Ensemble Similarity to itself not zero: {0:f}".format(result_value))

    @dec.slow
    def test_hes(self):
        results, details = encore.hes([self.ens1, self.ens2])
        result_value = results[0, 1]
        expected_value = 13946090.576
        assert_almost_equal(results[0, 1], expected_value, decimal=2,
                            err_msg="Unexpected value for Harmonic Ensemble Similarity: {0:f}. "
                                    "Expected {1:f}.".format(result_value, expected_value))

    @dec.slow
    def test_ces_to_self(self):
        results, details = encore.ces([self.ens1, self.ens1])
        result_value = results[0,0,1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value,
                            err_msg="ClusteringEnsemble Similarity to itself not zero: {0:f}".format(result_value))

    @dec.slow
    def test_ces(self):
        results, details = encore.ces([self.ens1, self.ens2])
        result_value = results[0,0,1]
        expected_value = 0.55392
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Unexpected value for Cluster Ensemble Similarity: {}. Expected {}.".format(result_value, expected_value))

    @dec.slow
    def test_dres_to_self(self):
        results, details = encore.dres([self.ens1, self.ens1])
        result_value = results[0,0,1]
        expected_value = 0.
        assert_almost_equal(result_value, expected_value, decimal=2,
                            err_msg="Dim. Reduction Ensemble Similarity to itself not zero: {0:f}"
                            .format(result_value))

    @dec.slow
    def test_dres(self):
        results, details = encore.dres([self.ens1, self.ens2])
        result_value = results[0,0,1]
        expected_value = 0.68
        assert_almost_equal(result_value, expected_value, decimal=1,
                            err_msg="Unexpected value for Dim. reduction Ensemble Similarity: {0:f}. "
                                    "Expected {1:f}.".format(result_value, expected_value))

