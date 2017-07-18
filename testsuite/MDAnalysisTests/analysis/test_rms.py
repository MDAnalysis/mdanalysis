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
from __future__ import print_function, division, absolute_import

import pytest
from _pytest.capture import capsys
from six.moves import range

import warnings
import os
import sys

import MDAnalysis
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

from numpy.testing import (TestCase, dec, assert_equal,
                           assert_almost_equal, raises, assert_,
                           assert_array_almost_equal)

import numpy as np

from MDAnalysis.exceptions import SelectionError, NoDataError
from MDAnalysisTests.datafiles import GRO, XTC, rmsfArray, PSF, DCD
from MDAnalysisTests import tempdir

# I want to catch all warnings in the tests. If this is not set at the start it
# could cause test that check for warnings to fail.
warnings.simplefilter('always')


class Testrmsd(TestCase):
    def setUp(self):
        shape = (5, 3)
        # vectors with length one
        ones = np.ones(shape) / np.sqrt(3)
        self.a = ones * np.arange(1, 6)[:, np.newaxis]
        self.b = self.a + ones

        self.u = mda.Universe(PSF, DCD)
        self.u2 = mda.Universe(PSF, DCD)

        self.p_first = self.u.select_atoms('protein')
        self.p_last = self.u2.select_atoms('protein')

        self.u.trajectory[2]
        self.u2.trajectory[-2]
        # reset coordinates
        self.u.trajectory[0]
        self.u2.trajectory[-1]

    def test_no_center(self):
        rmsd = rms.rmsd(self.a, self.b, center=False)
        assert_almost_equal(rmsd, 1.0)

    def test_center(self):
        rmsd = rms.rmsd(self.a, self.b, center=True)
        assert_almost_equal(rmsd, 0.0)

    def test_list(self):
        rmsd = rms.rmsd(self.a.tolist(),
                        self.b.tolist(),
                        center=False)
        assert_almost_equal(rmsd, 1.0)

    def test_superposition(self):
        bb = self.u.atoms.select_atoms('backbone')
        a = bb.positions.copy()
        self.u.trajectory[-1]
        b = bb.positions.copy()
        rmsd = rms.rmsd(a, b, superposition=True)
        assert_almost_equal(rmsd, 6.820321761927005)

    def test_weights(self):
        weights = np.zeros(len(self.a))
        weights[0] = 1
        weights[1] = 1
        weighted = rms.rmsd(self.a, self.b, weights=weights)
        firstCoords = rms.rmsd(self.a[:2], self.b[:2])
        assert_almost_equal(weighted, firstCoords)

    def test_weights_and_superposition_1(self):
        weights = np.ones(len(self.u.trajectory[0]))
        weighted = rms.rmsd(self.u.trajectory[0], self.u.trajectory[1],
                            weights=weights, superposition=True)
        firstCoords = rms.rmsd(self.u.trajectory[0], self.u.trajectory[1],
                               superposition=True)
        assert_almost_equal(weighted, firstCoords, decimal=5)

    def test_weights_and_superposition_2(self):
        weights = np.zeros(len(self.u.trajectory[0]))
        weights[:100] = 1
        weighted = rms.rmsd(self.u.trajectory[0], self.u.trajectory[-1],
                            weights=weights, superposition=True)
        firstCoords = rms.rmsd(self.u.trajectory[0][:100],
                               self.u.trajectory[-1][:100],
                               superposition=True)
        # very close to zero, change significant decimal places to 5
        assert_almost_equal(weighted, firstCoords, decimal=5)

    @staticmethod
    @raises(ValueError)
    def test_unequal_shape():
        a = np.ones((4, 3))
        b = np.ones((5, 3))
        rms.rmsd(a, b)

    @raises(ValueError)
    def test_wrong_weights(self):
        w = np.ones(2)
        rms.rmsd(self.a, self.b, w)

    def test_with_superposition_smaller(self):
        A = self.p_first.positions
        B = self.p_last.positions
        rmsd = rms.rmsd(A, B)
        rmsd_superposition = rms.rmsd(A, B, center=True, superposition=True)
        print(rmsd, rmsd_superposition)
        # by design the super positioned rmsd is smaller
        assert_(rmsd > rmsd_superposition)

    def test_with_superposition_equal(self):
        align.alignto(self.p_first, self.p_last)
        A = self.p_first.positions
        B = self.p_last.positions
        rmsd = rms.rmsd(A, B)
        rmsd_superposition = rms.rmsd(A, B, center=True, superposition=True)
        assert_almost_equal(rmsd, rmsd_superposition)


class TestRMSD(object):
    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(PSF, DCD)

    @pytest.fixture()
    def tempdir(self):
        return tempdir.TempDir()

    @pytest.fixture()
    def outfile(self, tempdir):
        return os.path.join(tempdir.name, 'rmsd.txt')

    @pytest.fixture()
    def correct_values(self):
        return [[0, 0, 0], [49, 48.9999, 4.68953]]

    @pytest.fixture()
    def correct_values_group(self):
        return [[0, 0, 0, 0, 0],
                 [49, 49, 4.7857, 4.7004,
                  4.68981]]

    def test_progress_meter(self, capsys, universe):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, verbose=True)
        RMSD.run()
        out, err = capsys.readouterr()
        expected = 'RMSD  6.93 A at frame    98/98  [100.0%]'
        actual = err.strip().split('\r')[-1]
        assert_equal(actual, expected)

    def test_rmsd(self, universe, correct_values):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, select='name CA',
                                            step=49)
        RMSD.run()
        assert_array_almost_equal(RMSD.rmsd, correct_values, 4,
                                  err_msg="error: rmsd profile should match" +
                                  "test values")

    def test_rmsd_single_frame(self, universe):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, select='name CA',
                                            start=5, stop=6).run()
        single_frame = [[5, 5, 0.91544906]]
        assert_array_almost_equal(RMSD.rmsd, single_frame, 4,
                                  err_msg="error: rmsd profile should match" +
                                  "test values")

    def test_mass_weighted_and_save(self, universe, outfile):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, select='name CA',
                                            step=49, weights='mass').run()
        RMSD.save(outfile)
        saved = np.loadtxt(outfile)
        assert_array_almost_equal(RMSD.rmsd, saved, 4,
                                  err_msg="error: rmsd profile should match " +
                                          "test values")

    def test_custom_weighted_and_save(self, universe, outfile):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, select='name CA',
                                            step=49, weights=universe.atoms.CA.masses).run()
        RMSD.save(outfile)
        saved = np.loadtxt(outfile)
        assert_array_almost_equal(RMSD.rmsd, saved, 4,
                                  err_msg="error: rmsd profile should match " +
                                          "test values")

    def test_mass_weighted_and_save_deprecated(self, universe, outfile):
        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter('always')
            RMSD = MDAnalysis.analysis.rms.RMSD(universe, select='name CA',
                                                step=49, mass_weighted=True).run()
        assert_equal(len(warn), 1)
        RMSD.save(outfile)
        saved = np.loadtxt(outfile)
        assert_array_almost_equal(RMSD.rmsd, saved, 4,
                                  err_msg="error: rmsd profile should match " +
                                          "test values")

    def test_rmsd_group_selections(self, universe, correct_values_group):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                            groupselections=
                                            ['backbone', 'name CA'],
                                            step=49).run()
        assert_array_almost_equal(RMSD.rmsd, correct_values_group, 4,
                                  err_msg="error: rmsd profile should match" +
                                  "test values")

    def test_ref_length_unequal_len(self, universe):
        reference = MDAnalysis.Universe(PSF, DCD)
        reference.atoms = reference.atoms[:-1]
        with pytest.raises(SelectionError):
            RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                                reference=reference)

    def test_mass_mismatches(self, universe):
        reference = MDAnalysis.Universe(PSF, DCD)
        reference.atoms.masses = 10
        with pytest.raises(SelectionError):
            RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                                reference=reference)

    def test_group_selections_unequal_len(self, universe):
        reference = MDAnalysis.Universe(PSF, DCD)
        reference.atoms[0].residue.resname='NOTMET'
        with pytest.raises(SelectionError):
            RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                                reference=reference,
                                                groupselections=
                                                ['resname MET','type NH3'])

    def test_save_before_run(self, universe):
        with pytest.raises(NoDataError):
            RMSD = MDAnalysis.analysis.rms.RMSD(universe)
            RMSD.save('blah')


class TestRMSF(TestCase):
    def setUp(self):
        self.universe = mda.Universe(GRO, XTC)
        self.tempdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tempdir.name, 'rmsf.xtc')

    def tearDown(self):
        del self.universe
        del self.tempdir

    def test_rmsf(self):
        rmsfs = rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run()
        test_rmsfs = np.load(rmsfArray)

        assert_almost_equal(rmsfs.rmsf, test_rmsfs, 5,
                            err_msg="error: rmsf profile should match test "
                            "values")

    def test_rmsf_single_frame(self):
        rmsfs = rms.RMSF(self.universe.select_atoms('name CA'), start=5, stop=6).run()

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be zero")

    def test_rmsf_old_run(self):
        rmsfs = rms.RMSF(self.universe.select_atoms('name CA')).run(start=5, stop=6)

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be zero")

    def test_rmsf_identical_frames(self):
        # write a dummy trajectory of all the same frame
        with mda.Writer(self.outfile, self.universe.atoms.n_atoms) as W:
            for _ in range(self.universe.trajectory.n_frames):
                W.write(self.universe)

        self.universe = mda.Universe(GRO, self.outfile)
        rmsfs = rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run()
        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be 0")

    def test_rmsf_weights(self):
        CA = self.universe.atoms.CA
        rmsfs_mass = rms.RMSF(CA, weights='mass')
        rmsfs_mass.run()

        rmsfs_custom = rms.RMSF(CA, weights=CA.masses)
        rmsfs_custom.run()

        assert_almost_equal(rmsfs_mass.rmsf, rmsfs_custom.rmsf, 5,
                            err_msg="error: rmsf profile should match test "
                            "values")
