# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
from __future__ import print_function

from six.moves import range

import MDAnalysis
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

from numpy.testing import (TestCase, dec, assert_equal,
                           assert_almost_equal, raises, assert_,
                           assert_array_almost_equal)

import numpy as np

import os
import sys

from MDAnalysis.exceptions import SelectionError, NoDataError
from MDAnalysisTests.datafiles import GRO, XTC, rmsfArray, PSF, DCD
from MDAnalysisTests import tempdir, parser_not_found


class Testrmsd(object):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def __init__(self):
        shape = (5, 3)
        # vectors with length one
        ones = np.ones(shape) / np.sqrt(3)
        self.a = ones * np.arange(1, 6)[:, np.newaxis]
        self.b = self.a + ones

        self.u = mda.Universe(PSF, DCD)
        self.u2 = mda.Universe(PSF, DCD)

        self.p_first = self.u.select_atoms('protein')
        self.p_last = self.u2.select_atoms('protein')

    def setUp(self):
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
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.tempdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tempdir.name, 'rmsd.txt')
        self.correct_values = [[0, 0, 0], [49, 48.9999, 4.68953]]
        self.correct_values_group = [[0, 0, 0, 0, 0],
                                     [49, 48.9999, 4.7857, 4.7002,
                                      4.68981]]

    def tearDown(self):
        del self.universe
        del self.tempdir

    def test_progress_meter(self):
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe, verbose=True)
        sys.stderr = sys.stdout
        RMSD.run()
        actual = sys.stderr.getvalue().strip().split('\r')[-1]
        expected = 'RMSD  6.93 A at frame    98/98  [100.0%]'
        sys.stderr = sys.__stderr__
        assert_equal(actual, expected)

    def test_rmsd(self):
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe, select='name CA',
                                            step=49)
        RMSD.run()
        assert_array_almost_equal(RMSD.rmsd, self.correct_values, 4,
                                  err_msg="error: rmsd profile should match" +
                                  "test values")

    def test_rmsd_single_frame(self):
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe, select='name CA',
                                            start=5, stop=6).run()
        single_frame = [[5, 5, 0.91544906]]
        assert_array_almost_equal(RMSD.rmsd, single_frame, 4,
                                  err_msg="error: rmsd profile should match" +
                                  "test values")

    def test_mass_weighted_and_save(self):
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe, select='name CA',
                                            step=49, weights='mass').run()
        RMSD.save(self.outfile)
        saved = np.loadtxt(self.outfile)
        assert_array_almost_equal(RMSD.rmsd, saved, 4,
                                  err_msg="error: rmsd profile should match " +
                                          "test values")

    def test_rmsd_group_selections(self):
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe,
                                            groupselections=
                                            ['backbone', 'name CA'],
                                            step=49).run()
        assert_array_almost_equal(RMSD.rmsd, self.correct_values_group, 4,
                                  err_msg="error: rmsd profile should match" +
                                  "test values")

    @raises(SelectionError)
    def test_ref_length_unequal_len(self):
        reference = MDAnalysis.Universe(PSF, DCD)
        reference.atoms = reference.atoms[:-1]
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe,
                                            reference=reference)

    @raises(SelectionError)
    def test_mass_mismatches(self):
        reference = MDAnalysis.Universe(PSF, DCD)
        reference.atoms.masses = 10
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe,
                                            reference=reference)


    @raises(SelectionError)
    def test_group_selections_unequal_len(self):
        reference = MDAnalysis.Universe(PSF, DCD)
        reference.atoms[0].residue.resname='NOTMET'
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe,
                                            reference=reference,
                                            groupselections=
                                            ['resname MET','type NH3'])
    @raises(NoDataError)
    def test_save_before_run(self):
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe)
        RMSD.save('blah')

class TestRMSF(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(GRO, XTC)
        self.tempdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tempdir.name, 'rmsf.xtc')

    def tearDown(self):
        del self.universe
        del self.tempdir

    def test_rmsf(self):
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(verbose=False)
        test_rmsfs = np.load(rmsfArray)

        assert_almost_equal(rmsfs.rmsf, test_rmsfs, 5,
                            err_msg="error: rmsf profile should match test " +
                            "values")

    def test_rmsf_single_frame(self):
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(start=5, stop=6, verbose=False)

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be zero")

    def test_rmsf_identical_frames(self):
        # write a dummy trajectory of all the same frame
        with MDAnalysis.Writer(self.outfile, self.universe.atoms.n_atoms) as W:
            for i in range(self.universe.trajectory.n_frames):
                W.write(self.universe)

        self.universe = MDAnalysis.Universe(GRO, self.outfile)
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(verbose=False)

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be 0")
