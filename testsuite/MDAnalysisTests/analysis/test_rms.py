# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
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

from six.moves import range

import MDAnalysis
import MDAnalysis.analysis.rms

from numpy.testing import TestCase, assert_almost_equal, assert_equal, raises
import numpy as np

import os
import tempdir

from MDAnalysisTests.datafiles import GRO, XTC, rmsfArray, PSF, DCD


class TestRMSD(object):
    def __init__(self):
        shape = (5, 3)
        self.a = np.arange(np.prod(shape)).reshape(shape)
        self.b = np.arange(np.prod(shape)).reshape(shape) + 1

    def test_no_center(self):
        rmsd = MDAnalysis.analysis.rms.rmsd(self.a, self.b, center=False)
        assert_equal(rmsd, 1.0)

    def test_center(self):
        rmsd = MDAnalysis.analysis.rms.rmsd(self.a, self.b, center=True)
        assert_equal(rmsd, 0.0)

    @staticmethod
    def test_list():
        a = [[0, 1, 2],
             [3, 4, 5]]
        b = [[1, 2, 3],
             [4, 5, 6]]
        rmsd = MDAnalysis.analysis.rms.rmsd(a, b, center=False)
        assert_equal(rmsd, 1.0)

    @staticmethod
    def test_superposition():
        u = MDAnalysis.Universe(PSF, DCD)
        bb = u.atoms.select_atoms('backbone')
        a = bb.positions.copy()
        u.trajectory[-1]
        b = bb.positions.copy()
        rmsd = MDAnalysis.analysis.rms.rmsd(a, b, superposition=True)
        assert_almost_equal(rmsd, 6.820321761927005)

    @staticmethod
    @raises(ValueError)
    def test_unequal_shape():
        a = np.ones((4, 3))
        b = np.ones((5, 3))
        MDAnalysis.analysis.rms.rmsd(a, b)

    @raises(ValueError)
    def test_wrong_weights(self):
        w = np.ones(2)
        MDAnalysis.analysis.rms.rmsd(self.a, self.b, w)


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
        rmsfs.run(quiet=True)
        test_rmsfs = np.load(rmsfArray)

        assert_almost_equal(rmsfs.rmsf, test_rmsfs, 5,
                            err_msg="error: rmsf profile should match test " +
                            "values")

    def test_rmsf_single_frame(self):
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(start=5, stop=6, quiet=True)

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be zero")

    def test_rmsf_identical_frames(self):
        # write a dummy trajectory of all the same frame
        with MDAnalysis.Writer(self.outfile, self.universe.atoms.n_atoms) as W:
            for i in range(self.universe.trajectory.n_frames):
                W.write(self.universe)

        self.universe = MDAnalysis.Universe(GRO, self.outfile)
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(quiet=True)

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be 0")
