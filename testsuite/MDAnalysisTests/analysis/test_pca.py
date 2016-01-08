#!/usr/bin/env python2.7

import MDAnalysis
import MDAnalysis.analysis.pca as pca
import numpy as np

from MDAnalysisTests import module_not_found
from numpy.testing import TestCase, assert_almost_equal, assert_equal, dec

from MDAnalysisTests.datafiles import GRO, XTC, evalsArray, evecsArray, scoresArray


class TestPca(TestCase):
    @dec.skipif(
        module_not_found('scipy'),
        "Test skipped because scipy is not available.")
    def setUp(self):
        self.universe = MDAnalysis.Universe(GRO, XTC)

    def tearDown(self):
        del self.universe

    def test_pca_numframes_1(self):
        p = pca.PCA(self.universe, atom_selection='name CA')
        assert_equal(p.numframes, 10)

    def test_pca_numframes_2(self):
        p = pca.PCA(self.universe, atom_selection='name CA', start=2)
        assert_equal(p.numframes, 8)

    def test_pca_natoms(self):
        p = pca.PCA(self.universe, atom_selection='name CA')
        assert_equal(p.natoms, 214)

    def test_pca_nvecs_1(self):
        p = pca.PCA(self.universe, atom_selection='name CA')
        assert_equal(p.nvecs, 642)

    def test_pca_nvecs_2(self):
        p = pca.PCA(self.universe,
                    atom_selection='name CA',
                    captured_variance=0.95)
        assert_equal(p.nvecs, 6)

    def test_pca_nvecs_3(self):
        p = pca.PCA(self.universe,
                    atom_selection='name CA',
                    captured_variance=4)
        assert_equal(p.nvecs, 4)

    def test_pca_evals(self):
        p = pca.PCA(self.universe,
                    atom_selection='name CA',
                    captured_variance=4)
        test_evals = np.load(evalsArray)
        assert_almost_equal(
            p.evals,
            test_evals,
            5,
            err_msg="error - eigenvalues do not match reference data")

    def test_pca_evs(self):
        p = pca.PCA(self.universe,
                    atom_selection='name CA',
                    captured_variance=4)
        test_evecs = np.load(evecsArray)
        assert_almost_equal(
            p.evecs[0, :4],
            test_evecs,
            5,
            err_msg="error - eigenvector does not match reference data")

    def test_pca_scores(self):
        p = pca.PCA(self.universe,
                    atom_selection='name CA',
                    captured_variance=4)
        test_scores = np.load(scoresArray)
        assert_almost_equal(
            p.scores(7),
            test_scores,
            5,
            err_msg="error - scores do not match reference data")

    def test_pca_transform(self):
        p = pca.PCA(self.universe, atom_selection='name CA')
        s = p.scores(5)
        x = p.transform(s, inverse=True)
        x2 = p._f.coords(5)
        assert_almost_equal(
            pca.rmsd(x, x2),
            0.0,
            4,
            err_msg=
            "error - transform/inverse transform sequence does not regenerate original coordinates")
