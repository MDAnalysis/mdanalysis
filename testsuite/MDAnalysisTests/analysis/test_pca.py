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
import numpy as np
import MDAnalysis
import MDAnalysis.analysis.pca as pca
from numpy.testing import (assert_almost_equal, assert_equal,
                           assert_array_almost_equal,raises)


from MDAnalysisTests.datafiles import PDB, XTC


class TestPCA(object):
    def setUp(self):
        self.u = MDAnalysis.Universe(PDB, XTC)
        self.pca = pca.PCA(self.u, select='backbone and name CA')
        self.cumulated_variance, self.pcs = self.pca.fit()
        self.n_atoms = self.u.select_atoms('backbone and name CA').n_atoms

    def test_cum_var(self):
        assert_almost_equal(self.cumulated_variance[-1], 1)
        # makes no sense to test values here, no physical meaning

    def test_pcs(self):
        assert_equal(self.pcs.shape, (self.n_atoms*3, self.n_atoms*3))

    def test_different_steps(self):
        cum_var, pcs = self.pca.fit(step=3)

    def test_transform(self):
        self.n_components = 1
        pca_space = self.pca.transform(n_components=self.n_components)
        assert_equal(pca_space.shape,
                     (self.u.trajectory.n_frames, self.n_components))
