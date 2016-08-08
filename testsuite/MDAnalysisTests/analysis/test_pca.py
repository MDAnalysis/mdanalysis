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
from MDAnalysis.analysis.align import _fit_to
from scipy.integrate import simps

from numpy.testing import (assert_almost_equal, assert_equal,
                           assert_array_almost_equal, raises)


from MDAnalysisTests.datafiles import PDB, XTC, RANDOM_WALK, waterPSF, waterDCD

class TestPCA(object):
    def setUp(self):
        self.u = MDAnalysis.Universe(PDB, XTC)
        self.pca = pca.PCA(self.u.atoms, select='backbone and name CA',
                           align=False)
        self.pca.run()
        self.n_atoms = self.u.select_atoms('backbone and name CA').n_atoms

    def test_cov(self):
        atoms = self.u.select_atoms('backbone and name CA')
        xyz = np.zeros((10, 642))
        for i, ts in enumerate(self.u.trajectory):
            xyz[i] = atoms.positions.ravel()

        cov = np.cov(xyz, rowvar=0)
        assert_array_almost_equal(self.pca.cov, cov, 4)

    def test_cum_var(self):
        assert_almost_equal(self.pca.cumulated_variance[-1], 1)

    def test_pcs(self):
        assert_equal(self.pca.p_components.shape,
                     (self.n_atoms*3, self.n_atoms*3))

    def test_different_steps(self):
        dot = self.pca.transform(self.u.select_atoms('backbone and name CA'),
                                 start=5, stop=7, step=1)
        assert_equal(dot.shape, (2, self.n_atoms*3))

    def test_transform(self):
        self.ag = self.u.select_atoms('backbone and name CA')
        self.pca_space = self.pca.transform(self.ag, n_components=1)
        assert_equal(self.pca_space.shape,
                     (self.u.trajectory.n_frames, 1))

    # Accepts universe as input, but shapes are not aligned due to n_atoms
    @raises(ValueError)
    def test_transform_mismatch(self):
        self.ag = self.u.select_atoms('backbone and name CA')
        self.pca_space = self.pca.transform(self.u, n_components=1)
        assert_equal(self.pca_space.shape,
                     (self.u.trajectory.n_frames, 1))

    def test_transform_universe(self):
        u1 = MDAnalysis.Universe(waterPSF, waterDCD)
        u2 = MDAnalysis.Universe(waterPSF, waterDCD)
        pca_test = pca.PCA(u1).run()
        pca_test.transform(u2)

    def test_cosine_content(self):
        rand = MDAnalysis.Universe(RANDOM_WALK)
        pca_random = pca.PCA(rand.atoms).run()
        dot = pca_random.transform(rand.atoms)
        content = cosine_content(dot, 0)
        assert_almost_equal(content, .99, 1)


def cosine_content(pc, i):
    t = np.arange(len(pc))
    T = len(pc)
    cos = np.cos(np.pi * t * (i + 1) / T)
    return (2.0 / T) * (simps(cos*pc[:, i])) ** 2 / simps(pc[:, i] ** 2)
