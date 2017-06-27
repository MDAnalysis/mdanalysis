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
from __future__ import print_function, absolute_import
import numpy as np
import MDAnalysis
import MDAnalysis.analysis.pca as pca

from numpy.testing import (assert_almost_equal, assert_equal, dec,
                           assert_array_almost_equal, raises)

from MDAnalysisTests.datafiles import (PSF, DCD, RANDOM_WALK, RANDOM_WALK_TOPO,
                                       waterPSF, waterDCD)
from MDAnalysisTests import module_not_found


class TestPCA(object):
    """ Test the PCA class """
    def setUp(self):
        self.u = MDAnalysis.Universe(PSF, DCD)
        self.u.transfer_to_memory()
        self.pca = pca.PCA(self.u, select='backbone and name CA',
                           align=False)
        self.pca.run()
        self.n_atoms = self.u.select_atoms('backbone and name CA').n_atoms

    def test_cov(self):
        atoms = self.u.select_atoms('backbone and name CA')
        xyz = np.zeros((self.pca.n_frames, self.pca._n_atoms*3))
        for i, ts in enumerate(self.u.trajectory):
            xyz[i] = atoms.positions.ravel()

        cov = np.cov(xyz, rowvar=0)
        assert_array_almost_equal(self.pca.cov, cov, 4)

    def test_cum_var(self):
        assert_almost_equal(self.pca.cumulated_variance[-1], 1)
        l = self.pca.cumulated_variance
        l = np.sort(l)
        assert_almost_equal(self.pca.cumulated_variance, l, 5)

    def test_pcs(self):
        assert_equal(self.pca.p_components.shape,
                     (self.n_atoms*3, self.n_atoms*3))

    def test_different_steps(self):
        dot = self.pca.transform(self.u.select_atoms('backbone and name CA'),
                                 start=5, stop=7, step=1)
        assert_equal(dot.shape, (2, self.n_atoms*3))

    def test_transform(self):
        ag = self.u.select_atoms('backbone and name CA')
        pca_space = self.pca.transform(ag, n_components=1)
        assert_equal(pca_space.shape,
                     (self.u.trajectory.n_frames, 1))

    # Accepts universe as input, but shapes are not aligned due to n_atoms
    @raises(ValueError)
    def test_transform_mismatch(self):
        pca_space = self.pca.transform(self.u, n_components=1)
        assert_equal(pca_space.shape,
                     (self.u.trajectory.n_frames, 1))

    @staticmethod
    def test_transform_universe():
        u1 = MDAnalysis.Universe(waterPSF, waterDCD)
        u2 = MDAnalysis.Universe(waterPSF, waterDCD)
        pca_test = pca.PCA(u1).run()
        pca_test.transform(u2)

    @staticmethod
    def test_cosine_content():
        rand = MDAnalysis.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)
        pca_random = pca.PCA(rand).run()
        dot = pca_random.transform(rand.atoms)
        content = pca.cosine_content(dot, 0)
        assert_almost_equal(content, .99, 1)
