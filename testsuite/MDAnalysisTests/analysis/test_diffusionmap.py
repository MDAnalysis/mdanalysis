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
import MDAnalysis.analysis.diffusionmap as diffusionmap
from numpy.testing import (assert_almost_equal, assert_equal,
                           assert_array_almost_equal,raises)


from MDAnalysisTests.datafiles import PDB, XTC


class TestDiffusionmap(object):
    def setUp(self):
        self.u = MDAnalysis.Universe(PDB, XTC)
        self.dist = diffusionmap.DistanceMatrix(self.u, select='backbone')
        self.dmap = diffusionmap.DiffusionMap(self.dist)
        self.dmap.run()
        self.eigvals = self.dmap.eigenvalues
        self.eigvects = self.dmap._eigenvectors

    def test_eg(self):
        # number of frames is trajectory is now 10 vs. 98
        assert_equal(self.eigvals.shape, (self.dist.n_frames, ))
        # makes no sense to test values here, no physical meaning

    def test_dist_weights(self):
        backbone = self.u.select_atoms('backbone')
        weights_atoms = np.ones(len(backbone.atoms))
        self.dist = diffusionmap.DistanceMatrix(self.u, select='backbone',
                                                weights=weights_atoms,
                                                step=3)
        self.dist.run()
        self.dmap = diffusionmap.DiffusionMap(self.dist)
        self.dmap.run()
        assert_array_almost_equal(self.dmap.eigenvalues, [1,1,1,1], 4)
        assert_array_almost_equal(self.dmap._eigenvectors,
                                  ([[ 0, 0, 1, 0],
                                    [ 0, 0, 0, 1],
                                    [ -.707,-.707, 0, 0],
                                    [  .707,-.707, 0, 0]]) ,2)

    def test_different_steps(self):
        self.dmap = diffusionmap.DiffusionMap(self.u, select='backbone', step=3)
        self.dmap.run()
        assert_equal(self.dmap._eigenvectors.shape, (4,4))

    def test_transform(self):
        self.n_eigenvectors = 4
        self.dmap = diffusionmap.DiffusionMap(self.u)
        self.dmap.run()
        diffusion_space = self.dmap.transform(self.n_eigenvectors,1)
        assert_equal(diffusion_space.shape,
                     (self.eigvects.shape[0],
                      self.n_eigenvectors))
