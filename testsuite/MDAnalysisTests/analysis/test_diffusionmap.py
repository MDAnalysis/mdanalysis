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
from numpy.testing import (assert_almost_equal, assert_equal, raises)


from MDAnalysisTests.datafiles import PDB, XTC


class TestDiffusionmap(object):
    def __init__(self):
        self.u = MDAnalysis.Universe(PDB, XTC)
        self.dist = diffusionmap.DistanceMatrix(self.u, select='backbone')
        self.dmap = diffusionmap.DiffusionMap(self.dist)
        self.dmap.run()
        self.eigvals = self.dmap.eigenvalues
        self.eigvects = self.dmap.eigenvectors
        self.weights_ker = np.ones((self.dist.nframes, ))

    def test_eg(self):
        # number of frames is trajectory is now 10 vs. 98
        assert_equal(self.eigvals.shape, (self.dist.nframes-1, ))
        # makes no sense to test values here, no physical meaning

    def test_dist_weights(self):
        backbone = self.u.select_atoms('backbone')
        weights_atoms = np.ones(len(backbone.atoms))
        dist = diffusionmap.DistanceMatrix(self.u, select='backbone',
                                           weights=weights_atoms)
        dist.run()

    def test_kernel_weights(self):
        dist = diffusionmap.DistanceMatrix(self.u, select='backbone')
        dmap = diffusionmap.DiffusionMap(dist,
                                         manifold_density=self.weights_ker)
        dmap.run()
        assert_almost_equal(self.eigvals, dmap.eigenvalues, decimal=5)
        assert_almost_equal(self.eigvects, dmap.eigenvectors, decimal=6)

    @raises(ValueError)
    def test_wrong_kernel_weights(self):
        dmap = diffusionmap.DiffusionMap(self.dist,
                                         manifold_density=np.ones((2,)))

    @raises(ValueError)
    def test_large_matrix_exception(self):
        dist = diffusionmap.DistanceMatrix(self.u, select='backbone')
        dist.nframes = 5000
        dmap = diffusionmap.DiffusionMap(dist)

    def test_timescaling(self):
        dmap = diffusionmap.DiffusionMap(self.dist, timescale=2)
        dmap.run()
        assert_equal(dmap.eigenvalues.shape, (self.dist.nframes-1, ))

    def test_different_steps(self):
        dist = diffusionmap.DistanceMatrix(self.u, select='backbone', step=3)
        dmap = diffusionmap.DiffusionMap(dist)
        dmap.run()

    def test_transform(self):
        self.num_eigenvectors = 4
        self.dmap.transform(self.num_eigenvectors)
        assert_equal(self.dmap.diffusion_space.shape,
                     (self.eigvects.shape[0],
                      self.num_eigenvectors))

    def test_universe_init(self):
        dmap = diffusionmap.DiffusionMap(self.u)
        dmap.run()
