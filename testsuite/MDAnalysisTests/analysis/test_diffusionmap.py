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
        # slow 6s test
        # u = MDAnalysis.Universe(PSF,DCD)
        # eg,ev=diffusionmap.diffusionmap(u)
        # assert_equal(self.eg.shape, (98,))
        # assert_almost_equal(self.eg[0], 1.0)
        # assert_almost_equal(self.eg[-1],0.03174182)
        # assert_equal(self.ev.shape, (98,98))
        # assert_almost_equal(self.ev[0,0], .095836037343022831)
        # faster
        self.u = MDAnalysis.Universe(PDB, XTC)
        self.dist = diffusionmap.DistanceMatrix(self.u, select='backbone')
        self.dist.run()
        self.epsilon = diffusionmap.EpsilonConstant(self.dist, 500)
        self.epsilon.determine_epsilon()
        self.dmap = diffusionmap.DiffusionMap(self.dist, self.epsilon)
        self.dmap.decompose_kernel()
        self.eigvals = self.dmap.eigenvalues
        self.eigvects = self.dmap.eigenvectors
        self.weights = np.ones((self.dist.nframes, ))

    def test_eg(self):
        # number of frames is trajectory is now 10 vs. 98
        assert_equal(self.eigvals.shape, (self.dist.nframes, ))
        # makes no sense to test values here, no physical meaning
        assert_almost_equal(self.eigvals[0], 1.0, decimal=5)

    def test_ev(self):
        #makes no sense to test values here, no physical meaning
        assert_equal(self.eigvects.shape, (self.dist.nframes, self.dist.nframes))

    def test_weights(self):
        dist = diffusionmap.DistanceMatrix(self.u, select='backbone')
        dist.run()
        dmap = diffusionmap.DiffusionMap(dist,self.epsilon,
                                          weights=self.weights)
        dmap.decompose_kernel()
        assert_almost_equal(self.eigvals, dmap.eigenvalues, decimal=5)
        assert_almost_equal(self.eigvects, dmap.eigenvectors, decimal=6)


    @raises(ValueError)
    def test_wrong_weights(self):
        dmap = diffusionmap.DiffusionMap(self.dist, self.epsilon,
                                          weights=np.ones((2,)))
                                          
    def test_timescaling(self):
        dmap = diffusionmap.DiffusionMap(self.dist, self.epsilon,
                                         timescale=2)
        dmap.decompose_kernel()
        assert_equal(dmap.eigenvalues.shape, (self.dist.nframes, ))
        # makes no sense to test values here, no physical meaning
        assert_almost_equal(self.eigvals[0], 1.0, decimal=5)


    def test_different_steps(self):
        dist = diffusionmap.DistanceMatrix(self.u, select='backbone',step=3)
        dist.run()
        dmap = diffusionmap.DiffusionMap(self.dist, self.epsilon)
        dmap.decompose_kernel()

    def test_transform(self):
        self.num_eigenvectors = 4
        self.dmap.transform(self.num_eigenvectors)
        assert_equal(self.dmap.diffusion_space.shape, 
                    (self.eigvects.shape[0],
                     self.num_eigenvectors))
