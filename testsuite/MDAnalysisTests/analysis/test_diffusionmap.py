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
        self.dmap = diffusionmap.DiffusionMap(self.u, select='backbone', k=5)
        self.dmap.run()
        self.eigvals = self.dmap.eigenvalues
        self.eigvects = self.dmap.eigenvectors
        self.weights = np.ones((self.dmap.nframes, ))

    def test_eg(self):
        # number of frames is trajectory is now 10 vs. 98
        assert_equal(self.eigvals.shape, (self.dmap.nframes, ))
        assert_almost_equal(self.eigvals[0], 1.0, decimal=5)
        assert_almost_equal(self.eigvals[-1], 0.0142, decimal=3)

    def test_ev(self):
        assert_equal(self.eigvects.shape, (self.dmap.nframes, self.dmap.nframes))
        assert_almost_equal(self.eigvects[0, 0], -0.3019, decimal=2)

    def test_weights(self):
        dmap2 = diffusionmap.DiffusionMap(self.u, select='backbone',
                                          weights=self.weights, k=5)
        dmap2.run()
        assert_almost_equal(self.eigvals, dmap2.eigenvalues, decimal=5)
        assert_almost_equal(self.eigvects, dmap2.eigenvectors, decimal=6)

    @raises(ValueError)
    def test_wrong_weights(self):
        dmap2 = diffusionmap.DiffusionMap(self.u, select='backbone',
                                          weights=np.ones((2,)),
                                          k=5)

    def test_constant_epsilon(self):
        dmap3 = diffusionmap.DiffusionMap(self.u, select='backbone', k=5,
                                          epsilon=4.62827648)
        dmap3.run()
        assert_almost_equal(self.eigvals, dmap3.eigenvalues, decimal=5)
