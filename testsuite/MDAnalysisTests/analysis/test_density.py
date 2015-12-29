# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
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
import os
import itertools
import tempdir

import MDAnalysis as mda
import MDAnalysis.analysis.density

from numpy.testing import TestCase, assert_equal, assert_almost_equal



from MDAnalysisTests.datafiles import TPR, XTC

class TestDensity(TestCase):
    nbins = 3, 4, 5
    counts = 100
    Lmax = 10.

    def setUp(self):
        self.bins = [np.linspace(0, self.Lmax, n+1) for n in self.nbins]
        h, edges = np.histogramdd(self.Lmax*np.random.random((self.counts, 3)), bins=self.bins)
        self.D = MDAnalysis.analysis.density.Density(h, edges,
                                                     parameters={'isDensity': False},
                                                     units={'length': 'A'})
        self.D.make_density()

    def test_shape(self):
        assert_equal(self.D.grid.shape, self.nbins)

    def test_edges(self):
        for dim, (edges, fixture) in enumerate(itertools.izip(
                self.D.edges, self.bins)):
            assert_almost_equal(edges, fixture,
                                err_msg="edges[{0}] mismatch".format(dim))

    def test_midpoints(self):
        midpoints = [0.5*(b[:-1] + b[1:]) for b in self.bins]
        for dim, (mp, fixture) in enumerate(itertools.izip(
                self.D.midpoints, midpoints)):
            assert_almost_equal(mp, fixture,
                                err_msg="midpoints[{0}] mismatch".format(dim))

    def test_delta(self):
        deltas = np.array([self.Lmax])/np.array(self.nbins)
        assert_almost_equal(self.D.delta, deltas)

    def test_grid(self):
        dV = self.D.delta.prod()  # orthorhombic grids only!
        # counts = (rho[0] * dV[0] + rho[1] * dV[1] ...) = sum_i rho[i] * dV
        assert_almost_equal(self.D.grid.sum() * dV, self.counts)

    def test_origin(self):
        midpoints = [0.5*(b[:-1] + b[1:]) for b in self.bins]
        origin = [m[0] for m in midpoints]
        assert_almost_equal(self.D.origin, origin)



class Test_density_from_Universe(TestCase):
    topology = TPR
    trajectory = XTC
    selection = "name OW"
    delta = 2.0
    meandensity = 0.016764271713091212

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tmpdir.name , 'density.dx')

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

    def test_density_from_Universe(self):
        u = mda.Universe(self.topology, self.trajectory)
        D = MDAnalysis.analysis.density.density_from_Universe(u, atomselection=self.selection,
                                                              delta=self.delta)
        assert_almost_equal(D.grid.mean(), self.meandensity,
                            err_msg="mean density does not match")

        D.export(self.outfile)

        D2 = MDAnalysis.analysis.density.Density(self.outfile)
        assert_almost_equal(D.grid, D2.grid)


