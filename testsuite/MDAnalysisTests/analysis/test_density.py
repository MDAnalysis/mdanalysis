# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
from __future__ import print_function

from six.moves import zip
import numpy as np
import os
import sys

from numpy.testing import (TestCase, assert_equal, assert_almost_equal, dec,
                           assert_raises)

import MDAnalysis as mda
# imported inside a skipif-protected method so that it can
# be tested in the absence of scipy
## import MDAnalysis.analysis.density

from MDAnalysisTests.datafiles import TPR, XTC
from MDAnalysisTests import module_not_found, tempdir
from mock import Mock, patch
from MDAnalysisTests.util import block_import


class TestDensity(TestCase):
    nbins = 3, 4, 5
    counts = 100
    Lmax = 10.

    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def setUp(self):
        import MDAnalysis.analysis.density

        self.bins = [np.linspace(0, self.Lmax, n+1) for n in self.nbins]
        h, edges = np.histogramdd(
            self.Lmax * np.sin(np.linspace(0, 1, self.counts * 3)).reshape(self.counts, 3),
            bins=self.bins)
        self.D = MDAnalysis.analysis.density.Density(h, edges,
                                                     parameters={'isDensity': False},
                                                     units={'length': 'A'})
        self.D.make_density()

    def test_shape(self):
        assert_equal(self.D.grid.shape, self.nbins)

    def test_edges(self):
        for dim, (edges, fixture) in enumerate(zip(
                self.D.edges, self.bins)):
            assert_almost_equal(edges, fixture,
                                err_msg="edges[{0}] mismatch".format(dim))

    def test_midpoints(self):
        midpoints = [0.5*(b[:-1] + b[1:]) for b in self.bins]
        for dim, (mp, fixture) in enumerate(zip(
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
    delta = 2.0
    selections = {'static': "name OW",
                  'dynamic': "name OW and around 4 (protein and resid 1-10)",
                  }
    references = {'static':
                      {'meandensity': 0.016764271713091212, },
                  'static_sliced':
                      {'meandensity': 0.016764270747693617, },
                  'dynamic':
                      {'meandensity': 0.00062423404854011104, },
                  }
    precision = 5

    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def setUp(self):
        self.outfile = 'density.dx'
        self.universe = mda.Universe(self.topology, self.trajectory)

    def tearDown(self):
        del self.universe

    def check_density_from_Universe(self, atomselection,
                                    ref_meandensity, **kwargs):
        import MDAnalysis.analysis.density

        with tempdir.in_tempdir():
            D = MDAnalysis.analysis.density.density_from_Universe(
                self.universe, atomselection=atomselection,
                delta=self.delta, **kwargs)
            assert_almost_equal(D.grid.mean(), ref_meandensity,
                                err_msg="mean density does not match")

            D.export(self.outfile)

            D2 = MDAnalysis.analysis.density.Density(self.outfile)
            assert_almost_equal(D.grid, D2.grid, decimal=self.precision,
                                err_msg="DX export failed: different grid sizes")


    def test_density_from_Universe(self):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['static']['meandensity'])

    def test_density_from_Universe_sliced(self):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['static_sliced']['meandensity'],
            start=1, stop=-1, step=2,
            )

    def test_density_from_Universe_update_selection(self):
        self.check_density_from_Universe(
            self.selections['dynamic'],
            self.references['dynamic']['meandensity'],
            update_selections=True)

class TestGridImport(TestCase):

    def setUp(self):
        sys.modules.pop('MDAnalysis.analysis.density', None)

    @block_import('gridData')
    def test_absence_griddata(self):
        # if gridData package is missing an ImportError should be raised
        # at the module level of MDAnalysis.analysis.density
        with assert_raises(ImportError):
            import MDAnalysis.analysis.density

    def test_presence_griddata(self):
        # no ImportError exception is raised when gridData is properly
        # imported by MDAnalysis.analysis.density

        # mock gridData in case there are testing scenarios where
        # it is not available
        mock = Mock()
        with patch.dict('sys.modules', {'gridData':mock}):
            try:
                import MDAnalysis.analysis.density
            except ImportError:
                self.fail('''MDAnalysis.analysis.density should not raise
                             an ImportError if gridData is available.''')

#issue 1231
def test_check_set_unit_keyerror():
    import MDAnalysis.analysis.density
    D = MDAnalysis.analysis.density.Density(np.zeros((2,2)), np.zeros((2,2)),
                                                    parameters={'isDensity': False},
                                                    units={'weight': 'A'})
    assert_raises(ValueError,D._check_set_unit,D.units)
    #D._check_set_unit(D.units)

def test_check_set_unit_attributeError():
    import MDAnalysis.analysis.density
    D = MDAnalysis.analysis.density.Density(np.zeros((2,2)), np.zeros((2,2)),
                                                    parameters={'isDensity': False},
                                                    units={'weight': 'A'})
    D.units = []
    assert_raises(ValueError,D._check_set_unit,D.units)