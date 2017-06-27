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
from __future__ import print_function, division, absolute_import

import MDAnalysis
from MDAnalysis.analysis import polymer
from MDAnalysis.exceptions import NoDataError

import numpy as np
import matplotlib

from numpy.testing import (
    assert_,
    assert_almost_equal,
    assert_raises,
    dec
)

from MDAnalysisTests.datafiles import Plength
from MDAnalysisTests import module_not_found


class TestPersistenceLength(object):
    def setUp(self):
        self.u = MDAnalysis.Universe(Plength)

    def tearDown(self):
        del self.u

    def test_ag_VE(self):
        ags = [self.u.atoms[:10], self.u.atoms[10:110]]
        assert_raises(ValueError, polymer.PersistenceLength, ags)

    def _make_p(self):
        ags = [r.atoms.select_atoms('name C* N*')
               for r in self.u.residues]

        p = polymer.PersistenceLength(ags)
        return p

    def test_run(self):
        p = self._make_p()
        p.run()

        assert_(len(p.results) == 280)
        assert_almost_equal(p.lb, 1.485, 3)

    def test_fit(self):
        p = self._make_p()
        p.run()
        p.perform_fit()

        assert_almost_equal(p.lp, 6.504, 3)
        assert_(len(p.fit) == len(p.results))

    def test_plot_ax_return(self):
        '''Ensure that a matplotlib axis object is
        returned when plot() is called.'''
        p = self._make_p()
        p.run()
        p.perform_fit()
        actual = p.plot()
        expected = matplotlib.axes.Axes
        assert_(isinstance(actual, expected))

    def test_raise_NoDataError(self):
        '''Ensure that a NoDataError is raised if
        perform_fit() is called before the run()
        method of AnalysisBase.'''
        p = self._make_p()
        assert_raises(NoDataError, p.perform_fit)

class TestFitExponential(object):
    def setUp(self):
        self.x = np.linspace(0, 250, 251)
        self.a_ref = 20.0
        self.y = np.exp(-self.x / self.a_ref)

    def tearDown(self):
        del self.x
        del self.a_ref
        del self.y

    def test_fit_simple(self):
        a = polymer.fit_exponential_decay(self.x, self.y)
        assert_(a == self.a_ref)

    def test_fit_noisy(self):
        noise = np.sin(self.x) * 0.01
        y2 = noise + self.y

        a = polymer.fit_exponential_decay(self.x, y2)

        assert_almost_equal(a, self.a_ref, decimal=3)
        #assert_(np.rint(a) == self.a_ref)
