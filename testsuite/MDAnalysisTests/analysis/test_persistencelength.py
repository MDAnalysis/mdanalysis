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

import MDAnalysis
from MDAnalysis.analysis import polymer
import numpy as np
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
        ags = [r.select_atoms('type C or type N')
               for r in self.u.residues]

        p = polymer.PersistenceLength(ags)
        return p

    def test_run(self):
        p = self._make_p()
        p.run()

        assert_(len(p.results) == 280)
        assert_almost_equal(p.lb, 1.485, 3)

    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def test_fit(self):
        p = self._make_p()
        p.run()
        p.perform_fit()

        assert_almost_equal(p.lp, 6.504, 3)
        assert_(len(p.fit) == len(p.results))


class TestFitExponential(object):
    def setUp(self):
        self.x = np.linspace(0, 250, 251)
        self.a_ref = 20.0
        self.y = np.exp(-self.x/self.a_ref)

    def tearDown(self):
        del self.x
        del self.a_ref
        del self.y

    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def test_fit_simple(self):
        a = polymer.fit_exponential_decay(self.x, self.y)
        assert_(a == self.a_ref)

    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def test_fit_noisy(self):
        y2 = self.y + (np.random.random(len(self.y)) - 0.5) * 0.05
        a = polymer.fit_exponential_decay(self.x, y2)
        assert_(np.rint(a) == self.a_ref)
