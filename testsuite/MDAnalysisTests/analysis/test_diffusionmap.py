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
import MDAnalysis
import MDAnalysis.analysis.diffusionmap as diffusionmap
from numpy.testing import (assert_almost_equal, assert_equal)


from MDAnalysisTests.datafiles import PDB, XTC


class TestDiffusionmap(object):
    def __init__(self):
        #slow 6s test
        #u = MDAnalysis.Universe(PSF,DCD)
        #eg,ev=diffusionmap.diffusionmap(u)
        #assert_equal(self.eg.shape, (98,))
        #assert_almost_equal(self.eg[0], 1.0)
        #assert_almost_equal(self.eg[-1],0.03174182)
        #assert_equal(self.ev.shape, (98,98))
        #assert_almost_equal(self.ev[0,0], .095836037343022831)
        #faster
        u = MDAnalysis.Universe(PDB, XTC)
        eg, ev = diffusionmap.diffusionmap(u, select = 'backbone', k = 5)
        self.eg = eg
        self.ev = ev
       

    def test_eg(self):
        assert_equal(self.eg.shape, (10, ))
        assert_almost_equal(self.eg[0], 1.0, decimal=5)
        assert_almost_equal(self.eg[-1], 0.0142, decimal = 3)


    def test_ev(self):
        assert_equal(self.ev.shape, (10, 10))
        assert_almost_equal(self.ev[0, 0], -0.3019, decimal = 2)



