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
import MDAnalysis as mda
import numpy as np

from MDAnalysisTests.datafiles import waterPSF, waterDCD
from MDAnalysis.analysis.lineardensity import LinearDensity
from numpy.testing import TestCase, assert_allclose

class TestLinearDensity(TestCase):
    def setUp(self):
        self.universe = mda.Universe(waterPSF, waterDCD)
        self.sel_string = 'all'
        self.selection = self.universe.select_atoms(self.sel_string)
        
        self.xpos = np.array([ 0., 0., 0., 0.0072334, 0.00473299, 0.,
                               0., 0., 0., 0. ])

    def test_serial(self):
        ld = LinearDensity(self.selection, binsize=5)
        ld.run()
        assert_allclose(self.xpos, ld.results['x']['pos'], rtol=1e-6, atol=0)

#    def test_parallel(self):
#        ld = LinearDensity(self.universe, self.selection, binsize=5)
#        ld.run(parallel=True)
#        assert_equal(self.xpos, ld.results['x']['pos'])


