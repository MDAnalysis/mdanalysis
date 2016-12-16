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

import numpy as np
from numpy.testing import (
    assert_,
    assert_array_equal,
    assert_equal,
    assert_raises,
)

from MDAnalysisTests.datafiles import TPR, XTC

import MDAnalysis as mda

class TestUpdatingSelection(object):
    def setUp(self):
        self.u = mda.Universe(TPR, XTC)
        self.ag = self.u.select_atoms(
            "prop x < 5 and prop y < 5 and prop z < 5")
        self.ag_updating = self.u.select_atoms(
            "prop x < 5 and prop y < 5 and prop z < 5", updating=True)
        self.ag_updating_compounded = self.u.select_atoms("around 2 group sele",
                                    sele=self.ag, updating=True)
        self.ag_updating_chained = self.u.select_atoms("around 2 group sele",
                                    sele=self.ag_updating, updating=True)
        self.ag_updating_chained2 = self.ag_updating.select_atoms("all",
                                                                updating=True)

    def test_update(self):
        assert_array_equal(self.ag_updating.indices, self.ag.indices)
        target_idxs = np.array([ 4469,  4470,  4472,  6289,  6290,  6291,
                                6292, 31313, 31314, 31315, 31316, 34661,
                                34663, 34664])
        self.u.trajectory.next()
        assert_equal(self.ag_updating._lastupdate, 0)
        assert_(not self.ag_updating.is_uptodate)
        assert_array_equal(self.ag_updating.indices, target_idxs)
        assert_(self.ag_updating.is_uptodate)
        self.ag_updating.is_uptodate = False
        assert_(self.ag_updating._lastupdate is None)

    def test_compounded_update(self):
        target_idxs0 = np.array([ 3650,  7406, 22703, 31426, 40357,
                                 40360, 41414])
        target_idxs1 = np.array([ 3650,  8146, 23469, 23472, 31426,
                                 31689, 31692, 34326, 41414])
        assert_array_equal(self.ag_updating_compounded.indices,
                           target_idxs0)
        self.u.trajectory.next()
        assert_array_equal(self.ag_updating_compounded.indices,
                           target_idxs1)

    def test_chained_update(self):
        target_idxs = np.array([ 4471,  7406, 11973, 11975, 34662, 44042])
        assert_array_equal(self.ag_updating_chained.indices,
                           self.ag_updating_compounded.indices)
        self.u.trajectory.next()
        assert_array_equal(self.ag_updating_chained.indices, target_idxs)

    def test_chained_update2(self):
        assert_array_equal(self.ag_updating_chained2.indices,
                           self.ag_updating.indices)
        self.u.trajectory.next()
        assert_array_equal(self.ag_updating_chained2.indices,
                           self.ag_updating.indices)

    def test_slice_is_static(self):
        ag_static1 = self.ag_updating[:] 
        ag_static2 = self.ag_updating.select_atoms("all") 
        assert_array_equal(ag_static1.indices, self.ag.indices)
        assert_array_equal(ag_static2.indices, self.ag.indices)
        self.u.trajectory.next()
        assert_array_equal(ag_static1.indices, self.ag.indices)
        assert_array_equal(ag_static2.indices, self.ag.indices)

    def test_kwarg_check(self):
        assert_raises(TypeError, self.u.select_atoms, "group updating",
                      {"updating":True})

class TestUpdatingSelectionNotraj(object):
    def setUp(self):
        self.u = mda.Universe(TPR)
        self.ag = self.u.select_atoms("name S*")
        self.ag_updating = self.u.select_atoms("name S*", updating=True)

    def test_update(self):
        assert_(self.ag_updating.is_uptodate)
        assert_array_equal(self.ag_updating.indices, self.ag.indices)
        assert_equal(self.ag_updating._lastupdate, -1)
        self.ag_updating.is_uptodate = False
        assert_(self.ag_updating._lastupdate is None)

