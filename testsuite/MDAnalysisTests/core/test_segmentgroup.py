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
from __future__ import absolute_import

from numpy.testing import (
    dec,
    assert_,
    assert_equal,
)
from unittest import skip

import MDAnalysis as mda

from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysisTests import parser_not_found


class TestSegmentGroup(object):
    # Legacy tests from before 363
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = mda.Universe(PSF, DCD)
        self.g = self.universe.atoms.segments

    def test_newSegmentGroup(self):
        """test that slicing a SegmentGroup returns a new SegmentGroup (Issue 135)"""
        g = self.universe.atoms.segments
        newg = g[:]
        assert_(isinstance(newg, mda.core.groups.SegmentGroup))
        assert_equal(len(newg), len(g))

    def test_n_atoms(self):
        assert_equal(self.g.n_atoms, 3341)

    def test_n_residues(self):
        assert_equal(self.g.n_residues, 214)

    def test_resids_dim(self):
        assert_equal(len(self.g.resids), len(self.g))
        for seg, resids in zip(self.g, self.g.resids):
            assert_(len(resids) == len(seg.residues))
            assert_equal(seg.residues.resids, resids)

    def test_resnums_dim(self):
        assert_equal(len(self.g.resnums), len(self.g))
        for seg, resnums in zip(self.g, self.g.resnums):
            assert_(len(resnums) == len(seg.residues))
            assert_equal(seg.residues.resnums, resnums)

    def test_segids_dim(self):
        assert_equal(len(self.g.segids), len(self.g))

    def test_set_segids(self):
        s = self.universe.select_atoms('all').segments
        s.segids = 'ADK'
        assert_equal(self.universe.segments.segids, ['ADK'],
                     err_msg="failed to set_segid on segments")

    def test_set_segid_updates_self(self):
        g = self.universe.select_atoms("resid 10:18").segments
        g.segids = 'ADK'
        assert_equal(g.segids, ['ADK'],
                     err_msg="old selection was not changed in place after set_segid")

    def test_atom_order(self):
        assert_equal(self.universe.segments.atoms.indices,
                     sorted(self.universe.segments.atoms.indices))


