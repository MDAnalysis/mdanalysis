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
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from numpy.testing import (
    assert_,
)

import MDAnalysis as mda

from MDAnalysisTests.datafiles import (
    GRO,
    two_water_gro_widebox,
)



class TestGROParser(object):
    def setUp(self):
        parser = mda.topology.GROParser.GROParser
        with parser(GRO) as p:
            s = p.parse()
        self.top = s

    def tearDown(self):
        del self.top

    def test_attributes(self):
        for attr in ['atomids', 'atomnames', 'resids', 'resnames']:
            assert_(hasattr(self.top, attr))

    def test_attr_size(self):
        for attr in ['atomids', 'atomnames']:
            assert_(len(self.top.atomids) == self.top.n_atoms)
            assert_(len(self.top.atomnames) == self.top.n_atoms)
        for attr in ['resids', 'resnames']:
            assert_(len(self.top.resids) == self.top.n_residues)
            assert_(len(self.top.resnames) == self.top.n_residues)

    def test_size(self):
        assert_(self.top.n_atoms == 47681)
        assert_(self.top.n_residues == 11302)
        assert_(self.top.n_segments == 0)

    def test_tt_size(self):
        assert_(self.top.tt.size == (47681, 11302, 0))


class TestGROWideBox(object):
    """Tests for Issue #548"""
    def test_atoms(self):
        parser = mda.topology.GROParser.GROParser
        with parser(two_water_gro_widebox) as p:
            s = p.parse()
        assert_(s.n_atoms == 6)
