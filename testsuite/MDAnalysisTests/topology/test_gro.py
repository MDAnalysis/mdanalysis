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
    assert_raises,
)

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    GRO,
    two_water_gro_widebox,
    GRO_empty_atom,
    GRO_missing_atomname,
)


class TestGROParser(ParserBase):
    parser = mda.topology.GROParser.GROParser
    filename = GRO
    expected_attrs = ['ids', 'names', 'resids', 'resnames', 'masses', 'types']
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 1

    def test_attr_size(self):
        for attr in ['ids', 'names']:
            assert_(len(self.top.ids) == self.top.n_atoms)
            assert_(len(self.top.names) == self.top.n_atoms)
        for attr in ['resids', 'resnames']:
            assert_(len(self.top.resids) == self.top.n_residues)
            assert_(len(self.top.resnames) == self.top.n_residues)


class TestGROWideBox(object):
    """Tests for Issue #548"""
    def test_atoms(self):
        parser = mda.topology.GROParser.GROParser
        with parser(two_water_gro_widebox) as p:
            s = p.parse()
        assert_(s.n_atoms == 6)


def test_parse_empty_atom_IOerror():
    parser = mda.topology.GROParser.GROParser
    with parser(GRO_empty_atom) as p:
      assert_raises(IOError, p.parse)


def test_parse_missing_atomname_IOerror():
    parser = mda.topology.GROParser.GROParser
    with parser(GRO_missing_atomname) as p:
      assert_raises(IOError, p.parse)
