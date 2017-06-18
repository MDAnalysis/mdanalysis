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
    GRO_residwrap,
    GRO_residwrap_0base,
    GRO_sameresid_diffresname,
)


class TestGROParser(ParserBase):
    parser = mda.topology.GROParser.GROParser
    filename = GRO
    expected_attrs = ['ids', 'names', 'resids', 'resnames', 'masses']
    guessed_attrs = ['masses', 'types']
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


class TestGroResidWrapping(object):
    # resid is 5 digit field, so is limited to 100k
    # check that parser recognises when resids have wrapped
    names = ['MET', 'ARG', 'ILE', 'ILE', 'LEU', 'LEU', 'GLY']
    lengths = [19, 24, 19, 19, 19, 19, 7]
    parser = mda.topology.GROParser.GROParser


    def test_wrapping_resids(self):
        parser = mda.topology.GROParser.GROParser
        with self.parser(GRO_residwrap) as p:
            top = p.parse()

        resids = [1, 99999, 100000, 100001, 199999, 200000, 200001]

        assert_(top.tt.size == (126, 7, 1))
        for i, (r, n, l) in enumerate(zip(resids, self.names, self.lengths)):
            assert_(top.resids.values[i] == r)
            assert_(top.resnames.values[i] == n)
            assert_(len(top.tt.residues2atoms_1d([i])) == l)

    def test_wrapping_resids_0base(self):
        with self.parser(GRO_residwrap_0base) as p:
            top = p.parse()

        resids = [0, 99999, 100000, 100001, 199999, 200000, 200001]

        assert_(top.tt.size == (126, 7, 1))
        for i, (r, n, l) in enumerate(zip(resids, self.names, self.lengths)):
            assert_(top.resids.values[i] == r)
            assert_(top.resnames.values[i] == n)
            assert_(len(top.tt.residues2atoms_1d([i])) == l)

def test_sameresid_diffresname():
    parser = mda.topology.GROParser.GROParser
    with parser(GRO_sameresid_diffresname) as p:
        top = p.parse()
    resids = [9, 9]
    resnames = ['GLN', 'POPC']
    for i, (resid, resname) in enumerate(zip(resids, resnames)):
        assert_(top.resids.values[i] == resid)
        assert_(top.resnames.values[i] == resname)
