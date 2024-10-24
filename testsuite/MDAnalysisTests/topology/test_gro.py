# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import pytest

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
from numpy.testing import assert_equal, assert_allclose


class TestGROParser(ParserBase):
    parser = mda.topology.GROParser.GROParser
    ref_filename = GRO
    expected_attrs = ['ids', 'names', 'resids', 'resnames']
    guessed_attrs = ['masses', 'types']
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 1

    def test_attr_size(self, top):
        assert len(top.ids) == top.n_atoms
        assert len(top.names) == top.n_atoms
        assert len(top.resids) == top.n_residues
        assert len(top.resnames) == top.n_residues

    def test_guessed_masses(self, filename):
        u = mda.Universe(filename)
        expected = [14.007,  1.008,  1.008,  1.008, 12.011,  1.008, 12.011]
        assert_allclose(u.atoms.masses[:7], expected)

    def test_guessed_types(self, filename):
        u = mda.Universe(filename)
        expected = ['N', 'H', 'H', 'H', 'C', 'H', 'C']
        assert_equal(u.atoms.types[:7], expected)


class TestGROWideBox(object):
    """Tests for Issue #548"""
    def test_atoms(self):
        parser = mda.topology.GROParser.GROParser
        with parser(two_water_gro_widebox) as p:
            s = p.parse()
        assert s.n_atoms == 6


def test_parse_empty_atom_IOerror():
    parser = mda.topology.GROParser.GROParser
    with parser(GRO_empty_atom) as p:
        with pytest.raises(IOError):
            p.parse()


def test_parse_missing_atomname_IOerror():
    parser = mda.topology.GROParser.GROParser
    with parser(GRO_missing_atomname) as p:
        with pytest.raises(IOError):
            p.parse()


class TestGroResidWrapping(object):
    # resid is 5 digit field, so is limited to 100k
    # check that parser recognises when resids have wrapped
    names = ['MET', 'ARG', 'ILE', 'ILE', 'LEU', 'LEU', 'GLY']
    lengths = [19, 24, 19, 19, 19, 19, 7]
    parser = mda.topology.GROParser.GROParser

    @pytest.mark.parametrize('parser, resids', (
        (GRO_residwrap, [1, 99999, 100000, 100001, 199999, 200000, 200001]),
        (GRO_residwrap_0base, [0, 99999, 100000, 100001, 199999, 200000,
                               200001])

    ))
    def test_wrapping_resids(self, parser, resids):
        with self.parser(parser) as p:
            top = p.parse()

        assert top.tt.size == (126, 7, 1)
        assert_equal(top.resids.values, resids)
        assert_equal(top.resids.values, resids)
        assert_equal(top.resnames.values, self.names)
        for i, l in enumerate(self.lengths):
            assert len(top.tt.residues2atoms_1d([i])) == l


def test_sameresid_diffresname():
    parser = mda.topology.GROParser.GROParser
    with parser(GRO_sameresid_diffresname) as p:
        top = p.parse()
    resids = [9, 9]
    resnames = ['GLN', 'POPC']
    for i, (resid, resname) in enumerate(zip(resids, resnames)):
        assert top.resids.values[i] == resid
        assert top.resnames.values[i] == resname
