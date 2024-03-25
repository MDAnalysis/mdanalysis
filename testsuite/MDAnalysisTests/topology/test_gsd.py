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
from MDAnalysis.topology.GSDParser import HAS_GSD

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import GSD
from MDAnalysisTests.datafiles import GSD_bonds
import os


@pytest.mark.skipif(not HAS_GSD, reason='gsd not installed')
class GSDBase(ParserBase):
    parser = mda.topology.GSDParser.GSDParser
    expected_attrs = ['ids', 'names', 'resids', 'resnames', 'masses',
                      'charges', 'radii', 'types',
                      'bonds', 'angles', 'dihedrals', 'impropers']
    expected_n_bonds = 0
    expected_n_angles = 0
    expected_n_dihedrals = 0
    expected_n_impropers = 0

    def test_attr_size(self, top):
        assert len(top.ids) == top.n_atoms
        assert len(top.names) == top.n_atoms
        assert len(top.resids) == top.n_residues
        assert len(top.resnames) == top.n_residues

    def test_atoms(self, top):
        assert top.n_atoms == self.expected_n_atoms

    def test_bonds(self, top):
        assert len(top.bonds.values) == self.expected_n_bonds
        if self.expected_n_bonds:
            assert isinstance(top.bonds.values[0], tuple)
        else:
            assert top.bonds.values == []

    def test_angles(self, top):
        assert len(top.angles.values) == self.expected_n_angles
        if self.expected_n_angles:
            assert isinstance(top.angles.values[0], tuple)
        else:
            assert top.angles.values == []

    def test_dihedrals(self, top):
        assert len(top.dihedrals.values) == self.expected_n_dihedrals
        if self.expected_n_dihedrals:
            assert isinstance(top.angles.values[0], tuple)
        else:
            assert top.dihedrals.values == []

    def test_impropers(self, top):
        assert len(top.impropers.values) == self.expected_n_impropers
        if self.expected_n_impropers:
            assert isinstance(top.angles.values[0], tuple)
        else:
            assert top.impropers.values == []


@pytest.mark.skipif(not HAS_GSD, reason='gsd not installed')
class TestGSDParser(GSDBase):
    ref_filename = GSD
    expected_n_atoms = 5832
    expected_n_residues = 648
    expected_n_segments = 1


@pytest.mark.skipif(not HAS_GSD, reason='gsd not installed')
class TestGSDParserBonds(GSDBase):
    ref_filename = GSD_bonds
    expected_n_atoms = 490
    expected_n_bonds = 441
    expected_n_angles = 392
    expected_n_dihedrals = 343
    expected_n_residues = 1
    expected_n_segments = 1

    def test_bonds_identity(self, top):
        vals = top.bonds.values
        for b in ((0, 1), (1, 2), (2, 3), (3, 4)):
            assert (b in vals) or (b[::-1] in vals)
        assert ((0, 450) not in vals)

    def test_angles_identity(self, top):
        vals = top.angles.values
        for b in ((0, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 5)):
            assert (b in vals) or (b[::-1] in vals)
        assert ((0, 350, 450) not in vals)

    def test_dihedrals_identity(self, top):
        vals = top.dihedrals.values
        for b in ((0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5), (3, 4, 5, 6)):
            assert (b in vals) or (b[::-1] in vals)
        assert ((0, 250, 350, 450) not in vals)
