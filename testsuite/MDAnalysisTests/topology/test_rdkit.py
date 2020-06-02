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
from __future__ import absolute_import

import pytest

import MDAnalysis as mda
from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import mol2_molecule, PDB_helix

Chem = pytest.importorskip('rdkit.Chem')

class RDKitParserBase(ParserBase):
    parser = mda.topology.RDKitParser.RDKitParser
    expected_attrs = ['ids', 'names', 'elements', 'masses',
                      'resids', 'resnums',
                      'segids',
                      'bonds',
                     ]
    
    expected_n_atoms = 0
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 0

    def test_creates_universe(self, filename):
        u = mda.Universe(filename, format='RDKIT')
        assert isinstance(u, mda.Universe)

    def test_bonds_total_counts(self, top):
        assert len(top.bonds.values) == self.expected_n_bonds


class TestRDKitParserMOL2(RDKitParserBase):
    ref_filename = mol2_molecule

    expected_n_atoms = 49
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 51

    @pytest.fixture
    def filename(self):
        return Chem.MolFromMol2File(self.ref_filename, removeHs=False)

    def test_bond_orders(self, top, filename):
        # not exactly what's written in the mol2 file but the 2 double bonds in 
        # the thiophene moiety should be written as aromatic, and they should 
        # be perceived as such by RDKit
        expected = [1.5, 1, 1, 2, 1, 2, 1, 1, 1.5, 1, 1, 1.5, 2, 2, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1, 1.5, 1, 1.5, 1, 1.5, 1,
        1, 1, 1.5, 1, 1, 1, 1, 1.5, 1, 1, 2, 1, 1]
        expected = [float(i) for i in expected]
        assert sorted(top.bonds.order) == sorted(expected)


class TestRDKitParserPDB(RDKitParserBase):
    ref_filename = PDB_helix

    expected_attrs = RDKitParserBase.expected_attrs + ['resnames']
    guessed_attrs = ['types']
    
    expected_n_atoms = 137
    expected_n_residues = 13
    expected_n_segments = 1
    expected_n_bonds = 137

    @pytest.fixture
    def filename(self):
        return Chem.MolFromPDBFile(self.ref_filename, removeHs=False)