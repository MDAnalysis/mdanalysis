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
from MDAnalysisTests.datafiles import (
    ITP_ala,  # GROMACS itp
)

class TestITP(ParserBase):
    parser = mda.topology.ITPParser.ITPParser
    ref_filename = ITP_ala
    expected_attrs = ['ids', 'names', 'types', 'masses',
                      'charges', 'chargegroups',
                      'resids', 'resnames',
                      'segids', 'moltypes',
                      'bonds', 'angles', 'dihedrals', 'impropers']
    expected_n_atoms = 63
    expected_n_residues = 10
    expected_n_segments = 1

    @pytest.fixture
    def universe(self, filename):
        return mda.Universe(filename)

    def test_bonds_total_counts(self, top):
        assert len(top.bonds.values) == 62
    
    def test_bonds_atom_counts(self, universe):
        assert len(universe.atoms[[0]].bonds) == 3
        assert len(universe.atoms[[42]].bonds) == 1

    def test_bonds_identity(self, top):
        vals = top.bonds.values
        for b in ((0, 1), (0, 2), (0, 3), (3, 4)):
            assert b in vals
        
    def test_bonds_type(self, top):
        assert top.bonds.types[0] == '2'

    def test_angles_total_counts(self, top):
        assert len(top.angles.values) == 91

    def test_angles_atom_counts(self, universe):
        assert len(universe.atoms[[0]].angles) == 5
        assert len(universe.atoms[[42]].angles) == 2

    def test_angles_identity(self, top):
        vals = top.angles.values
        for b in ((1, 0, 2), (1, 0, 3), (2, 0, 3)):
            assert (b in vals) or (b[::-1] in vals)
    
    def test_angles_type(self, top):
        assert top.angles.types[0] == '2'

    def test_dihedrals_total_counts(self, top):
        assert len(top.dihedrals.values) == 30

    def test_dihedrals_atom_counts(self, universe):
        assert len(universe.atoms[[0]].dihedrals) == 2

    def test_dihedrals_identity(self, top):
        vals = top.dihedrals.values
        for b in ((1, 0, 3, 5), (0, 3, 5, 7)):
            assert (b in vals) or (b[::-1] in vals)
    
    def test_dihedrals_identity(self, top):
        assert top.dihedrals.types[0] == '1'
    
    def test_impropers_total_counts(self, top):
        assert len(top.impropers.values) == 29

    def test_impropers_atom_counts(self, universe):
        assert len(universe.atoms[[0]].impropers) == 1

    def test_impropers_identity(self, top):
        vals = top.impropers.values
        for b in ((3, 0, 5, 4), (5, 3, 7, 6)):
            assert (b in vals) or (b[::-1] in vals)
    
    def test_impropers_identity(self, top):
        assert top.impropers.types[0] == '2'
    