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
import numpy as np
from simtk.openmm import app

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import PDB, PDB_conect

class OpenMMTopologyBase(ParserBase):
    parser = mda.topology.OpenMMParser.OpenMMTopologyParser
    expected_attrs = ['ids', 'names', 'resids', 'resnames', 'masses',
                      'bonds']
    expected_n_bonds = 0

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

    def test_resids(self, top):
        assert len(top.resids.values) == self.expected_n_residues
        if self.expected_n_residues:
            assert isinstance(top.resids.values, np.ndarray)
        else:
            assert top.resids.values == []

    def test_resnames(self, top):
        assert len(top.resnames.values) == self.expected_n_residues
        if self.expected_n_residues:
            assert isinstance(top.resnames.values, np.ndarray)
        else:
            assert top.resnames.values == []

    def test_resnums(self, top):
        assert len(top.resnums.values) == self.expected_n_residues
        if self.expected_n_residues:
            assert isinstance(top.resnums.values, np.ndarray)
        else:
            assert top.resnums.values == []

    def test_segids(self, top):
        assert len(top.segids.values) == self.expected_n_segments
        if self.expected_n_segments:
            assert isinstance(top.segids.values, np.ndarray)
        else:
            assert top.segids.values == []


class TestOpenMMTopologyParser(OpenMMTopologyBase):
    ref_filename = app.PDBFile(PDB).topology
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 1

class TestOpenMMTopologyParserBonds(OpenMMTopologyBase):
    ref_filename = app.PDBFile(PDB_conect).topology
    expected_n_atoms = 10
    expected_n_bonds = 8
    expected_n_residues = 1
    expected_n_segments = 1


