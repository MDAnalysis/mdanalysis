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

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import CONECT, PDBX, PDB
from copy import deepcopy

try:
    from openmm import app
except ImportError:
    try:
        from simtk.openmm import app
    except ImportError:
        pytest.skip(allow_module_level=True)


class OpenMMTopologyBase(ParserBase):
    parser = mda.converters.OpenMMParser.OpenMMTopologyParser
    expected_attrs = [
        "ids",
        "names",
        "resids",
        "resnames",
        "masses",
        "bonds",
        "chainIDs",
        "elements",
    ]
    expected_n_bonds = 0

    def test_creates_universe(self, filename):
        """Check that Universe works with this Parser"""
        u = mda.Universe(filename, topology_format="OPENMMTOPOLOGY")
        assert isinstance(u, mda.Universe)

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
        assert all(isinstance(segid, str) for segid in top.segids.values)
        if self.expected_n_segments:
            assert isinstance(top.segids.values, np.ndarray)
        else:
            assert top.segids.values == []

    def test_elements(self, top):
        if 'elements' in self.expected_attrs:
            assert len(top.elements.values) == self.expected_n_atoms
            assert isinstance(top.elements.values, np.ndarray)
            assert all(isinstance(elem, str) for elem in top.elements.values)

    def test_atomtypes(self, top):
        if 'atomtypes' in self.expected_attrs:
            assert len(top.atomtypes.values) == self.expected_n_atoms
            assert isinstance(top, np.ndarray)
            assert all(isinstance(atomtype, str)
                       for atomtype in top.atomtypes.values)

    def test_masses(self, top):
        assert len(top.masses.values) == self.expected_n_atoms
        assert isinstance(top.masses.values, np.ndarray)
        assert all(isinstance(mass, np.float64)
                   for mass in top.masses.values)


class OpenMMAppTopologyBase(OpenMMTopologyBase):
    parser = mda.converters.OpenMMParser.OpenMMAppTopologyParser
    expected_attrs = [
        "ids",
        "names",
        "resids",
        "resnames",
        "masses",
        "bonds",
        "chainIDs",
        "elements",
    ]
    expected_n_bonds = 0

    def test_creates_universe(self, filename):
        """Check that Universe works with this Parser"""
        u = mda.Universe(filename, topology_format="OPENMMAPP")
        assert isinstance(u, mda.Universe)


class TestOpenMMTopologyParser(OpenMMTopologyBase):
    ref_filename = app.PDBFile(CONECT).topology
    expected_n_atoms = 1890
    expected_n_residues = 199
    expected_n_segments = 3
    expected_n_bonds = 1922


class TestOpenMMTopologyParserWithNoElements(OpenMMTopologyBase):
    ref_filename = app.PDBFile(PDB).topology
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 1
    expected_n_bonds = 25533

    def test_with_partial_elements(self, top):
        if 'elements' in self.expected_attrs:
            omm_topology = self.ref_filename
            wmsg = ("Unknown element None found for some atoms. "
                    "These have been given an empty element record "
                    "with their atomtype set to 'X' "
                    "and their mass set to 0.0. "
                    "If needed they can be guessed using "
                    "MDAnalysis.topology.guessers.")
            with pytest.warns(UserWarning, match=wmsg):
                mda_top = self.parser(PDB) \
                          ._mda_topology_from_omm_topology(omm_topology)

    def test_with_no_elements(self, top):
        self.expected_attrs.remove('elements')
        omm_topology = deepcopy(self.ref_filename)
        for a in omm_topology.atoms():
            a.element = None

        wmsg = ("Element information is missing, elements attribute "
                "will not be populated. "
                "Atomtype attribute will be guessed using atom "
                "name and mass will be guessed using atomtype."
                "See MDAnalysis.topology.guessers."
                )
        with pytest.warns(UserWarning, match=wmsg):
            mda_top = self.parser(PDB) \
                      ._mda_topology_from_omm_topology(omm_topology)
        self.expected_attrs.append('elements')


class TestOpenMMPDBFileParser(OpenMMAppTopologyBase):
    ref_filename = app.PDBFile(CONECT)
    expected_n_atoms = 1890
    expected_n_residues = 199
    expected_n_segments = 3
    expected_n_bonds = 1922


class TestOpenMMPDBxFileParser(OpenMMAppTopologyBase):
    ref_filename = app.PDBxFile(PDBX)
    expected_n_atoms = 60
    expected_n_residues = 7
    expected_n_segments = 1
    expected_n_bonds = 62
