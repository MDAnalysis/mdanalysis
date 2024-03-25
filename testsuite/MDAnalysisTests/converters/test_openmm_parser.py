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


try:
    from openmm.app import Element, Topology
    from openmm.unit import daltons
    from openmm import app
except ImportError:
    try:
        from simtk.openmm.app import Element, Topology
        from simtk.unit import daltons
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
        "types"
    ]
    expected_n_bonds = 0

    @pytest.fixture()
    def top(self, filename):
        with self.parser(filename) as p:
            yield p.parse()

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
        else:
            assert not hasattr(top, 'elements')

    def test_atomtypes(self, top):
        assert len(top.types.values) == self.expected_n_atoms
        if self.expected_n_atoms:
            assert isinstance(top.types.values, np.ndarray)
        else:
            assert top.types.values == []

    def test_masses(self, top):
        assert len(top.masses.values) == self.expected_n_atoms
        if self.expected_n_atoms:
            assert isinstance(top.masses.values, np.ndarray)
            assert all(isinstance(mass, np.float64)
                       for mass in top.masses.values)
        else:
            assert top.masses.values == []

    def test_guessed_attributes(self, filename):
        u = mda.Universe(filename, topology_format="OPENMMTOPOLOGY")
        u_guessed_attrs = [attr.attrname for attr
                           in u._topology.guessed_attributes]
        for attr in self.guessed_attrs:
            assert hasattr(u.atoms, attr)
            assert attr in u_guessed_attrs


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
        "types"
    ]
    expected_n_bonds = 0

    @pytest.fixture()
    def top(self, filename):
        with self.parser(filename) as p:
            yield p.parse()

    def test_creates_universe(self, filename):
        """Check that Universe works with this Parser"""
        u = mda.Universe(filename, topology_format="OPENMMAPP")
        assert isinstance(u, mda.Universe)

    def test_guessed_attributes(self, filename):
        u = mda.Universe(filename, topology_format="OPENMMAPP")
        for attr in self.guessed_attrs:
            assert hasattr(u.atoms, attr)


class TestOpenMMTopologyParser(OpenMMTopologyBase):
    ref_filename = app.PDBFile(CONECT).topology
    expected_n_atoms = 1890
    expected_n_residues = 199
    expected_n_segments = 3
    expected_n_bonds = 1922


class TestOpenMMTopologyParserWithPartialElements(OpenMMTopologyBase):
    parser = mda.converters.OpenMMParser.OpenMMTopologyParser
    ref_filename = app.PDBFile(PDB).topology
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 1
    expected_n_bonds = 25533

    def test_with_partial_elements(self):
        wmsg1 = ("Element information missing for some atoms. "
                 "These have been given an empty element record ")

        wmsg2 = (
            "For absent elements, atomtype has been  "
            "set to 'X' and mass has been set to 0.0. "
            "If needed these can be guessed using "
            "universe.guess_TopologyAttrs(to_guess=['masses', 'types']). "
            "(for MDAnalysis version 2.x this is done automatically,"
            " but it will be removed in future versions).")

        with pytest.warns(UserWarning) as warnings:
            mda_top = self.parser(self.ref_filename).parse()
            assert mda_top.types.values[3344] == 'X'
            assert mda_top.types.values[3388] == 'X'
            assert mda_top.elements.values[3344] == ''
            assert mda_top.elements.values[3388] == ''
            assert mda_top.masses.values[3344] == 0.0
            assert mda_top.masses.values[3388] == 0.0

            assert len(warnings) == 2
            assert str(warnings[0].message) == wmsg1
            assert str(warnings[1].message) == wmsg2


def test_no_elements_warn():
    parser = mda.converters.OpenMMParser.OpenMMTopologyParser
    omm_top = app.PDBFile(CONECT).topology
    for a in omm_top.atoms():
        a.element = None

    wmsg = (
        "Element information is missing for all the atoms. "
        "Elements attribute will not be populated. "
        "Atomtype attribute will be guessed using atom "
        "name and mass will be guessed using atomtype."
        "for MDAnalysis version 2.x this is done automatically, "
        "but it will be removed in future versions. "
        "These can be guessed using "
        "universe.guess_TopologyAttrs(to_guess=['masses', 'types']) "
        "See MDAnalysis.guessers.")

    with pytest.warns(UserWarning) as warnings:
        mda_top = parser(omm_top).parse()
        assert str(warnings[0].message) == wmsg


def test_invalid_element_symbols():
    parser = mda.converters.OpenMMParser.OpenMMTopologyParser
    omm_top = Topology()
    bad1 = Element(0, "*", "*", 0*daltons)
    bad2 = Element(0, "?", "?", 6*daltons)
    bad3 = None
    silver = Element.getBySymbol('Ag')
    chain = omm_top.addChain()
    res = omm_top.addResidue('R', chain)
    omm_top.addAtom(name='Ag', element=silver, residue=res)
    omm_top.addAtom(name='Bad1', element=bad1, residue=res)
    omm_top.addAtom(name='Bad2', element=bad2, residue=res)
    omm_top.addAtom(name='Bad3', element=bad3, residue=res)
    mda_top = parser(omm_top).parse()

    assert mda_top.types.values[0] == 'Ag'
    assert mda_top.types.values[1] == '*'
    assert mda_top.types.values[2] == '?'
    assert mda_top.types.values[3] == 'X'
    assert mda_top.elements.values[0] == 'Ag'
    assert mda_top.elements.values[1] == ''
    assert mda_top.elements.values[2] == ''
    assert mda_top.elements.values[3] == ''
    assert mda_top.masses.values[0] == 107.86822
    assert mda_top.masses.values[1] == 0.0
    assert mda_top.masses.values[2] == 6.0
    assert mda_top.masses.values[3] == 0.0


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
