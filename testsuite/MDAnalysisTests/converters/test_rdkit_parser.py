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

import warnings
import pytest
import numpy as np
from numpy.testing import assert_equal

import MDAnalysis as mda
from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import mol2_molecule, PDB_helix, SDF_molecule

# TODO: remove these shims when RDKit
# has a release supporting NumPy 2
Chem = pytest.importorskip('rdkit.Chem')
AllChem = pytest.importorskip('rdkit.Chem.AllChem')


class RDKitParserBase(ParserBase):
    parser = mda.converters.RDKitParser.RDKitParser
    expected_attrs = ['ids', 'names', 'elements', 'masses', 'aromaticities',
                      'resids', 'resnums', 'chiralities',
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

    expected_attrs = RDKitParserBase.expected_attrs + ['charges']

    expected_n_atoms = 49
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 51

    @pytest.fixture
    def filename(self):
        return Chem.MolFromMol2File(self.ref_filename, removeHs=False)

    def _create_mol_gasteiger_charges(self):
        mol = Chem.MolFromMol2File(self.ref_filename, removeHs=False)
        AllChem.ComputeGasteigerCharges(mol)
        return mol

    def _remove_tripos_charges(self, mol):
        for atom in mol.GetAtoms():
            atom.ClearProp("_TriposPartialCharge")
    
    @pytest.fixture
    def top_gas_tripos(self):
        mol = self._create_mol_gasteiger_charges()
        return self.parser(mol).parse()

    @pytest.fixture
    def filename_gasteiger(self):
        mol = self._create_mol_gasteiger_charges()
        self._remove_tripos_charges(mol)
        return mol

    @pytest.fixture
    def top_gasteiger(self):
        mol = self._create_mol_gasteiger_charges()
        self._remove_tripos_charges(mol)
        return self.parser(mol).parse()
    

    def test_bond_orders(self, top, filename):
        expected = [bond.GetBondTypeAsDouble() for bond in filename.GetBonds()]
        assert top.bonds.order == expected
    
    def test_multiple_charge_priority(self, 
        top_gas_tripos, filename_gasteiger):
        expected = np.array([
            a.GetDoubleProp('_GasteigerCharge') for a in 
            filename_gasteiger.GetAtoms()], dtype=np.float32)
        assert_equal(expected, top_gas_tripos.charges.values)

    def test_multiple_charge_props_warning(self):
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            mol = self._create_mol_gasteiger_charges()
            # Trigger a warning.
            top = self.parser(mol).parse()
            # Verify the warning
            assert len(w) == 1
            assert "_GasteigerCharge and _TriposPartialCharge" in str(
                w[-1].message)

    def test_gasteiger_charges(self, top_gasteiger, filename_gasteiger):
        expected = np.array([
            a.GetDoubleProp('_GasteigerCharge') for a in 
            filename_gasteiger.GetAtoms()], dtype=np.float32)
        assert_equal(expected, top_gasteiger.charges.values)

    def test_tripos_charges(self, top, filename):
        expected = np.array([
            a.GetDoubleProp('_TriposPartialCharge') for a in filename.GetAtoms()
            ], dtype=np.float32)
        assert_equal(expected, top.charges.values)

    def test_aromaticity(self, top, filename):
        expected = np.array([
            atom.GetIsAromatic() for atom in filename.GetAtoms()])
        assert_equal(expected, top.aromaticities.values)


class TestRDKitParserPDB(RDKitParserBase):
    ref_filename = PDB_helix

    expected_attrs = RDKitParserBase.expected_attrs + [
        'resnames', 'altLocs', 'chainIDs', 'occupancies', 'icodes',
        'tempfactors']
    guessed_attrs = ['types']
    
    expected_n_atoms = 137
    expected_n_residues = 13
    expected_n_segments = 1
    expected_n_bonds = 137

    @pytest.fixture
    def filename(self):
        return Chem.MolFromPDBFile(self.ref_filename, removeHs=False)

    def test_partial_residueinfo_raise_error(self, filename):
        mol = Chem.RemoveHs(filename)
        mh = Chem.AddHs(mol)
        with pytest.raises(ValueError,
                           match="ResidueInfo is only partially available"):
            mda.Universe(mh)
        mh = Chem.AddHs(mol, addResidueInfo=True)
        mda.Universe(mh)
    

class TestRDKitParserSMILES(RDKitParserBase):
    ref_filename = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    guessed_attrs = ['types']

    expected_n_atoms = 24
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 25

    @pytest.fixture
    def filename(self):
        mol = Chem.MolFromSmiles(self.ref_filename)
        mol = Chem.AddHs(mol)
        return mol


class TestRDKitParserSDF(RDKitParserBase):
    ref_filename = SDF_molecule

    guessed_attrs = ['types']

    expected_n_atoms = 49
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 49

    @pytest.fixture
    def filename(self):
        return Chem.SDMolSupplier(SDF_molecule, removeHs=False)[0]

    def test_bond_orders(self, top, filename):
        expected = [bond.GetBondTypeAsDouble() for bond in filename.GetBonds()]
        assert top.bonds.order == expected
