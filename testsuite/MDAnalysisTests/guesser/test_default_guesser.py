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
from pytest import approx
import MDAnalysis as mda

from numpy.testing import assert_equal, assert_allclose
import numpy as np
from MDAnalysis.core.topologyattrs import (Angles, Atomtypes, Atomnames)
from MDAnalysis.guesser.default_guesser import DefaultGuesser
from MDAnalysis.core.topology import Topology
from MDAnalysisTests import make_Universe
from MDAnalysisTests.core.test_fragments import make_starshape
import MDAnalysis.tests.datafiles as datafiles

from MDAnalysisTests.util import import_not_available

try:
    from rdkit import Chem
    from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
except ImportError:
    pass

requires_rdkit = pytest.mark.skipif(import_not_available("rdkit"),
                                    reason="requires RDKit")


@pytest.fixture
def default_guesser():
    return DefaultGuesser(None)


class TestGuessMasses(object):
    def test_guess_masses_from_universe(self):
        topology = Topology(3, attrs=[Atomtypes(['C', 'C', 'H'])])
        u = mda.Universe(topology)
        u.guess_TopologyAttributes(to_guess=['masses'])

        assert isinstance(u.atoms.masses, np.ndarray)
        assert_allclose(u.atoms.masses, np.array(
            [12.011, 12.011, 1.008]), rtol=1e-3, atol=0)

    def test_guess_masses_from_guesser_object(self, default_guesser):
        elements = ['H', 'Ca', 'Am']
        values = np.array([1.008, 40.08000, 243.0])
        assert_allclose(default_guesser.guess_masses(
            elements), values, rtol=1e-3, atol=0)

    def test_guess_masses_warn(self):
        topology = Topology(1, attrs=[Atomtypes(['X'])])
        with pytest.warns(UserWarning,
                          match="Failed to guess the mass for "
                                "the following atom type"):
            mda.Universe(topology, to_guess=['masses'])

    def test_guess_masses_miss(self):
        topology = Topology(2, attrs=[Atomtypes(['X', 'Z'])])
        u = mda.Universe(topology)
        assert_allclose(u.atoms.masses, np.array(
            [0.0, 0.0]), rtol=1e-3, atol=0)

    @pytest.mark.parametrize('element, value', (('H', 1.008), ('XYZ', 0.0), ))
    def test_get_atom_mass(self, element, value, default_guesser):
        default_guesser.get_atom_mass(element) == approx(value)

    def test_guess_atom_mass(self, default_guesser):
        assert default_guesser.guess_atom_mass('1H') == approx(1.008)

    def test_guess_masses_with_no_reference_elements(self):
        u = mda.Universe.empty(3)
        with pytest.raises(ValueError,
                           match=('there is no reference attributes ')):
            u.guess_TopologyAttributes('default', ['masses'])


class TestGuessTypes(object):

    def test_guess_types(self):
        topology = Topology(2, attrs=[Atomnames(['MG2+', 'C12'])])
        u = mda.Universe(topology)
        u.guess_TopologyAttributes(to_guess=['types'])
        values = DefaultGuesser(None).guess_types(u.atoms.names)
        assert isinstance(u.atoms.types, np.ndarray)
        assert_equal(u.atoms.types, np.array(['MG', 'C'], dtype=object))
        assert_equal(values, np.array(['MG', 'C'], dtype=object))

    def test_guess_atom_element(self, default_guesser):
        assert default_guesser.guess_atom_element('MG2+') == 'MG'

    def test_guess_atom_element_empty(self, default_guesser):
        assert default_guesser.guess_atom_element('') == ''

    def test_guess_atom_element_singledigit(self, default_guesser):
        assert default_guesser.guess_atom_element('1') == '1'

    def test_guess_atom_element_1H(self, default_guesser):
        assert default_guesser.guess_atom_element('1H') == 'H'
        assert default_guesser.guess_atom_element('2H') == 'H'

    @pytest.mark.parametrize('name, element', (
        ('AO5*', 'O'),
        ('F-', 'F'),
        ('HB1', 'H'),
        ('OC2', 'O'),
        ('1he2', 'H'),
        ('3hg2', 'H'),
        ('OH-', 'O'),
        ('HO', 'H'),
        ('he', 'H'),
        ('zn', 'ZN'),
        ('Ca2+', 'CA'),
        ('CA', 'C'),
    ))
    def test_guess_element_from_name(self, name, element, default_guesser):
        assert default_guesser.guess_atom_element(name) == element


def test_guess_charge(default_guesser):
    # this always returns 0.0
    assert default_guesser.guess_atom_charge('this') == approx(0.0)


def test_guess_bonds_Error():
    u = make_Universe(trajectory=True)
    with pytest.raises(ValueError):
        u.guess_TopologyAttributes(to_guess=['bonds'])


def test_guess_bond_coord_error():
    u = mda.Universe(datafiles.PDB)
    with pytest.raises(ValueError):
        DefaultGuesser(u).guess_bonds(u.atoms, [[1, 2, 2]])


def test_guess_impropers(default_guesser):
    u = make_starshape()

    ag = u.atoms[:5]
    default_guesser.guess_angles(ag.bonds)
    u.add_TopologyAttr(Angles(default_guesser.guess_angles(ag.bonds)))

    vals = default_guesser.guess_improper_dihedrals(ag.angles)
    assert_equal(len(vals), 12)


def test_guess_dihedrals_with_no_bonds():
    "Test guessing bonds and angles first before guessing "
    "dihedral for a univers with no bonds and angles "
    u = mda.Universe(datafiles.two_water_gro)
    u.guess_TopologyAttributes(to_guess=['dihedrals'])
    assert hasattr(u, 'bonds')
    assert hasattr(u, 'angles')


def bond_sort(arr):
    # sort from low to high, also within a tuple
    # e.g. ([5, 4], [0, 1], [0, 3]) -> ([0, 1], [0, 3], [4, 5])
    out = []
    for (i, j) in arr:
        if i > j:
            i, j = j, i
        out.append((i, j))
    return sorted(out)


def test_guess_bonds_water(default_guesser):
    u = mda.Universe(datafiles.two_water_gro)
    bonds = bond_sort(
        default_guesser.guess_bonds(
            u.atoms,
            u.atoms.positions,
            u.dimensions))
    assert_equal(bonds, ((0, 1),
                         (0, 2),
                         (3, 4),
                         (3, 5)))


def test_guess_bonds_adk():
    u = mda.Universe(datafiles.PSF, datafiles.DCD)
    u.guess_TopologyAttributes(force_guess=['types'])
    guesser = DefaultGuesser(None)
    bonds = bond_sort(guesser.guess_bonds(u.atoms, u.atoms.positions))
    assert_equal(np.sort(u.bonds.indices, axis=0),
                 np.sort(bonds, axis=0))


def test_guess_bonds_peptide():
    u = mda.Universe(datafiles.PSF_NAMD, datafiles.PDB_NAMD)
    u.guess_TopologyAttributes(force_guess=['types'])
    guesser = DefaultGuesser(None)
    bonds = bond_sort(guesser.guess_bonds(u.atoms, u.atoms.positions))
    assert_equal(np.sort(u.bonds.indices, axis=0),
                 np.sort(bonds, axis=0))


@pytest.mark.parametrize("smi", [
    "c1ccccc1",
    "C1=CC=CC=C1",
    "CCO",
    "c1ccccc1Cc1ccccc1",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
])
@requires_rdkit
def test_guess_aromaticities(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    expected = np.array([atom.GetIsAromatic() for atom in mol.GetAtoms()])
    u = mda.Universe(mol)
    guesser = DefaultGuesser(None)
    values = guesser.guess_aromaticities(u.atoms)
    u.guess_TopologyAttributes(to_guess=['aromaticities'])
    assert_equal(values, expected)
    assert_equal(u.atoms.aromaticities, expected)


@pytest.mark.parametrize("smi", [
    "c1ccccc1",
    "C1=CC=CC=C1",
    "CCO",
    "c1ccccc1Cc1ccccc1",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
])
@requires_rdkit
def test_guess_gasteiger_charges(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    ComputeGasteigerCharges(mol, throwOnParamFailure=True)
    expected = np.array([atom.GetDoubleProp("_GasteigerCharge")
                         for atom in mol.GetAtoms()], dtype=np.float32)
    u = mda.Universe(mol)
    guesser = DefaultGuesser(None)
    values = guesser.guess_gasteiger_charges(u.atoms)
    assert_equal(values, expected)
