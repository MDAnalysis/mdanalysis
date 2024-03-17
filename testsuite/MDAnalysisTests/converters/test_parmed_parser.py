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
from numpy.testing import assert_equal

import MDAnalysis as mda
from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    PSF_NAMD_GBIS,
    PRM
)

pmd = pytest.importorskip('parmed')


class BaseTestParmedParser(ParserBase):
    parser = mda.converters.ParmEdParser.ParmEdParser
    expected_attrs = ['ids', 'names', 'types', 'masses',
                      'charges', 'altLocs', 'occupancies',
                      'tempfactors', 'gbscreens', 'solventradii',
                      'nbindices', 'rmins', 'epsilons', 'rmin14s',
                      'epsilon14s', 'elements', 'chainIDs',
                      'resids', 'resnames', 'resnums',
                      'segids',
                      'bonds', 'ureybradleys', 'angles',
                      'dihedrals', 'impropers', 'cmaps']

    expected_n_atoms = 0
    expected_n_residues = 1
    expected_n_segments = 1
    expected_n_bonds = 0
    expected_n_angles = 0
    expected_n_dihedrals = 0
    expected_n_impropers = 0
    expected_n_cmaps = 0
    expected_n_ureybradleys = 0

    @pytest.fixture
    def filename(self):
        return pmd.load_file(self.ref_filename)

    @pytest.fixture
    def universe(self, filename):
        return mda.Universe(filename)

    def test_creates_universe(self, filename):
        u = mda.Universe(filename)
        assert isinstance(u, mda.Universe)

    def test_bonds_total_counts(self, top, filename):
        unique = set([(a.atom1.idx, a.atom2.idx)
                      for a in filename.bonds])
        assert len(top.bonds.values) == len(unique)

    def test_angles_total_counts(self, top, filename):
        unique = set([(a.atom1.idx, a.atom2.idx, a.atom3.idx)
                      for a in filename.angles])
        assert len(top.angles.values) == len(unique)

    def test_dihedrals_total_counts(self, top, filename):
        unique = set([(a.atom1.idx, a.atom2.idx, a.atom3.idx, a.atom4.idx)
                      for a in filename.dihedrals])
        assert len(top.dihedrals.values) == len(unique)

    def test_impropers_total_counts(self, top, filename):
        unique = set([(a.atom1.idx, a.atom2.idx, a.atom3.idx, a.atom4.idx)
                      for a in filename.impropers])
        assert len(top.impropers.values) == len(unique)

    def test_cmaps_total_counts(self, top, filename):
        unique = set([(a.atom1.idx, a.atom2.idx, a.atom3.idx,
                       a.atom4.idx, a.atom5.idx)
                      for a in filename.cmaps])
        assert len(top.cmaps.values) == len(unique)

    def test_ureybradleys_total_counts(self, top, filename):
        unique = set([(a.atom1.idx, a.atom2.idx)
                      for a in filename.urey_bradleys])
        assert len(top.ureybradleys.values) == len(unique)

    def test_elements(self, top):
        for erange, evals in zip(self.elems_ranges, self.expected_elems):
            assert_equal(top.elements.values[erange[0]:erange[1]], evals,
                         "unexpected element match")


class TestParmedParserPSF(BaseTestParmedParser):
    """
    PSF with CMAPs
    """

    ref_filename = PSF_NAMD_GBIS

    expected_n_atoms = 3341
    expected_n_residues = 214
    expected_n_bonds = 3365
    expected_n_angles = 6123
    expected_n_dihedrals = 8921
    expected_n_impropers = 541
    expected_n_cmaps = 212
    elems_ranges = ((100, 120),)
    # No atomic numbers set by parmed == no elements
    expected_elems = (np.array(
        ['N', 'H', 'C', 'H', 'C', 'H', 'H', 'C', 'H', 'C', 'H', 'H', 'H', 'C',
         'H', 'H', 'H', 'C', 'O', 'N',], dtype=object),)

    def test_bonds_atom_counts(self, universe):
        assert len(universe.atoms[[0]].bonds) == 4
        assert len(universe.atoms[[42]].bonds) == 1

    @pytest.mark.parametrize('value', (
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        ))
    def test_bonds_identity(self, top, value):
        vals = top.bonds.values
        assert value in vals or value[::-1] in vals

    def test_bond_types(self, universe):
        b1 = universe.bonds[0]

        assert b1.order == 1.0
        assert tuple(b1.indices) == (0, 1)
        assert b1.type.type is None

    def test_angles_atom_counts(self, universe):
        assert len(universe.atoms[[0]].angles), 9
        assert len(universe.atoms[[42]].angles), 2

    @pytest.mark.parametrize('value', (
        (1, 0, 2),
        (1, 0, 3),
        (1, 0, 4),
        ))
    def test_angles_identity(self, top, value):
        vals = top.angles.values
        assert value in vals or value[::-1] in vals

    def test_dihedrals_atom_counts(self, universe):
        assert len(universe.atoms[[0]].dihedrals) == 14

    @pytest.mark.parametrize('value', (
        (0, 4, 6, 7),
        (0, 4, 6, 8),
        (0, 4, 6, 9),
        (0, 4, 17, 18),
        ))
    def test_dihedrals_identity(self, top, value):
        vals = top.dihedrals.values
        assert value in vals or value[::-1] in vals

    @pytest.mark.parametrize('value', (
        (17, 19, 21, 41, 43),
        (60, 62, 64, 79, 81)
    ))
    def test_cmaps_identity(self, top, value):
        vals = top.cmaps.values
        assert value in vals or value[::-1] in vals


class TestParmedParserPRM(BaseTestParmedParser):
    """
    PRM
    """

    ref_filename = PRM

    expected_n_atoms = 252
    expected_n_residues = 14
    elems_ranges = ((0, 8), (30, 38))
    expected_elems = (np.array(['N', 'H', 'H', 'H', 'C', 'H', 'C', 'H'],
                               dtype=object),
                      np.array(['H', 'C', 'H', 'H', 'C', 'C', 'H', 'C'],
                               dtype=object))

    def test_bonds_atom_counts(self, universe):
        assert len(universe.atoms[[0]].bonds) == 4
        assert len(universe.atoms[[42]].bonds) == 1

    @pytest.mark.parametrize('value', (
        (10, 11),
        (10, 12),
        (4, 6),
        (4, 10),
        ))
    def test_bonds_identity(self, top, value):
        vals = top.bonds.values
        assert value in vals or value[::-1] in vals

    def test_bond_types(self, universe):
        b1 = universe.bonds[0]

        assert b1.order == 1.0
        assert tuple(b1.indices) == (0, 1)
        assert b1.type.type.k == 434
        assert b1.type.type.req == 1.010

    def test_angles_atom_counts(self, universe):
        assert len(universe.atoms[[0]].angles), 9
        assert len(universe.atoms[[42]].angles), 2

    @pytest.mark.parametrize('value', (
        (11, 10, 12),
        (10, 12, 14),
        (6, 4, 10),
        (4, 10, 11),
        ))
    def test_angles_identity(self, top, value):
        vals = top.angles.values
        assert value in vals or value[::-1] in vals

    def test_dihedrals_atom_counts(self, universe):
        assert len(universe.atoms[[0]].dihedrals) == 14

    @pytest.mark.parametrize('value', (
        (11, 10, 12, 14),
        (10, 12, 14, 16),
        ))
    def test_dihedrals_identity(self, top, value):
        vals = top.dihedrals.values
        assert value in vals or value[::-1] in vals

    def test_dihedral_types(self, universe):
        ag = universe.atoms[[10, 12, 14, 16]]
        dih = universe.dihedrals.atomgroup_intersection(ag,
                                                        strict=True)[0]
        assert len(dih.type) == 4
        for i, (phi_k, per) in enumerate((
            (2.0, 1),
            (2.0, 2),
            (0.4, 3),
            (0.0, 4),
        )):
            assert dih.type[i].type.phi_k == phi_k
            assert dih.type[i].type.per == per


def test_old_import_warning():
    wmsg = "Please import the ParmEd classes from MDAnalysis.converters"
    with pytest.warns(DeprecationWarning, match=wmsg):
        import MDAnalysis.topology.ParmEdParser
