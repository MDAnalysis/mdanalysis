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
from numpy.testing import assert_equal, assert_allclose
import numpy as np
from io import StringIO

import MDAnalysis as mda
from MDAnalysisTests.topology.base import ParserBase
from MDAnalysis.tests.datafiles import (
    LAMMPSdata,
    LAMMPScnt,
    LAMMPScnt2,
    LAMMPShyd,
    LAMMPShyd2,
    LAMMPSdata_deletedatoms,
    LAMMPSDUMP,
    LAMMPSDUMP_long,
    LAMMPSdata_PairIJ,
)


class LammpsBase(ParserBase):
    parser = mda.topology.LAMMPSParser.DATAParser
    expected_n_segments = 1
    expected_attrs = ['types', 'resids', 'masses', 'charges',
                      'bonds', 'angles', 'dihedrals', 'impropers']

    def test_n_atom_types(self, top):
        assert_equal(len(set(top.types.values)), self.expected_n_atom_types)

    def test_n_bonds(self, top):
        if self.ref_n_bonds:
            assert_equal(len(top.bonds.values), self.ref_n_bonds)
        else:
            assert not hasattr(top, 'bonds')

    def test_bond_member(self, top):
        if self.ref_n_bonds:
            assert self.ref_bond in top.bonds.values

    def test_n_angles(self, top):
        if self.ref_n_angles:
            assert_equal(len(top.angles.values), self.ref_n_angles)
        else:
            assert not hasattr(self.top, 'angles')

    def test_angle_member(self, top):
        if self.ref_n_angles:
            assert self.ref_angle in top.angles.values

    def test_n_dihedrals(self, top):
        if self.ref_n_dihedrals:
            assert_equal(len(top.dihedrals.values), self.ref_n_dihedrals)
        else:
            assert not hasattr(self.top, 'dihedrals')

    def test_dihedral_member(self, top):
        if self.ref_n_dihedrals:
            assert self.ref_dihedral in top.dihedrals.values

    def test_n_impropers(self, top):
        if self.ref_n_impropers:
            assert_equal(len(top.impropers.values), self.ref_n_impropers)
        else:
            assert not hasattr(self.top, 'impropers')

    def test_improper_member(self, top):
        if self.ref_n_impropers:
            assert self.ref_improper in top.impropers.values

    def test_creates_universe(self, filename):
        u = mda.Universe(filename, format='DATA')

    def test_guessed_attributes(self, filename):
        u = mda.Universe(filename, format='DATA')
        for attr in self.guessed_attrs:
            assert hasattr(u.atoms, attr)


class TestLammpsData(LammpsBase):
    """Tests the reading of lammps .data topology files.

    The reading of coords and velocities is done separately in
    test_coordinates
    """
    ref_filename = LAMMPSdata
    expected_n_atoms = 18364
    expected_n_atom_types = 10
    expected_n_residues = 25
    ref_n_bonds = 18336
    ref_bond = (12, 14)
    ref_n_angles = 29904
    ref_angle = (3, 6, 9)
    ref_n_dihedrals = 5712
    ref_dihedral = (82, 85, 88, 89)
    ref_n_impropers = 0


class TestLAMMPSCNT(LammpsBase):
    expected_n_atoms = 604
    expected_n_atom_types = 1
    expected_n_residues = 1
    ref_n_bonds = 906
    ref_bond = (9, 467)
    ref_n_angles = 1812
    ref_angle = (17, 16, 31)
    ref_n_dihedrals = 3624
    ref_dihedral = (22, 39, 40, 41)
    ref_n_impropers = 604
    ref_improper = (210, 159, 212, 566)

    @pytest.fixture(params=[LAMMPScnt, LAMMPScnt2])
    def filename(self, request):
        return request.param


class TestLAMMPSHYD(LammpsBase):
    expected_n_atoms = 2
    expected_n_atom_types = 1
    expected_n_residues = 1
    ref_n_bonds = 1
    ref_bond = (0, 1)
    ref_n_angles = 0
    ref_n_dihedrals = 0
    ref_n_impropers = 0

    @pytest.fixture(params=[LAMMPShyd, LAMMPShyd2])
    def filename(self, request):
        return request.param


class TestLAMMPSDeletedAtoms(LammpsBase):
    ref_filename = LAMMPSdata_deletedatoms
    expected_n_atoms = 10
    expected_n_atom_types = 2
    expected_n_residues = 1

    ref_n_bonds = 9
    ref_bond = (0, 3)
    ref_n_angles = 0
    ref_n_dihedrals = 0
    ref_n_impropers = 0

    def test_atom_ids(self, filename):
        u = mda.Universe(filename)

        assert_equal(u.atoms.ids,
                     [1, 10, 1002, 2003, 2004, 2005, 2006, 2007, 2008, 2009])

    def test_traj(self, filename):
        u = mda.Universe(filename)

        assert_equal(u.atoms.positions,
                     np.array([[11.8998565674, 48.4455718994, 19.0971984863],
                               [14.5285415649, 50.6892776489, 19.9419136047],
                               [12.8466796875, 48.1473007202, 18.6461906433],
                               [11.0093536377, 48.7145767212, 18.5247917175],
                               [12.4033203125, 49.2582168579, 20.2825050354],
                               [13.0947723389, 48.8437194824, 21.0175533295],
                               [11.540184021,  49.6138534546, 20.8459072113],
                               [13.0085144043, 50.6062469482, 19.9141769409],
                               [12.9834518433, 51.1562423706, 18.9713554382],
                               [12.6588821411, 51.4160842896, 20.5548400879]],
                              dtype=np.float32))


class TestLammpsDataPairIJ(LammpsBase):
    """Tests the reading of lammps .data topology file with a
    PairIJ Coeffs section
    """

    expected_attrs = ['types', 'resids', 'masses',
                      'bonds', 'angles', 'dihedrals', 'impropers']
    ref_filename = LAMMPSdata_PairIJ
    expected_n_atoms = 800
    expected_n_atom_types = 2
    expected_n_residues = 1
    ref_n_bonds = 799
    ref_bond = (397, 398)
    ref_n_angles = 390
    ref_angle = (722, 723, 724)
    ref_n_dihedrals = 385
    ref_dihedral = (722, 723, 724, 725)
    ref_n_impropers = 0


LAMMPS_NORESID = """\
LAMMPS data file via write_data, version 11 Aug 2017, timestep = 0

1 atoms
1 atom types

0.0000000000000000e+00 4.3008000000000003e+01 xlo xhi
0.0000000000000000e+00 4.3008000000000003e+01 ylo yhi
0.0000000000000000e+00 4.3008000000000003e+01 zlo zhi

Masses

1 28

Pair Coeffs # lj/cut

1 0.1892 3.75

Atoms # atomic

1 1 3.7151744275286681e+01 1.8684434743140471e+01 1.9285127961842125e+01 0 0 0
"""

def test_noresid():
    u = mda.Universe(StringIO(LAMMPS_NORESID), format='data',
                     atom_style='id type x y z')
    assert len(u.atoms) == 1

    assert_equal(u.atoms[0].mass, 28.0)
    assert_equal(u.atoms.positions,
                 np.array([[3.7151744275286681e+01,
                            1.8684434743140471e+01,
                            1.9285127961842125e+01]], dtype=np.float32))

def test_noresid_failure():
    with pytest.raises(
            ValueError,
            match='.+?You can supply a description of the atom_style.+?',
    ):
        u = mda.Universe(StringIO(LAMMPS_NORESID), format='data')


def test_interpret_atom_style():
    style = mda.topology.LAMMPSParser.DATAParser._interpret_atom_style(
        'id charge type z y x')

    assert isinstance(style, dict)
    assert style['id'] == 0
    assert style['type'] == 2
    assert style['charge'] == 1
    assert style['x'] == 5
    assert style['y'] == 4
    assert style['z'] == 3


def test_interpret_atom_style_missing():
    with pytest.raises(ValueError,
                       match='atom_style string missing required.+?'):
        style = mda.topology.LAMMPSParser.DATAParser._interpret_atom_style(
            'id charge z y x')


class TestDumpParser(ParserBase):
    expected_attrs = ['types', 'masses']
    expected_n_atoms = 24
    expected_n_residues = 1
    expected_n_segments = 1

    parser = mda.topology.LAMMPSParser.LammpsDumpParser
    ref_filename = LAMMPSDUMP

    def test_creates_universe(self):
        u = mda.Universe(self.ref_filename, format='LAMMPSDUMP')

        assert isinstance(u, mda.Universe)
        assert len(u.atoms) == 24

    def test_masses_warning(self):
        # masses are mandatory, but badly guessed
        # check that user is alerted
        with self.parser(self.ref_filename) as p:
            with pytest.warns(UserWarning, match='Guessed all Masses to 1.0'):
                p.parse()

    def test_guessed_attributes(self, filename):
        u = mda.Universe(filename, format='LAMMPSDUMP')
        for attr in self.guessed_attrs:
            assert hasattr(u.atoms, attr)

    def test_id_ordering(self):
        # ids are nonsequential in file, but should get rearranged
        u = mda.Universe(self.ref_filename, format='LAMMPSDUMP')
        # the 4th in file has id==13, but should have been sorted
        assert u.atoms[3].id == 4

    def test_guessed_masses(self, filename):
        u = mda.Universe(filename, format='LAMMPSDUMP')
        expected = [1., 1., 1., 1., 1., 1., 1.]
        assert_allclose(u.atoms.masses[:7], expected)

    def test_guessed_types(self, filename):
        u = mda.Universe(filename, format='LAMMPSDUMP')
        expected = ['2', '1', '1', '2', '1', '1', '2']
        assert (u.atoms.types[:7] == expected).all()

# this tests that topology can still be constructed if non-standard or uneven
# column present.
class TestDumpParserLong(TestDumpParser):

    ref_filename = LAMMPSDUMP_long
