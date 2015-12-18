# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from numpy.testing import (
    assert_,
    assert_equal,
)

import MDAnalysis as mda
from MDAnalysisTests.topology.base import ParserBase
from MDAnalysis.tests.datafiles import (
    LAMMPSdata,
    LAMMPScnt, LAMMPScnt2,
    LAMMPShyd, LAMMPShyd2,
)



class LammpsBase(ParserBase):
    parser = mda.topology.LAMMPSParser.DATAParser
    expected_n_segments = 1
    expected_attrs = ['types', 'resids', 'masses', 'charges']

    def test_n_bonds(self):
        if self.ref_n_bonds:
            assert_equal(len(self.top.bonds.values),
                         self.ref_n_bonds)
        else:
            assert_(not hasattr(self.top, 'bonds'))

    def test_bond_member(self):
        if self.ref_n_bonds:
            assert_(self.ref_bond in self.top.bonds.values)

    def test_n_angles(self):
        if self.ref_n_angles:
            assert_equal(len(self.top.angles.values),
                         self.ref_n_angles)
        else:
            assert_(not hasattr(self.top, 'angles'))

    def test_angle_member(self):
        if self.ref_n_angles:
            assert_(self.ref_angle in self.top.angles.values)

    def test_n_dihedrals(self):
        if self.ref_n_dihedrals:
            assert_equal(len(self.top.dihedrals.values),
                         self.ref_n_dihedrals)
        else:
            assert_(not hasattr(self.top, 'dihedrals'))

    def test_dihedral_member(self):
        if self.ref_n_dihedrals:
            assert_(self.ref_dihedral in self.top.dihedrals.values)

    def test_n_impropers(self):
        if self.ref_n_impropers:
            assert_equal(len(self.top.impropers.values),
                         self.ref_n_impropers)
        else:
            assert_(not hasattr(self.top, 'impropers'))

    def test_improper_member(self):
        if self.ref_n_impropers:
            assert_(self.ref_improper in self.top.impropers.values)


class TestLammpsData(LammpsBase):
    """Tests the reading of lammps .data topology files.

    The reading of coords and velocities is done separately in
    test_coordinates
    """
    filename = LAMMPSdata
    expected_n_atoms = 18360
    expected_n_residues = 24
    ref_n_bonds = 18336
    ref_bond = (12, 14)
    ref_n_angles = 29904
    ref_angle = (3, 6, 9)
    ref_n_dihedrals = 5712
    ref_dihedral = (82, 85, 88, 89)
    ref_n_impropers = 0


class TestLAMMPSCNT(LammpsBase):
    filename = LAMMPScnt
    expected_n_atoms = 604
    expected_n_residues = 1
    ref_n_bonds = 906
    ref_bond = (9, 467)
    ref_n_angles = 1812
    ref_angle = (17, 16, 31)
    ref_n_dihedrals = 3624
    ref_dihedral = (22, 39, 40, 41)
    ref_n_impropers = 604
    ref_improper = (210, 159, 212, 566)


class TestLAMMPSCNT2(TestLAMMPSCNT):
    filename = LAMMPScnt2


class TestLAMMPSHYD(LammpsBase):
    filename = LAMMPShyd
    expected_n_atoms = 2
    expected_n_residues = 1
    ref_n_bonds = 1
    ref_bond = (0, 1)
    ref_n_angles = 0
    ref_n_dihedrals = 0
    ref_n_impropers = 0


class TestLAMMPSHYD2(TestLAMMPSHYD):
    filename = LAMMPShyd2
