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

import MDAnalysis
from MDAnalysis.tests.datafiles import (
    LAMMPSdata,
    LAMMPScnt, LAMMPScnt2,
    LAMMPShyd, LAMMPShyd2,
)

from MDAnalysisTests.test_topology import _TestTopology


class _TestLammpsTop(_TestTopology):
    def test_n_bonds(self):
        if self.ref_n_bonds:
            assert_equal(len(self.universe._topology['bonds']),
                         self.ref_n_bonds)
        else:
            assert_('bonds' not in self.universe._topology)

    def test_bond_member(self):
        if self.ref_n_bonds:
            assert_(self.ref_bond in self.universe._topology['bonds'])

    def test_n_angles(self):
        if self.ref_n_angles:
            assert_equal(len(self.universe._topology['angles']),
                         self.ref_n_angles)
        else:
            assert_('angles' not in self.universe._topology)

    def test_angle_member(self):
        if self.ref_n_angles:
            assert_(self.ref_angle in self.universe._topology['angles'])

    def test_n_dihedrals(self):
        if self.ref_n_dihedrals:
            assert_equal(len(self.universe._topology['dihedrals']),
                         self.ref_n_dihedrals)
        else:
            assert_('dihedrals' not in self.universe._topology)

    def test_dihedral_member(self):
        if self.ref_n_dihedrals:
            assert_(self.ref_dihedral in self.universe._topology['dihedrals'])

    def test_n_impropers(self):
        if self.ref_n_impropers:
            assert_equal(len(self.universe._topology['impropers']),
                         self.ref_n_impropers)
        else:
            assert_('impropers' not in self.universe._topology)

    def test_improper_member(self):
        if self.ref_n_impropers:
            assert_(self.ref_improper in self.universe._topology['impropers'])


class TestLammpsData(_TestLammpsTop):
    """Tests the reading of lammps .data topology files.

    The reading of coords and velocities is done separately in test_coordinates
    """
    topology = LAMMPSdata
    parser = MDAnalysis.topology.LAMMPSParser.DATAParser
    ref_n_atoms = 18360
    ref_numresidues = 24
    ref_n_bonds = 18336
    ref_bond = (12, 14)
    ref_n_angles = 29904
    ref_angle = (3, 6, 9)
    ref_n_dihedrals = 5712
    ref_dihedral = (82, 85, 88, 89)
    ref_n_impropers = 0

    def test_charge(self):
        # No charges were supplied, should default to 0.0
        assert_equal(self.universe.atoms[0].charge, 0.0)

    def test_resid(self):
        assert_equal(len(self.universe.residues[0]), 765)

    # Testing _psf prevent building TGs
    # test length and random item from within
    def test_masses(self):
        assert_equal(self.universe.atoms[0].mass, 0.012)


class TestLAMMPSCNT(_TestLammpsTop):
    topology = LAMMPScnt
    parser = MDAnalysis.topology.LAMMPSParser.DATAParser
    ref_n_atoms = 604
    ref_numresidues = 1
    ref_n_bonds = 906
    ref_bond = (9, 467)
    ref_n_angles = 1812
    ref_angle = (17, 16, 31)
    ref_n_dihedrals = 3624
    ref_dihedral = (22, 39, 40, 41)
    ref_n_impropers = 604
    ref_improper = (210, 159, 212, 566)


class TestLAMMPSCNT2(TestLAMMPSCNT):
    topology = LAMMPScnt2
    def setUp(self):
        self.universe = MDAnalysis.Universe(self.topology,
                                            format='data')

    def test_correct_parser(self):
        pass


class TestLAMMPSHYD(_TestLammpsTop):
    topology = LAMMPShyd
    parser = MDAnalysis.topology.LAMMPSParser.DATAParser
    ref_n_atoms = 2
    ref_numresidues = 1
    ref_n_bonds = 1
    ref_bond = (0, 1)
    ref_n_angles = 0
    ref_n_dihedrals = 0
    ref_n_impropers = 0


class TestLAMMPSHYD2(TestLAMMPSHYD):
    topology = LAMMPShyd2
    def setUp(self):
        self.universe = MDAnalysis.Universe(self.topology,
                                            format='data')
    def test_correct_parser(self):
        pass

