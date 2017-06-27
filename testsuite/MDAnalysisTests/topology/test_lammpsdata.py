# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

from numpy.testing import (
    assert_,
    assert_equal,
)
import numpy as np

import MDAnalysis as mda
from MDAnalysisTests.topology.base import ParserBase
from MDAnalysis.tests.datafiles import (
    LAMMPSdata,
    LAMMPScnt, LAMMPScnt2,
    LAMMPShyd, LAMMPShyd2,
    LAMMPSdata_deletedatoms,
)



class LammpsBase(ParserBase):
    parser = mda.topology.LAMMPSParser.DATAParser
    expected_n_segments = 1
    expected_attrs = ['types', 'resids', 'masses', 'charges']

    def test_n_atom_types(self):
        assert_equal(len(set(self.top.types.values)), self.expected_n_atom_types)

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

    def test_creates_universe(self):
        u = mda.Universe(self.filename, format='DATA')


class TestLammpsData(LammpsBase):
    """Tests the reading of lammps .data topology files.

    The reading of coords and velocities is done separately in
    test_coordinates
    """
    filename = LAMMPSdata
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
    filename = LAMMPScnt
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


class TestLAMMPSCNT2(TestLAMMPSCNT):
    filename = LAMMPScnt2


class TestLAMMPSHYD(LammpsBase):
    filename = LAMMPShyd
    expected_n_atoms = 2
    expected_n_atom_types = 1
    expected_n_residues = 1
    ref_n_bonds = 1
    ref_bond = (0, 1)
    ref_n_angles = 0
    ref_n_dihedrals = 0
    ref_n_impropers = 0


class TestLAMMPSHYD2(TestLAMMPSHYD):
    filename = LAMMPShyd2


class TestLAMMPSDeletedAtoms(LammpsBase):
    filename = LAMMPSdata_deletedatoms

    expected_n_atoms = 10
    expected_n_atom_types = 2
    expected_n_residues = 1

    ref_n_bonds = 9
    ref_bond = (0, 3)
    ref_n_angles = 0
    ref_n_dihedrals = 0
    ref_n_impropers = 0

    def test_atom_ids(self):
        u = mda.Universe(self.filename)

        assert_equal(u.atoms.ids,
                     [1, 10, 1002, 2003, 2004, 2005, 2006, 2007, 2008, 2009])

    def test_traj(self):
        u = mda.Universe(self.filename)

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

