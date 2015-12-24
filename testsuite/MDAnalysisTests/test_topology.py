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

import MDAnalysis
from MDAnalysis.core.AtomGroup import AtomGroup
from MDAnalysis.lib.distances import calc_bonds, calc_angles, calc_dihedrals
from MDAnalysis.lib.util import guess_format
from MDAnalysis.topology.core import (
    guess_atom_type, guess_atom_element, get_atom_mass,
    guess_bonds, guess_angles, guess_dihedrals, guess_improper_dihedrals,
    get_parser_for)
from MDAnalysis.core.topologyobjects import (
    TopologyGroup, TopologyObject, TopologyDict,
    Bond, Angle, Dihedral, ImproperDihedral)
from MDAnalysis.tests.datafiles import (
    PSF, DCD,
)
from MDAnalysisTests.plugins.knownfailure import knownfailure

from numpy.testing import *
from nose.plugins.attrib import attr
import numpy as np


def check_atom_type(atype, aname):
    assert_equal(guess_atom_type(aname), atype)


def check_atom_element(element, aname):
    assert_equal(guess_atom_element(aname), element)

def test_atom_element_IE():
    """Issue #476
    guess_atom_element raises IndexError when given
    name that is a single digit ('1')
    """
    assert_equal(guess_atom_element('1'), '1')


class _TestGuessAtomType(object):
    atype = None
    testnames = []
    mass = None
    element = None

    def test_guess_atom_type(self):
        for aname in self.testnames:
            yield check_atom_type, self.atype, aname

    def test_guess_atom_mass(self):
        assert_equal(get_atom_mass(self.atype), self.mass)

    def test_guess_atom_element(self):
        for aname in self.testnames:
            yield check_atom_element, self.element, aname


class TestHydrogen(_TestGuessAtomType):
    atype = 'H'
    element = 'H'
    mass = 1.008
    testnames = ['H', 'HZ', '1HZ', '2HW', 'HE']


class TestCarbon(_TestGuessAtomType):
    atype = 'C'
    element = 'C'
    mass = 12.0110
    testnames = ['C', 'CA']


class TestSodium(_TestGuessAtomType):
    atype = 'NA'
    element = 'NA'
    mass = 22.989770
    testnames = ['NA', 'NA+', 'SOD', 'QN']


class TestPotassium(_TestGuessAtomType):
    atype = 'K'
    element = 'K'
    mass = 39.102
    testnames = ['K', 'K+', 'POT', 'QK']


class TestChloride(_TestGuessAtomType):
    atype = 'CL'
    element = 'CL'
    mass = 35.450
    testnames = ['CL', 'CL-', 'CLA', 'CLAL']


class TestNitrogen(_TestGuessAtomType):
    atype = 'N'
    element = 'N'
    mass = 14.007
    testnames = ['N', 'NZ', 'NE']


class TestPhosphorous(_TestGuessAtomType):
    atype = 'P'
    element = 'P'
    mass = 30.974000
    testnames = ['P', 'PL', 'PO']


class TestSulfur(_TestGuessAtomType):
    atype = 'S'
    element = 'S'
    mass = 32.06000
    testnames = ['S', 'SG']


class TestOxygen(_TestGuessAtomType):
    atype = 'O'
    element = 'O'
    mass = 15.99900
    testnames = ['O', 'OA', 'OXT', '1OG', 'OW']


class TestCalcium(_TestGuessAtomType):
    atype = 'CA'
    element = 'CA'
    mass = 40.080000
    testnames = ['CAL', 'CA2+', 'C0']


class TestMagnesium(_TestGuessAtomType):
    atype = 'MG'
    element = 'MG'
    mass = 24.305000
    testnames = ['MG', 'MG2+']


class TestTopologyGroup_Cython(TestCase):
    """
    Check that the shortcut to all cython functions:
     - work (return proper values)
     - catch errors
    """

    def setUp(self):
        self.u = MDAnalysis.Universe(PSF, DCD)
        # topologygroups for testing
        # bond, angle, dihedral, improper
        ag = self.u.atoms[:5]
        self.bgroup = ag.bonds
        self.agroup = ag.angles
        self.tgroup = ag.dihedrals
        self.igroup = ag.impropers

    def tearDown(self):
        del self.u
        del self.bgroup
        del self.agroup
        del self.tgroup
        del self.igroup

    # bonds
    def test_wrong_type_bonds(self):
        for tg in [self.agroup, self.tgroup, self.igroup]:
            assert_raises(TypeError, tg.bonds)

    def test_right_type_bonds(self):
        assert_equal(self.bgroup.bonds(),
                     calc_bonds(self.bgroup.atom1.positions,
                                self.bgroup.atom2.positions))
        assert_equal(self.bgroup.bonds(pbc=True),
                     calc_bonds(self.bgroup.atom1.positions,
                                self.bgroup.atom2.positions,
                                box=self.u.dimensions))
        assert_equal(self.bgroup.values(),
                     calc_bonds(self.bgroup.atom1.positions,
                                self.bgroup.atom2.positions))
        assert_equal(self.bgroup.values(pbc=True),
                     calc_bonds(self.bgroup.atom1.positions,
                                self.bgroup.atom2.positions,
                                box=self.u.dimensions))

    # angles
    def test_wrong_type_angles(self):
        for tg in [self.bgroup, self.tgroup, self.igroup]:
            assert_raises(TypeError, tg.angles)

    def test_right_type_angles(self):
        assert_equal(self.agroup.angles(),
                     calc_angles(self.agroup.atom1.positions,
                                 self.agroup.atom2.positions,
                                 self.agroup.atom3.positions))
        assert_equal(self.agroup.angles(pbc=True),
                     calc_angles(self.agroup.atom1.positions,
                                 self.agroup.atom2.positions,
                                 self.agroup.atom3.positions,
                                 box=self.u.dimensions))
        assert_equal(self.agroup.values(),
                     calc_angles(self.agroup.atom1.positions,
                                 self.agroup.atom2.positions,
                                 self.agroup.atom3.positions))
        assert_equal(self.agroup.values(pbc=True),
                     calc_angles(self.agroup.atom1.positions,
                                 self.agroup.atom2.positions,
                                 self.agroup.atom3.positions,
                                 box=self.u.dimensions))

    # dihedrals & impropers
    def test_wrong_type_dihedrals(self):
        for tg in [self.bgroup, self.agroup]:
            assert_raises(TypeError, tg.dihedrals)

    def test_right_type_dihedrals(self):
        assert_equal(self.tgroup.dihedrals(),
                     calc_dihedrals(self.tgroup.atom1.positions,
                                   self.tgroup.atom2.positions,
                                   self.tgroup.atom3.positions,
                                   self.tgroup.atom4.positions))
        assert_equal(self.tgroup.dihedrals(pbc=True),
                     calc_dihedrals(self.tgroup.atom1.positions,
                                   self.tgroup.atom2.positions,
                                   self.tgroup.atom3.positions,
                                   self.tgroup.atom4.positions,
                                   box=self.u.dimensions))
        assert_equal(self.tgroup.values(),
                     calc_dihedrals(self.tgroup.atom1.positions,
                                   self.tgroup.atom2.positions,
                                   self.tgroup.atom3.positions,
                                   self.tgroup.atom4.positions))
        assert_equal(self.tgroup.values(pbc=True),
                     calc_dihedrals(self.tgroup.atom1.positions,
                                   self.tgroup.atom2.positions,
                                   self.tgroup.atom3.positions,
                                   self.tgroup.atom4.positions,
                                   box=self.u.dimensions))

    def test_right_type_impropers(self):
        assert_equal(self.igroup.dihedrals(),
                     calc_dihedrals(self.igroup.atom1.positions,
                                   self.igroup.atom2.positions,
                                   self.igroup.atom3.positions,
                                   self.igroup.atom4.positions))
        assert_equal(self.igroup.dihedrals(pbc=True),
                     calc_dihedrals(self.igroup.atom1.positions,
                                   self.igroup.atom2.positions,
                                   self.igroup.atom3.positions,
                                   self.igroup.atom4.positions,
                                   box=self.u.dimensions))
        assert_equal(self.igroup.values(),
                     calc_dihedrals(self.igroup.atom1.positions,
                                   self.igroup.atom2.positions,
                                   self.igroup.atom3.positions,
                                   self.igroup.atom4.positions))
        assert_equal(self.igroup.values(pbc=True),
                     calc_dihedrals(self.igroup.atom1.positions,
                                   self.igroup.atom2.positions,
                                   self.igroup.atom3.positions,
                                   self.igroup.atom4.positions,
                                   box=self.u.dimensions))


class TestTopologyGuessers(TestCase):
    """Test the various ways of automating topology creation in the Universe

    guess_bonds
    guess_angles
    guess_dihedrals
    guess_improper_dihedrals
    """

    def setUp(self):
        self.u = MDAnalysis.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    # guess_bonds
    def test_guess_bonds_wronginput(self):  # give too few coords and watch it explode
        assert_raises(ValueError, guess_bonds, self.u.atoms, self.u.atoms.positions[:-2])

    def test_guess_bonds_badtype(self):
        # rename my carbons and watch it get confused about missing types
        self.u.select_atoms('type C').set_types('QQ')
        assert_raises(ValueError, guess_bonds, self.u.atoms, self.u.atoms.positions)

    def test_guess_bonds_withag(self):
        # here's one I prepared earlier
        bondlist = (
            (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3),
            (1, 4), (2, 3), (2, 4), (2, 8), (3, 4), (4, 5), (4, 6),
            (6, 7), (6, 8), (6, 9), (7, 8))
        user_vdw = {'56': 1.5, '2': 1.5, '22': 1.5, '6': 1.5, '23': 1.5, '3': 1.5}

        assert_equal(guess_bonds(self.u.atoms[:10],
                                 self.u.atoms.positions[:10],
                                 vdwradii=user_vdw),
                     bondlist)

    def test_guess_bonds_withlist(self):
        bondlist = (
            (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3),
            (1, 4), (2, 3), (2, 4), (2, 8), (3, 4), (4, 5), (4, 6),
            (6, 7), (6, 8), (6, 9), (7, 8))
        user_vdw = {'56': 1.5, '2': 1.5, '22': 1.5, '6': 1.5, '23': 1.5, '3': 1.5}

        assert_equal(guess_bonds(list(self.u.atoms[:10]),
                                 self.u.atoms.positions[:10],
                                 vdwradii=user_vdw),
                     bondlist)

    def test_guess_angles(self):
        ag = self.u.atoms[:10]

        # Guess angles based on bond information
        guessed_angs = guess_angles(ag.bonds)

        # This gets a TG to write results back out in tuple of tuples formats
        dump_result = ag.angles.dump_contents()

        assert_equal(set(guessed_angs), set(dump_result))

    def test_guess_dihedrals(self):
        ag = self.u.atoms[:10]
        ag.bonds
        ag.angles

        guessed_tors = guess_dihedrals(ag.angles)

        dump_result = ag.dihedrals.dump_contents()

        assert_equal(set(guessed_tors), set(dump_result))

    def test_guess_improper_dihedrals(self):
        ag = self.u.atoms[:5]
        ag.bonds
        ag.angles

        # Group of angles to work off
        angs = self.u.atoms.angles.atomgroup_intersection(ag, strict=True)

        # Pre calculated reference result
        result = (
            (0, 4, 3, 2), (0, 3, 1, 4), (0, 4, 1, 3),
            (0, 2, 1, 4), (0, 4, 2, 1), (0, 4, 1, 2),
            (0, 3, 2, 4), (0, 3, 2, 1), (0, 4, 3, 1),
            (0, 3, 1, 2), (0, 4, 2, 3), (0, 2, 1, 3))

        imps = guess_improper_dihedrals(angs)

        assert_equal(set(result), set(imps))

    # Functional testing for these functions
    # Test that Universe accepts their output as input
    def test_guess_angles_set(self):
        self.u.atoms.angles = guess_angles(self.u.atoms.bonds)

        assert_equal(len(self.u.atoms.angles), 6123)

    def test_guess_dihedrals_set(self):
        self.u.dihedrals = guess_dihedrals(self.u.atoms.angles)

        assert_equal(len(self.u.dihedrals), 8921)

    def test_guess_impropers_set(self):
        self.u.impropers = guess_improper_dihedrals(self.u.atoms.angles)

        assert_equal(len(self.u.impropers), 10314)
