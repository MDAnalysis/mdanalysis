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
    PRMpbc, PRM12, PSF, PSF_NAMD, PSF_nosegid, DMS, PDB_small, DCD,
    LAMMPSdata, trz4data, TPR, PDB, XYZ_mini, GMS_SYMOPT, GMS_ASYMSURF,
    DLP_CONFIG, DLP_CONFIG_order, DLP_CONFIG_minimal,
    DLP_HISTORY, DLP_HISTORY_order, DLP_HISTORY_minimal, HoomdXMLdata)
from MDAnalysisTests.plugins.knownfailure import knownfailure

from numpy.testing import *
from nose.plugins.attrib import attr
import numpy as np


def check_atom_type(atype, aname):
    assert_equal(guess_atom_type(aname), atype)


def check_atom_element(element, aname):
    assert_equal(guess_atom_element(aname), element)


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


# add more...

# specific topology readers
# add more!

# CHARMM and NAMD PSF

class _TestTopology(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(self.topology)

    def tearDown(self):
        del self.universe

    def test_correct_parser(self):
        """Check that get_parser returns the intended parser"""
        try:
            perm = self.perm
        except AttributeError:
            perm = False
        ret = get_parser_for(self.topology, permissive=perm)

        assert_equal(self.parser, ret)

    def test_parser(self):
        """Check that the parser works as intended,
        and that the returned value is a dictionary
        """
        with self.parser(self.topology) as p:
            ret = p.parse()
        assert_equal(type(ret), type(dict()))

    #    def test_parser_raises_IOE(self):
    #        """Check that when given junk input, they raise IOError"""
    #        p = self.parser(trz4data)
    #        assert_raises(IOError, p.parse) or assert_raises(ValueError, p.parse)

    def test_parser_atoms(self):
        assert_equal(self.universe.atoms.n_atoms,
                     self.ref_n_atoms,
                     "wrong number of atoms in topology")
        assert_equal(self.universe.atoms.n_residues,
                     self.ref_numresidues,
                     "wrong number of residues in topology")

    def test_atom_number(self):
        assert_equal(self.universe.atoms[0].index, 0,
                     "first atom should have Atom.index 0")
        assert_equal(self.universe.atoms[-1].index, self.ref_n_atoms - 1,
                     "last atom has wrong Atom.index")


class RefAdKSmall(object):
    """Mixin class to provide comparison numbers.

    Based on small PDB with AdK (:data:`PDB_small`).
    """
    topology = PSF
    parser = MDAnalysis.topology.PSFParser.PSFParser
    ref_n_atoms = 3341
    ref_numresidues = 214


class TestPSF_CHARMM_STANDARD(_TestTopology, RefAdKSmall):
    """Testing CHARMM standard PSF file format"""


class RefNAMD_CGENFF(object):
    """Testfiles provided by JiyongPark77.

    NAMD/VMD XPLOR-style PSF file (using CGENFF residues/atoms).

    https://github.com/MDAnalysis/mdanalysis/issues/107
    """
    topology = PSF_NAMD
    parser = MDAnalysis.topology.PSFParser.PSFParser
    ref_n_atoms = 130
    ref_numresidues = 6

class RefHoomdXML(object):
    """Hoomd XML test data and answers
    """
    topology = HoomdXMLdata
    parser = MDAnalysis.topology.HoomdXMLParser.HoomdXMLParser
    ref_n_atoms = 769
    ref_numresidues = 1

class TestHoomdXML(_TestTopology, RefHoomdXML):
    """Testing Hoomd XML file"""
    pass

class TestPSF_NAMD_CGENFF(_TestTopology, RefNAMD_CGENFF):
    """Testing NAMD PSF file (with CGENFF atom types, Issue 107)"""


class TestPSF_Issue121(TestCase):
    @attr('issue')
    def test_nosegid(self):
        try:
            u = MDAnalysis.Universe(PSF_nosegid)
        except IndexError:
            raise AssertionError("Issue 121 not fixed: cannot load PSF with empty SEGID")
        assert_equal(u.atoms.n_atoms, 98)
        assert_equal(u.segments.segids, ["SYSTEM"])


class TestPSF_bonds(TestCase):
    """Tests reading of bonds angles and dihedrals in psf files"""

    def setUp(self):
        topology = PSF
        self.universe = MDAnalysis.Universe(topology)
        self.universe.build_topology()

    def tearDown(self):
        del self.universe

    def test_bonds_counts(self):
        assert_equal(len(self.universe._topology['bonds']), 3365)
        assert_equal(len(self.universe.atoms[0].bonds), 4)
        assert_equal(len(self.universe.atoms[42].bonds), 1)

    def test_bonds_identity(self):
        a1 = self.universe.atoms[0]
        a2 = self.universe.atoms[1]
        a3 = self.universe.atoms[2]
        a4 = self.universe.atoms[3]
        a5 = self.universe.atoms[4]
        a42 = self.universe.atoms[41]
        # Bonds might change order, so use any checks through bond list
        assert_equal(all([a1 in b for b in a1.bonds]), True)  # check all bonds have this atom
        assert_equal(any([a2 in b for b in a1.bonds]), True)  # then check certain atoms are present
        assert_equal(any([a3 in b for b in a1.bonds]), True)
        assert_equal(any([a4 in b for b in a1.bonds]), True)
        assert_equal(any([a5 in b for b in a1.bonds]), True)
        assert_equal(any([a42 in b for b in a1.bonds]), False)  # and check everything isn't True

    def test_angles_counts(self):
        assert_equal(len(self.universe._topology['angles']), 6123)
        assert_equal(len(self.universe.atoms[0].angles), 9)
        assert_equal(len(self.universe.atoms[42].angles), 2)

    def test_angles_identity(self):
        a1 = self.universe.atoms[0]
        a2 = self.universe.atoms[1]
        a3 = self.universe.atoms[2]
        a5 = self.universe.atoms[4]
        a6 = self.universe.atoms[5]
        a42 = self.universe.atoms[41]
        assert_equal(all([a1 in b for b in a1.angles]), True)
        assert_equal(any([a2 in b and a3 in b for b in a1.angles]), True)
        assert_equal(any([a5 in b and a6 in b for b in a1.angles]), True)
        assert_equal(any([a42 in b for b in a1.angles]), False)
        assert_equal(any([a2 in b and a6 in b for b in a1.angles]),
                     False)  # both a2 and a6 feature, but never simultaneously

    def test_dihedrals_counts(self):
        assert_equal(len(self.universe._topology['dihedrals']), 8921)
        assert_equal(len(self.universe.atoms[0].dihedrals), 14)

    def test_dihedrals_identity(self):
        a1 = self.universe.atoms[0]
        a2 = self.universe.atoms[1]
        a3 = self.universe.atoms[2]
        a5 = self.universe.atoms[4]
        a6 = self.universe.atoms[5]
        a7 = self.universe.atoms[6]
        a42 = self.universe.atoms[41]
        assert_equal(all([a1 in b for b in a1.dihedrals]), True)
        assert_equal(any([a2 in b and a5 in b and a6 in b for b in a1.dihedrals]), True)
        assert_equal(any([a2 in b and a5 in b and a7 in b for b in a1.dihedrals]), True)
        assert_equal(any([a42 in b for b in a1.dihedrals]), False)
        assert_equal(any([a2 in b and a3 in b and a6 in b for b in a1.dihedrals]), False)


class TestTopologyObjects(TestCase):
    """Test the base TopologyObject funtionality

    init
    repr
    eq
    ne
    iter
    len
    """

    def setUp(self):
        self.precision = 3  # rather lenient but see #271
        self.u = MDAnalysis.Universe(PSF, DCD)
        self.u.build_topology()
        self.a1 = list(self.u.atoms[1:3])
        self.a2 = list(self.u.atoms[3:5])

        self.TO1 = TopologyObject(self.a1)
        self.TO2 = TopologyObject(self.a2)
        # this atom only has one bond, so the order bonds are done (random because of sets)
        # won't come back to bite us
        self.b = self.u.atoms[12].bonds[0]

    def tearDown(self):
        del self.u
        del self.a1
        del self.a2
        del self.TO1
        del self.TO2
        del self.b

    def test_repr(self):
        assert_equal(repr(self.TO1),
                     '<TopologyObject between: Atom 2 (HT1 of MET-1), Atom 3 (HT2 of MET-1)>')

    def test_eq(self):
        TO1_b = TopologyObject(self.a1)

        assert_equal(self.TO1 == TO1_b, True)
        assert_equal(self.TO1 == self.TO2, False)

    def test_ne(self):
        TO1_b = TopologyObject(self.a1)

        assert_equal(self.TO1 != TO1_b, False)
        assert_equal(self.TO1 != self.TO2, True)

    def test_gt(self):
        assert_equal(self.TO1 > self.TO2, False)

    def test_lt(self):
        assert_equal(self.TO1 < self.TO2, True)

    def test_iter(self):
        assert_equal(self.a1, list(self.TO1))

    def test_len(self):
        assert_equal(len(self.a1), 2)

    def test_indices(self):
        assert_equal(self.b.indices, tuple([b.index for b in self.b.atoms]))

    # Bond class checks
    def test_partner(self):
        a1, a2 = self.b
        a3 = self.u.atoms[0]

        assert_equal(self.b.partner(a1), a2)
        assert_equal(self.b.partner(a2), a1)
        assert_raises(ValueError, self.b.partner, a3)

    def test_isguessed(self):
        b = Bond([self.u.atoms[0], self.u.atoms[12]])
        b.is_guessed = True

        assert_equal(b.is_guessed, True)

    def test_bondlength(self):
        assert_almost_equal(self.b.length(), 1.7661301556941993, self.precision)

    def test_bondrepr(self):
        assert_equal(repr(self.b),
                     '<Bond between: Atom 10 (CG of MET 1 None) and Atom 13 (SD of MET1 None), length 1.77 A>')

    # Angle class checks
    def test_angle(self):
        angle = self.u.atoms[210].angles[0]

        assert_almost_equal(angle.angle(), 107.20893, self.precision)
        assert_almost_equal(angle.value(), 107.20893, self.precision)

    # Dihedral class check
    def test_dihedral(self):
        dihedral = self.u.atoms[14].dihedrals[0]

        assert_almost_equal(dihedral.dihedral(), 18.317778, self.precision)
        assert_almost_equal(dihedral.value(), 18.317778, self.precision)

    # Improper_Dihedral class check
    def test_improper(self):
        imp = self.u.atoms[4].impropers[0]

        assert_almost_equal(imp.improper(), -3.8370631, self.precision)
        assert_almost_equal(imp.value(), -3.8370631, self.precision)


class TestTopologyGroup(TestCase):
    """Tests TopologyDict and TopologyGroup classes with psf input"""

    def setUp(self):
        topology = PSF
        self.universe = MDAnalysis.Universe(topology)
        # force the loading of topology
        self.universe.build_topology()
        self.res1 = self.universe.residues[0]
        self.res2 = self.universe.residues[1]
        # topologydicts for testing
        self.b_td = self.universe.bonds.topDict
        self.a_td = self.universe.angles.topDict
        self.t_td = self.universe.dihedrals.topDict

    def tearDown(self):
        del self.universe
        del self.res1
        del self.res2
        del self.b_td
        del self.a_td
        del self.t_td

    # Checking TopologyDict functionality
    # * check that enough types are made
    # * check the identity of one of the keys
    # * check uniqueness of keys (ie reversed doesnt exist)
    # * then check same key reversed is accepted
    # * then select based on key
    # * then select based on reversed key and check output is the same
    # all for Bonds Angles and Dihedrals
    def test_td_len(self):
        assert_equal(len(self.b_td), 57)

    def test_td_iter(self):
        assert_equal(list(self.b_td), list(self.b_td.dict.keys()))

    def test_td_keyerror(self):
        assert_raises(KeyError, self.b_td.__getitem__, ('something', 'stupid'))

    def test_bonds_types(self):
        """Tests TopologyDict for bonds"""
        assert_equal(len(self.universe.atoms.bonds.types()), 57)
        assert_equal(len(self.res1.bonds.types()), 12)

    def test_bonds_contains(self):
        assert_equal(('57', '2') in self.b_td, True)

    def test_bond_uniqueness(self):
        bondtypes = self.universe.bonds.types()
        # check that a key doesn't appear in reversed format in keylist
        # have to exclude case of b[::-1] == b as this is false positive
        assert_equal(any([b[::-1] in bondtypes for b in bondtypes if b[::-1] != b]),
                     False)

    def test_bond_reversal(self):
        bondtypes = self.universe.bonds.types()
        b = bondtypes[1]
        assert_equal(all([b in self.b_td, b[::-1] in self.b_td]), True)

        tg1 = self.b_td[b]
        tg2 = self.b_td[b[::-1]]
        assert_equal(tg1 == tg2, True)

    def test_angles_types(self):
        """TopologyDict for angles"""
        assert_equal(len(self.universe.atoms.angles.types()), 130)

    def test_angles_contains(self):
        assert_equal(('23', '73', '1') in self.a_td, True)

    def test_angles_uniqueness(self):
        bondtypes = self.a_td.keys()
        assert_equal(any([b[::-1] in bondtypes for b in bondtypes if b[::-1] != b]),
                     False)

    def test_angles_reversal(self):
        bondtypes = self.a_td.keys()
        b = bondtypes[1]
        assert_equal(all([b in self.a_td, b[::-1] in self.a_td]), True)

        tg1 = self.a_td[b]
        tg2 = self.a_td[b[::-1]]
        assert_equal(tg1 == tg2, True)

    def test_dihedrals_types(self):
        """TopologyDict for dihedrals"""
        assert_equal(len(self.universe.atoms.dihedrals.types()), 220)

    def test_dihedrals_contains(self):
        assert_equal(('30', '29', '20', '70') in self.t_td, True)

    def test_dihedrals_uniqueness(self):
        bondtypes = self.t_td.keys()
        assert_equal(any([b[::-1] in bondtypes for b in bondtypes if b[::-1] != b]),
                     False)

    def test_dihedrals_reversal(self):
        bondtypes = self.t_td.keys()
        b = bondtypes[1]
        assert_equal(all([b in self.t_td, b[::-1] in self.t_td]), True)

        tg1 = self.t_td[b]
        tg2 = self.t_td[b[::-1]]
        assert_equal(tg1 == tg2, True)

    def test_bad_creation(self):
        """Test making a TopologyDict out of nonsense"""
        inputlist = ['a', 'b', 'c']
        assert_raises(TypeError, TopologyDict, inputlist)

    def test_bad_creation_TG(self):
        """Test making a TopologyGroup out of nonsense"""
        inputlist = ['a', 'b', 'c']
        assert_raises(TypeError, TopologyGroup, inputlist)

    def test_TG_equality(self):
        """Make two identical TGs,
        * check they're equal
        * change one very slightly and see if they notice
        """
        tg = self.universe.bonds.selectBonds(('23', '3'))
        tg2 = self.universe.bonds.selectBonds(('23', '3'))

        assert_equal(tg == tg2, True)

        tg3 = self.universe.bonds.selectBonds(('81', '10'))
        assert_equal(tg == tg3, False)
        assert_equal(tg != tg3, True)

    def test_create_TopologyGroup(self):
        res1_tg = self.res1.bonds.selectBonds(('23', '3'))  # make a tg
        assert_equal(len(res1_tg), 4)  # check size of tg
        testbond = self.universe.atoms[7].bonds[0]
        assert_equal(testbond in res1_tg, True)  # check a known bond is present

        res1_tg2 = self.res1.bonds.selectBonds(('23', '3'))
        assert_equal(res1_tg == res1_tg2, True)

    def test_create_empty_TG(self):
        tg = TopologyGroup([])

        def check(a):
            if a:
                return True
            else:
                return False

        assert_equal(check(tg), False)
        assert_equal(len(tg), 0)
        assert_equal(tg.toptype, None)

    # Loose TG intersection
    def test_TG_loose_intersection(self):
        """Pull bonds from a TG which are at least partially in an AG"""

        def check_loose_intersection(topg, atomg):
            return all([any([a in ag for a in b.atoms]) for b in topg.bondlist])

        def manual(topg, atomg):
            man = []
            for b in topg.bondlist:
                if any([a in atomg for a in b.atoms]):
                    man.append(b)
            return TopologyGroup(man)

        u = self.universe
        ag = self.universe.atoms[10:60]

        # Check that every bond has at least one atom in the atomgroup
        for TG in [u.bonds, u.angles, u.dihedrals, u.impropers]:
            newTG = TG.atomgroup_intersection(ag)

            assert_equal(check_loose_intersection(newTG, ag), True,
                         err_msg="Loose intersection failed with: " + TG.toptype)
            assert_equal(manual(TG, ag), newTG,
                         err_msg="Loose intersection failed with: " + TG.toptype)

    # Strict TG intersection
    def test_TG_strict_intersection(self):
        """Pull bonds from TG which are fully in an AG"""

        def check_strict_intersection(topg, atomg):
            new_topg = topg.atomgroup_intersection(atomg, strict=True)

            return all([all([a in atomg for a in b.atoms]) for b in new_topg.bondlist])

        def manual(topg, atomg):
            if len(atomg) == 1:  # hack for Atom input
                atomg = [atomg]
            man = []
            for b in topg.bondlist:
                if all([a in atomg for a in b.atoms]):
                    man.append(b)

            if len(man) > 0:
                return TopologyGroup(man)
            else:
                return None

        u = self.universe
        testinput = self.universe.atoms[10:60]

        # bonds
        assert_equal(check_strict_intersection(u.bonds, testinput), True)
        assert_equal(manual(u.bonds, testinput),
                     u.bonds.atomgroup_intersection(testinput, strict=True))
        # angles
        assert_equal(check_strict_intersection(u.angles, testinput), True)
        assert_equal(manual(u.angles, testinput),
                     u.angles.atomgroup_intersection(testinput, strict=True))
        # dihedrals
        assert_equal(check_strict_intersection(u.dihedrals, testinput), True)
        assert_equal(manual(u.dihedrals, testinput),
                     u.dihedrals.atomgroup_intersection(testinput, strict=True))

    def test_verticalTG(self):
        b1 = self.universe.atoms[0].dihedrals[0]
        b2 = self.universe.atoms[20].dihedrals[0]

        TG = TopologyGroup([b1, b2])

        forwards = [AtomGroup([b1[i], b2[i]]) for i in range(4)]
        backwards = [AtomGroup([b2[i], b1[i]]) for i in range(4)]

        verts = [TG.atom1, TG.atom2, TG.atom3, TG.atom4]

        # the lists might be in one of two formats, but always in a strict order
        # ie any(1234 or 4321) but not (1324)
        assert_equal(any([
            all([list(x) == list(y) for x, y in zip(forwards, verts)]),
            all([list(x) == list(y) for x, y in zip(backwards, verts)])]), True)

    def test_add_TopologyGroups(self):
        res1_tg = self.res1.bonds.selectBonds(('23', '3'))
        res2_tg = self.res2.bonds.selectBonds(('23', '3'))

        combined_tg = res1_tg + res2_tg  # add tgs together
        assert_equal(len(combined_tg), 10)

        big_tg = self.universe.atoms.bonds.selectBonds(('23', '3'))

        big_tg += combined_tg  # try and add some already included bonds
        assert_equal(len(big_tg), 494)  # check len doesn't change

    def test_add_singleitem(self):
        tg = self.universe.bonds[:10]
        to = self.universe.bonds[55]

        assert_equal(len(tg + to), 11)

    def test_add_wrongtype_TopologyGroup(self):
        def adder(a, b):
            return a + b

        tg = self.universe.bonds[:10]  # TG of bonds

        ang = self.universe.angles[10]  # single angle
        angg = self.universe.angles[:10]  # TG of angles

        for other in [ang, angg]:
            assert_raises(TypeError, adder, tg, other)

    def test_bad_add_TopologyGroup(self):
        def adder(a, b):
            return a + b

        tg = self.universe.bonds[:10]  # TopologyGroup

        ag = self.universe.atoms[:10]  # AtomGroup

        assert_raises(TypeError, adder, tg, ag)

    def test_TG_getitem_single(self):
        tg = self.universe.bonds[:10]

        bondlist = list(tg)
        bond = tg[0]

        assert_equal(bond, bondlist[0])

    def test_TG_getitem_slice(self):
        tg = self.universe.bonds[:10]

        tg2 = tg[1:4]

        assert_equal(list(tg2), tg.bondlist[1:4])

    def test_TG_getitem_fancy(self):
        tg = self.universe.bonds[:10]

        tg2 = tg[[1, 4, 5]]

        manual = TopologyGroup([tg[i] for i in [1, 4, 5]])

        assert_equal(list(tg2), list(manual))

    def test_TG_getitem_bool(self):
        # Issue #282
        sel = np.array([True, False, True])
        tg = self.universe.bonds[10:30]
        tg2 = tg[sel]
        assert_equal(len(tg2), 2)
        for b in [tg[0], tg[2]]:
            assert_equal(b in tg2, True)

    def test_TG_getitem_bool_IE(self):
        sel = []
        tg = self.universe.bonds[10:13]
        tg2 = tg[sel]
        assert_equal(len(tg2), 0)

    def test_TG_dumpconts(self):
        """
        this function tries to spit back out the tuple of tuples that made it

        because bonds are sorted before initialisation, sometimes this list
        has bonds that are backwards compared to the input, hence all the hacking
        in this test.
        """
        inpt = self.universe._topology['bonds']  # what we started with
        inpt = [tuple(sorted(a)) for a in inpt]  # make sure each entry is sorted

        dump = self.universe.bonds.dump_contents()

        assert_equal(set(dump), set(inpt))

    def test_TG_indices_creation(self):
        """Create a TG from indices"""
        bonds = [(0, 1), (1, 2)]

        tg = TopologyGroup.from_indices(bonds, self.universe.atoms,
                                        bondclass=Bond)

        assert_equal(len(tg), 2)
        b1 = self.universe.atoms[[0, 1]].bond
        b2 = self.universe.atoms[[1, 2]].bond
        assert_(b1 in tg)
        assert_(b2 in tg)
        assert_equal(bonds, tg.to_indices())

    def test_TG_from_indices_roundtrip(self):
        """Round trip check of dumping indices then recreating"""
        tg = self.universe.bonds[:10]
        idx = tg.to_indices()

        tg2 = TopologyGroup.from_indices(idx, self.universe.atoms,
                                         bondclass=Bond)

        # This doesn't work as it uses set operation and .from_indices
        # has created new Bond instances
        # assert_equal(tg, tg2)
        # instead...
        assert_equal(len(tg), len(tg2))
        assert_equal(tg.to_indices(), tg2.to_indices())

    def test_TG_without_bondclass(self):
        ag = self.universe.atoms

        tg = TopologyGroup.from_indices([(0, 1), (1, 2), (2, 3)], ag)
        assert_(len(tg) == 3)
        assert_(type(tg[0]) == Bond)

        tg = TopologyGroup.from_indices([(0, 1, 2), (1, 2, 3), (2, 3, 4)],
                                        ag)
        assert_(len(tg) == 3)
        assert_(type(tg[0]) == Angle)

        tg = TopologyGroup.from_indices([(0, 1, 2, 3), (1, 2, 3, 4)],
                                        ag)
        assert_(len(tg) == 2)
        assert_(type(tg[0]) == Dihedral)

    def test_force_bondclass(self):
        # Make a TG of improper dihedral
        tg = TopologyGroup.from_indices([(0, 1, 2, 3), (2, 3, 4, 6)],
                                        self.universe.atoms,
                                        bondclass=ImproperDihedral)
        assert_(type(tg[0]) == ImproperDihedral)

    def test_from_indices_nonglobal_idx(self):
        idx = [(0, 1), (4, 5)]

        ag = self.universe.atoms[100:]

        tg = TopologyGroup.from_indices(idx, ag)

        b1 = self.universe.atoms[[100, 101]].bond
        b2 = self.universe.atoms[[104, 105]].bond

        assert_(b1 in tg)
        assert_(b2 in tg)

    def test_from_indices_VE(self):
        idx = [(0, 1, 2, 3, 4), (2, 3, 4, 5, 6)]

        assert_raises(ValueError, TopologyGroup.from_indices,
                      idx, self.universe.atoms)


class TestTopologyGroup_Cython(TestCase):
    """
    Check that the shortcut to all cython functions:
     - work (return proper values)
     - catch errors
    """

    def setUp(self):
        self.u = MDAnalysis.Universe(PSF, DCD)
        self.u.build_topology()
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


# AMBER
class RefCappedAla(object):
    """Mixin class to provide comparison numbers.

    Capped Ala in water
    """
    topology = PRMpbc
    parser = MDAnalysis.topology.TOPParser.TOPParser
    ref_n_atoms = 5071
    ref_numresidues = 1686
    ref_proteinatoms = 22


class RefAMBER12(object):
    """Fixture data for testing AMBER12 reading (Issue 100)"""
    topology = PRM12
    parser = MDAnalysis.topology.TOPParser.TOPParser
    ref_n_atoms = 8923
    ref_numresidues = 2861
    ref_proteinatoms = 0


class TestAMBER(_TestTopology, RefCappedAla):
    """Testing AMBER PRMTOP parser (Issue 76)"""
    # note: hard to test the issue because one needs a very specifi datafile
    #       so this test really checks that we did not break the parser for the
    #       existing test cases


class TestAMBER12(_TestTopology, RefAMBER12):
    """Testing AMBER 12 PRMTOP parser (Issue 100)"""


# PDB

class RefPDB(object):
    topology = PDB_small
    parser = MDAnalysis.topology.PDBParser.PDBParser
    ref_n_atoms = 3341
    ref_numresidues = 214


class RefPDB_Perm(RefPDB):
    perm = True
    parser = MDAnalysis.topology.PrimitivePDBParser.PrimitivePDBParser


class TestPDB(_TestTopology, RefPDB):
    """Testing PDB topology parsing (PrimitivePDB)"""
    pass


class TestPDB_Perm(_TestTopology, RefPDB_Perm):
    pass


class RefXPDB(object):
    topology = PDB
    parser = MDAnalysis.topology.ExtendedPDBParser.ExtendedPDBParser
    ref_n_atoms = 47681
    ref_numresidues = 11302


# DESRES
class RefDMS(object):
    topology = DMS
    parser = MDAnalysis.topology.DMSParser.DMSParser
    ref_n_atoms = 3341
    ref_numresidues = 214


class TestDMSReader(_TestTopology, RefDMS):
    def test_number_of_bonds(self):
        # Desired value taken from VMD
        #      Info)    Atoms: 3341
        assert_equal(len(self.universe.bonds), 3365)

    def test_bond_order(self):
        pass

    def test_segid(self):
        segid = set([a.segid for a in self.universe.atoms])
        assert_equal(segid, set(("4AKE",)))

    def test_atomsels(self):
        # Desired value taken from VMD atomsel
        s0 = self.universe.select_atoms("name CA")
        assert_equal(len(s0), 214)

        s1 = self.universe.select_atoms("resid 33")
        assert_equal(len(s1), 12)

        s2 = self.universe.select_atoms("segid 4AKE")
        assert_equal(len(s2), 3341)

        s3 = self.universe.select_atoms("resname ALA")
        assert_equal(len(s3), 190)

    def test_atom_number(self):
        assert_equal(self.universe.atoms[0].index, 0,
                     "first atom should have Atom.index 0")
        assert_equal(self.universe.atoms[-1].index, 3341 - 1,
                     "last atom has wrong Atom.index")


# GROMACS TPR
# see also test_tprparser
class RefTPR(object):
    parser = MDAnalysis.topology.TPRParser.TPRParser
    topology = TPR
    ref_n_atoms = 47681
    ref_numresidues = 11302


class TestTPRParser(_TestTopology, RefTPR):
    pass


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
        angs = self.u.angles.atomgroup_intersection(ag, strict=True)

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
        self.u.angles = guess_angles(self.u.bonds)

        assert_equal(len(self.u.angles), 6123)

    def test_guess_dihedrals_set(self):
        self.u.dihedrals = guess_dihedrals(self.u.angles)

        assert_equal(len(self.u.dihedrals), 8921)

    def test_guess_impropers_set(self):
        self.u.impropers = guess_improper_dihedrals(self.u.angles)

        assert_equal(len(self.u.impropers), 10314)


class RefLammpsData(object):
    topology = LAMMPSdata
    parser = MDAnalysis.topology.LAMMPSParser.DATAParser
    ref_n_atoms = 18360
    ref_numresidues = 24


class TestLammpsData(_TestTopology, RefLammpsData):
    """Tests the reading of lammps .data topology files.

    The reading of coords and velocities is done separately in test_coordinates
    """

    def test_charge(self):
        # No charges were supplied, should default to 0.0
        assert_equal(self.universe.atoms[0].charge, 0.0)

    def test_resid(self):
        assert_equal(len(self.universe.residues[0]), 765)

    # Testing _psf prevent building TGs
    # test length and random item from within
    def test_bonds(self):
        assert_equal(len(self.universe._topology['bonds']), 18336)
        assert_equal((5684, 5685) in self.universe._topology['bonds'], True)

    def test_angles(self):
        assert_equal(len(self.universe._topology['angles']), 29904)
        assert_equal((7575, 7578, 7579) in self.universe._topology['angles'], True)

    def test_dihedrals(self):
        assert_equal(len(self.universe._topology['dihedrals']), 5712)
        assert_equal((3210, 3212, 3215, 3218) in self.universe._topology['dihedrals'],
                     True)

    def test_masses(self):
        assert_equal(self.universe.atoms[0].mass, 0.012)


class RefXYZ(object):
    topology = XYZ_mini
    parser = MDAnalysis.topology.XYZParser.XYZParser
    ref_n_atoms = 3
    ref_numresidues = 1


class TestXYZTopology(RefXYZ, _TestTopology):
    def test_segments(self):
        assert_equal(len(self.universe.segments), 1)


class RefGMSsym(object):
    topology = GMS_SYMOPT
    parser = MDAnalysis.topology.GMSParser.GMSParser
    ref_n_atoms = 4
    ref_numresidues = 1


class TestGMS_withSymmetry(_TestTopology, RefGMSsym):
    """Testing GAMESS output file format"""


class RefGMSasym(object):
    topology = GMS_ASYMSURF
    parser = MDAnalysis.topology.GMSParser.GMSParser
    ref_n_atoms = 6
    ref_numresidues = 1


class TestGMS_noSymmetry(_TestTopology, RefGMSasym):
    """Testing GAMESS output file format"""


class _DLPolyParser(object):
    """Test of real data"""
    def tearDown(self):
        del self.p
        del self.f

    def test_usage(self):
        with self.p(self.f) as parser:
            struc = parser.parse()

        assert_equal('atoms' in struc, True)
        assert_equal(len(struc['atoms']), 216)

    def test_names(self):
        with self.p(self.f) as parser:
            struc = parser.parse()

        atoms = struc['atoms']

        assert_equal(atoms[0].name, 'K+')
        assert_equal(atoms[4].name, 'Cl-')


class TestDLPolyConfigParser(_DLPolyParser):
    def setUp(self):
        self.p = MDAnalysis.topology.DLPolyParser.ConfigParser
        self.f = DLP_CONFIG


class TestDLPolyHistoryParser(_DLPolyParser):
    def setUp(self):
        self.p = MDAnalysis.topology.DLPolyParser.HistoryParser
        self.f = DLP_HISTORY


# Artificial DL_Poly data for testing limits
class _DLPoly(object):
    def tearDown(self):
        del self.p
        del self.f

    def test_usage(self):
        with self.p(self.f) as parser:
            struc = parser.parse()

        assert_equal('atoms' in struc, True)
        assert_equal(len(struc['atoms']), 3)

    def test_names(self):
        with self.p(self.f) as parser:
            struc = parser.parse()

        atoms = struc['atoms']

        assert_equal(atoms[0].name, 'C')
        assert_equal(atoms[1].name, 'B')
        assert_equal(atoms[2].name, 'A')


class TestDLPolyConfigOrder(_DLPoly):
    def setUp(self):
        self.p = MDAnalysis.topology.DLPolyParser.ConfigParser
        self.f = DLP_CONFIG_order


class TestDLPolyConfigMinimal(_DLPoly):
    def setUp(self):
        self.p = MDAnalysis.topology.DLPolyParser.ConfigParser
        self.f = DLP_CONFIG_minimal


class TestDLPolyHistoryOrder(_DLPoly):
    def setUp(self):
        self.p = MDAnalysis.topology.DLPolyParser.HistoryParser
        self.f = DLP_HISTORY_order


class TestDLPolyHistoryMinimal(_DLPoly):
    def setUp(self):
        self.p = MDAnalysis.topology.DLPolyParser.HistoryParser
        self.f = DLP_HISTORY_minimal
