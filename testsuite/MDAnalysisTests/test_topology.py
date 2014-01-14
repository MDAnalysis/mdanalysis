# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

import MDAnalysis
from MDAnalysis.topology.core import guess_atom_type, guess_atom_element, get_atom_mass
from MDAnalysis.tests.datafiles import PRMpbc, PRM12, PSF, PSF_NAMD, PSF_nosegid, DMS, PDB_small

from numpy.testing import *
from nose.plugins.attrib import attr

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
    testnames = ['CAL','CA2+', 'C0']

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

    def test_parser(self):
        assert_equal(self.universe.atoms.numberOfAtoms(),
                     self.ref_numatoms,
                     "wrong number of atoms in topology")
        assert_equal(self.universe.atoms.numberOfResidues(),
                     self.ref_numresidues,
                     "wrong number of residues in topology")

    def test_atom_number(self):
        assert_equal(self.universe.atoms[0].number, 0,
                     "first atom should have Atom.number 0")
        assert_equal(self.universe.atoms[-1].number, self.ref_numatoms - 1,
                     "last atom has wrong Atom.number")

class RefAdKSmall(object):
    """Mixin class to provide comparison numbers.

    Based on small PDB with AdK (:data:`PDB_small`).
    """
    topology = PSF
    ref_numatoms = 3341
    ref_numresidues = 214

class TestPSF_CHARMM_STANDARD(_TestTopology, RefAdKSmall):
    """Testing CHARMM standard PSF file format"""

class RefNAMD_CGENFF(object):
    """Testfiles provided by JiyongPark77.

    NAMD/VMD XPLOR-style PSF file (using CGENFF residues/atoms).

    http://code.google.com/p/mdanalysis/issues/detail?id=107
    """
    topology = PSF_NAMD
    ref_numatoms = 130
    ref_numresidues = 6

class TestPSF_NAMD_CGENFF(_TestTopology, RefNAMD_CGENFF):
    """Testing NAMD PSF file (with CGENFF atom types, Issue 107)"""

class TestPSF_Issue121(TestCase):
    @attr('issue')
    def test_nosegid(self):
        try:
            u = MDAnalysis.Universe(PSF_nosegid)
        except IndexError:
            raise AssertionError("Issue 121 not fixed: cannot load PSF with empty SEGID")
        assert_equal(u.atoms.numberOfAtoms(), 98)
        assert_equal(u.atoms.segids(), ["SYSTEM"])

class TestPSF_bonds(TestCase):
    """Tests reading of bonds angles and torsions in psf files"""

    def setUp(self):
        topology = PSF
        self.universe = MDAnalysis.Universe(topology)

    def tearDown(self):
        del self.universe

    def test_bonds_counts(self):
        assert_equal(len(self.universe._psf['_bonds']), 3365)
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
        assert_equal(all([a1 in b for b in a1.bonds]), True) # check all bonds have this atom
        assert_equal(any([a2 in b for b in a1.bonds]), True) # then check certain atoms are present
        assert_equal(any([a3 in b for b in a1.bonds]), True)
        assert_equal(any([a4 in b for b in a1.bonds]), True)
        assert_equal(any([a5 in b for b in a1.bonds]), True)
        assert_equal(any([a42 in b for b in a1.bonds]), False) # and check everything isn't True

    def test_angles_counts(self):
        assert_equal(len(self.universe._psf['_angles']), 6123)
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
        assert_equal(any([a2 in b and a6 in b for b in a1.angles]), False) #both a2 and a6 feature, but never simultaneously 

    def test_torsions_counts(self):
        assert_equal(len(self.universe._psf['_dihe']), 8921)
        assert_equal(len(self.universe.atoms[0].torsions), 14)

    def test_torsions_identity(self):
        a1 = self.universe.atoms[0]
        a2 = self.universe.atoms[1]
        a3 = self.universe.atoms[2]
        a5 = self.universe.atoms[4]
        a6 = self.universe.atoms[5]
        a7 = self.universe.atoms[6]
        a42 = self.universe.atoms[41]  
        assert_equal(all([a1 in b for b in a1.torsions]), True)
        assert_equal(any([a2 in b and a5 in b and a6 in b for b in a1.torsions]), True)
        assert_equal(any([a2 in b and a5 in b and a7 in b for b in a1.torsions]), True)
        assert_equal(any([a42 in b for b in a1.torsions]), False)
        assert_equal(any([a2 in b and a3 in b and a6 in b for b in a1.torsions]), False)

class TestPSF_TopologyGroup(TestCase):
    """Tests TopologyDict and TopologyGroup classes with psf input"""
    def setUp(self):
        topology = PSF
        self.universe = MDAnalysis.Universe(topology)
        self.res1 = self.universe.residues[0]
        self.res2 = self.universe.residues[1]

    def tearDown(self):
        del self.universe

    def test_bonds(self):
        assert_equal(len(self.universe.atoms.bondDict), 57)
        assert_equal(self.res1.numberOfBondTypes(), 12)

    def test_angles(self):
        assert_equal(len(self.universe.atoms.angleDict), 130)

    def test_torsions(self):
        assert_equal(len(self.universe.atoms.torsionDict), 220)

    def test_create_TopologyGroup(self):
        res1_tg = self.res1.bondDict['23', '3'] # make a tg
        assert_equal(len(res1_tg), 4)

        testbond = self.universe.atoms[7].bonds[0]
        assert_equal(testbond in res1_tg, True) # check a known bond is present

    def test_add_TopologyGroups(self):
        res1_tg = self.res1.bondDict['23', '3']
        res2_tg = self.res2.selectBonds(('23', '3'))
        assert_equal(len(res2_tg), 6)

        combined_tg = res1_tg + res2_tg # add tgs together
        assert_equal(len(combined_tg), 10)

        big_tg = self.universe.atoms.bondDict['23', '3']
        assert_equal(len(big_tg), 494)

        big_tg += combined_tg # try and add some already included bonds
        assert_equal(len(big_tg), 494) # check len doesn't change

# AMBER

class RefCappedAla(object):
    """Mixin class to provide comparison numbers.

    Capped Ala in water
    """
    topology = PRMpbc
    ref_numatoms = 5071
    ref_numresidues = 1686
    ref_proteinatoms = 22

class RefAMBER12(object):
    """Fixture data for testing AMBER12 reading (Issue 100)"""
    topology = PRM12
    ref_numatoms = 8923
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
    ref_numatoms = 3341
    ref_numresidues = 214

class TestPDB(_TestTopology, RefPDB):
    """Testing PDB topology parsing (PrimitivePDB)"""

# DESRES

class TestDMSReader(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(DMS)
        self.ts = self.universe.trajectory.ts

    def tearDown(self):
        del self.universe
        del self.ts

    def test_number_of_bonds(self):
        # Desired value taken from VMD
        #      Info)    Atoms: 3341
        assert_equal(len(self.universe.bonds),3365)

    def test_bond_order(self):
        pass

    def test_segid(self):
        segid = set([a.segid for a in self.universe.atoms])
        assert_equal(segid, set(("4AKE",) ))

    def test_atomsels(self):
        # Desired value taken from VMD atomsel
        s0 = self.universe.selectAtoms("name CA")
        assert_equal(len(s0), 214)

        s1 = self.universe.selectAtoms("resid 33")
        assert_equal(len(s1), 12)

        s2 = self.universe.selectAtoms("segid 4AKE")
        assert_equal(len(s2), 3341)

        s3 = self.universe.selectAtoms("resname ALA")
        assert_equal(len(s3), 190)

    def test_atom_number(self):
        assert_equal(self.universe.atoms[0].number, 0,
                     "first atom should have Atom.number 0")
        assert_equal(self.universe.atoms[-1].number, 3341 - 1,
                     "last atom has wrong Atom.number")

# GROMACS TPR
# see test_tprparser
