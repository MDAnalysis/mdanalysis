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
from MDAnalysis.tests.datafiles import PRMpbc, PRM12, PSF, PSF_NAMD, PSF_nosegid

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
    def test_parser(self):
        U = MDAnalysis.Universe(self.topology)
        assert_equal(U.atoms.numberOfAtoms(), self.ref_numatoms, "wrong number of atoms in topology")
        assert_equal(U.atoms.numberOfResidues(), self.ref_numresidues, "wrong number of residues in topology")

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

