import MDAnalysis
from MDAnalysis.topology.core import guess_atom_type, get_atom_mass

from numpy.testing import *

def check_atom_type(atype, aname):
    assert_equal(guess_atom_type(aname), atype)

class _TestGuessAtomType(object):
    atype = None
    testnames = []
    mass = None
    def test_guess_atom_type(self):        
        for aname in self.testnames:
            yield check_atom_type, self.atype, aname

    def test_guess_atom_mass(self):        
        assert_equal(get_atom_mass(self.atype), self.mass)


class TestHydrogen(_TestGuessAtomType):
    atype = 'H'
    mass = 1.008
    testnames = ['H', 'HZ', '1HZ', '2HW']

class TestCarbon(_TestGuessAtomType):
    atype = 'C'
    mass = 12.0110
    testnames = ['C', 'CA']

class TestSodium(_TestGuessAtomType):
    atype = 'NA'
    mass = 22.989770
    testnames = ['NA', 'NA+', 'SOD', 'QN']

class TestPotassium(_TestGuessAtomType):
    atype = 'K'
    mass = 39.102
    testnames = ['K', 'K+', 'POT', 'QK']

class TestChloride(_TestGuessAtomType):
    atype = 'CL'
    mass = 35.450
    testnames = ['CL', 'CL-', 'CLA', 'CLAL']

class TestNitrogen(_TestGuessAtomType):
    atype = 'N'
    mass = 14.007
    testnames = ['N', 'NZ']

class TestPhosphorous(_TestGuessAtomType):
    atype = 'P'
    mass = 30.974000
    testnames = ['P', 'PL', 'PO']

class TestSulfur(_TestGuessAtomType):
    atype = 'S'
    mass = 32.06000
    testnames = ['S', 'SG']

class TestOxygen(_TestGuessAtomType):
    atype = 'O'
    mass = 15.99900
    testnames = ['O', 'OA', 'OXT', '1OG', 'OW']

class TestCalcium(_TestGuessAtomType):
    atype = 'CA'
    mass = 40.080000
    testnames = ['CAL','CA2+', 'C0'] 

class TestMagnesium(_TestGuessAtomType):
    atype = 'MG'
    mass = 24.305000
    testnames = ['MG', 'MG2+']

# add more...
