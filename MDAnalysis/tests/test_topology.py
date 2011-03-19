import MDAnalysis
from MDAnalysis.topology.core import guess_atom_type

from numpy.testing import *

def check_atom_type(atype, aname):
    assert_equal(guess_atom_type(aname), atype)

class _TestGuessAtomType(object):
    atype = None
    testnames = []
    def test_guess_atom_type(self):        
        for aname in self.testnames:
            yield check_atom_type, self.atype, aname

class TestHydrogen(_TestGuessAtomType):
    atype = 'H'
    testnames = ['H', 'HZ', '1HZ', '2HW']

class TestCarbon(_TestGuessAtomType):
    atype = 'C'
    testnames = ['C', 'CA']

class TestSodium(_TestGuessAtomType):
    atype = 'NA'
    testnames = ['NA', 'NA+', 'SOD', 'QN']

class TestPotassium(_TestGuessAtomType):
    atype = 'K'
    testnames = ['K', 'K+', 'POT', 'QK']

class TestChloride(_TestGuessAtomType):
    atype = 'CL'
    testnames = ['CL', 'CL-', 'CLA', 'CLAL']

class TestNitrogen(_TestGuessAtomType):
    atype = 'N'
    testnames = ['N', 'NZ']

class TestPhosphate(_TestGuessAtomType):
    atype = 'P'
    testnames = ['P', 'PO']

class TestSulphur(_TestGuessAtomType):
    atype = 'S'
    testnames = ['S', 'SG']

class TestOxygen(_TestGuessAtomType):
    atype = 'O'
    testnames = ['O', 'OA', 'OXT', '1OG', 'OW']

class TestCalcium(_TestGuessAtomType):
    atype = 'CA'
    testnames = ['CAL','CA2+', 'C0'] 

# add more...
