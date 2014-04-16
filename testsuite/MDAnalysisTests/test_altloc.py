from MDAnalysis import Universe
import tempfile, os
from numpy.testing import *
from .datafiles import altloc

class TestAltloc(TestCase):
    def setUp(self):
        fd, self.outfile = tempfile.mkstemp(suffix=".pdb")  # output is always same as input (=DCD)
        os.close(fd)
    def tearDown(self):
        try:
            os.unlink(self.outfile)  
        except OSError:
            pass
    def test_atomgroups(self):
        u = Universe(altloc)
        segidB0 = len(u.selectAtoms("segid B and (not altloc B)"))
        segidB1 = len(u.selectAtoms("segid B and (not altloc A)"))
        assert_equal(segidB0, segidB1)
        altlocB0 = len(u.selectAtoms("segid B and (altloc A)"))
        altlocB1 = len(u.selectAtoms("segid B and (altloc B)"))
        assert_equal(altlocB0, altlocB1)
        sum = len(u.selectAtoms("segid B"))
        assert_equal(sum, segidB0 + altlocB0)
    
    def test_bonds(self): 
        u = Universe(altloc, bonds=True)
        bonds0 = u.selectAtoms("segid B and (altloc A)")[0].bonds
        bonds1 = u.selectAtoms("segid B and (altloc B)")[0].bonds
        assert_equal(len(bonds0), len(bonds1))
        
    def test_write_read(self):
        u = Universe(altloc)
        u.selectAtoms("all").write(self.outfile)
        u2 = Universe(self.outfile)
        assert_equal( len(u.atoms), len(u2.atoms))