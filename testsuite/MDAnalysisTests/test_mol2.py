from numpy.testing import *
import os, tempfile
from .datafiles import mol2_molecules, mol2_molecule, mol2_broken_molecule
from MDAnalysis import Universe

class TestMol2(TestCase):
    def setUp(self):
        fd, self.outfile = tempfile.mkstemp(suffix=".mol2")
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

    def test_read(self):
        u = Universe(mol2_molecules)
        assert_equal(len(u.atoms), 49)
        assert_equal(u.trajectory.numframes, 200)
        
        u.trajectory[199]
        assert_array_almost_equal(u.atoms.positions[0], [1.7240, 11.2730, 14.1200])
        
    def test_write(self):
        ref = Universe(mol2_molecules)
        ref.atoms.write(self.outfile)
        u = Universe(self.outfile)
        assert_equal(len(u.atoms), len(ref.atoms))
        assert_equal(len(u.trajectory), len(ref.trajectory))
        assert_array_equal(u.atoms.positions, ref.atoms.positions)
        u.trajectory[199]; ref.trajectory[199]
        assert_array_equal(u.atoms.positions, ref.atoms.positions)
        
    def test_write_selection(self):
        ref = Universe(mol2_molecule)
        gr0 = ref.selectAtoms("name C*")
        gr0.write(self.outfile)
        u = Universe(self.outfile)
        gr1 = u.selectAtoms("name C*")
        assert_equal(len(gr0), len(gr1))
      
    def test_broken_molecule(self):
        with self.assertRaises(Exception) as context:
            u = Universe(mol2_broken_molecule)
        self.assertEqual(context.exception.message, "The mol2 block (BrokenMolecule.mol2:0) has no atoms")