import MDAnalysis
from MDAnalysis.tests.datafiles import PSF,DCD

from numpy.testing import *
from nose.plugins.attrib import attr

import tempfile

class _WriteAtoms(TestCase):
    """Set up the standard AdK system in implicit solvent."""
    ext = None   # override to test various output writers
    precision = 3

    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        suffix = '.' + self.ext
        fd, self.outfile = tempfile.mkstemp(suffix=suffix)
        
    def universe_from_tmp(self):
        return MDAnalysis.Universe(self.outfile)

    def test_write_atoms(self):
        """Test writing all atoms to external file."""
        self.universe.atoms.write(self.outfile)
        u2 = self.universe_from_tmp()
        assert_array_almost_equal(self.universe.atoms.coordinates(), u2.atoms.coordinates(), self.precision,
                                  err_msg="atom coordinate mismatch between original and %s file" % self.ext)

    def test_write_selection(self):
        """Test writing CA atoms to external file."""
        CA = self.universe.selectAtoms('name CA')
        CA.write(self.outfile)
        u2 = self.universe_from_tmp()
        CA2 = u2.selectAtoms('name CA')

        assert_equal(len(u2.atoms), len(CA), "written CA selection does not match original selection")
        assert_almost_equal(CA2.coordinates(), CA.coordinates(), self.precision,
                            err_msg="CA coordinates do not agree with original")

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass

class TestWritePDB(_WriteAtoms):
    ext = "pdb"
    precision = 3

import MDAnalysis.coordinates
class TestWriteCRD(_WriteAtoms):
    ext = "crd"
    precision = 2

    # replicated methods here so that I can decorate as known failures;
    # remove once Issue 40 is resolved.
    @dec.knownfailureif(not 'CRD' in MDAnalysis.coordinates._trajectory_readers, 
                        "CRD reader is not implemented yet (see Issue 40)")
    def test_write_atoms(self):
        super(TestWriteCRD, self).test_write_atoms()

    @dec.knownfailureif(not 'CRD' in MDAnalysis.coordinates._trajectory_readers, 
                        "CRD reader is not implemented yet (see Issue 40)")
    def test_write_selection(self):
        super(TestWriteCRD, self).test_write_selection()

class TestWriteGRO(_WriteAtoms):
    ext = "gro"
    precision = 2


