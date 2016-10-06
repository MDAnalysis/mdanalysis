from six.moves import zip
import os
from numpy.testing import (
    assert_,
    assert_array_equal,
)

import MDAnalysis as mda

from MDAnalysisTests.datafiles import CRD
from MDAnalysisTests import tempdir

class TestCRDWriter(object):
    def setUp(self):
        self.u = mda.Universe(CRD)
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/out.crd'

    def tearDown(self):
        del self.u
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.tmpdir
        del self.outfile

    def test_write_atoms(self):
        # Test that written file when read gives same coordinates
        self.u.atoms.write(self.outfile)

        u2 = mda.Universe(self.outfile)

        assert_array_equal(self.u.atoms.positions,
                           u2.atoms.positions)

    def test_roundtrip(self):
        # Write out a copy of the Universe, and compare this against the original
        # This is more rigorous than simply checking the coordinates as it checks
        # all formatting
        self.u.atoms.write(self.outfile)
        
        def CRD_iter(fn):
            with open(fn, 'r') as inf:
                for line in inf:
                    if not line.startswith('*'):
                        yield line

        for ref, other in zip(CRD_iter(CRD), CRD_iter(self.outfile)):
            assert_(ref == other)

    def test_write_EXT(self):
        # TODO: Write tests that use EXT output format
        # Must have *lots* of atoms, maybe fake the system
        # to make tests faster
        pass
