from six.moves import zip
import os
from numpy.testing import (
    assert_,
    assert_equal,
    assert_array_equal,
    assert_warns,
)

import MDAnalysis as mda

from MDAnalysisTests.datafiles import CRD
from MDAnalysisTests import tempdir
from MDAnalysisTests.core.groupbase import make_Universe

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

class TestCRDWriterMissingAttrs(object):
    # All required attributes with the default value
    req_attrs = {'resnames': 'UNK',
                 'resids': 1,
                 'names': 'X',
                 'tempfactors': 0.0,
                 }

    def _check_warns(self, missing_attr):
        attrs = list(self.req_attrs.keys())
        attrs.remove(missing_attr)
        u = make_Universe(attrs, trajectory=True)

        tmpdir = tempdir.TempDir()
        outfile = tmpdir.name + '/out.crd'
        assert_warns(UserWarning,
                     u.atoms.write, outfile)

    def _check_write(self, missing_attr):
        attrs = list(self.req_attrs.keys())
        attrs.remove(missing_attr)
        u = make_Universe(attrs, trajectory=True)

        tmpdir = tempdir.TempDir()
        outfile = tmpdir.name + '/out.crd'
        u.atoms.write(outfile)
        u2 = mda.Universe(outfile)

        # Check all other attrs aren't disturbed
        for attr in attrs:
            assert_equal(getattr(u.atoms, attr),
                         getattr(u2.atoms, attr))
        # Check missing attr is as expected
        assert_equal(getattr(u2.atoms, missing_attr),
                     self.req_attrs[missing_attr])

    def test_crdwriter(self):
        for attr in self.req_attrs:
            yield self._check_warns, attr
            yield self._check_write, attr
            
