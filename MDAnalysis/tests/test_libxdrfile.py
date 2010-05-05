import MDAnalysis
from numpy.testing import *
from MDAnalysis.tests.datafiles import XTC,TRR

import MDAnalysis.coordinates.xdrfile.libxdrfile as xdr

class TestLib(TestCase):
    def test_constants(self):
        assert_equal(xdr.DIM, 3, "xdr library not compiled for DIM=3 ?!?")

    def test_xdropen(self):
        XDR = xdr.xdrfile_open(XTC, 'r')
        assert_(XDR != None, "Failed top open xtc file")
        rc = xdr.xdrfile_close(XDR)
        assert_equal(rc, 0, "Failed to close xtc file")  # this can segfault 


class TestXTC(TestCase):
    def test_numatoms(self):
        natoms = xdr.read_xtc_natoms(XTC)
        assert_equal(natoms, 47681, "Number of atoms in XTC frame")

    def test_numframes(self):
        numframes = xdr.read_xtc_numframes(XTC)
        assert_equal(numframes, 10, "Number of frames in XTC trajectory")


class TestTRR(TestCase):
    def test_numatoms(self):
        natoms = xdr.read_trr_natoms(TRR)
        assert_equal(natoms, 47681, "Number of atoms in TRR frame")

    def test_numframes(self):
        numframes = xdr.read_trr_numframes(TRR)
        assert_equal(numframes, 10, "Number of frames in TRR trajectory")

