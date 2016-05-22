import MDAnalysis as mda
import os

from numpy.testing import (assert_almost_equal, assert_equal)
from unittest import TestCase

from MDAnalysisTests.coordinates.reference import (RefAdKSmall)
from MDAnalysisTests.coordinates.base import _SingleFrameReader
from MDAnalysisTests.datafiles import (PQR)
from MDAnalysisTests import tempdir


class TestPQRReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(PQR)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3

    def test_total_charge(self):
        assert_almost_equal(
            self.universe.atoms.total_charge(), self.ref_charmm_totalcharge, 3,
            "Total charge (in CHARMM) does not match expected value.")

    def test_hydrogenCharges(self):
        assert_almost_equal(self.universe.atoms.H.charges,
                            self.ref_charmm_Hcharges, 3,
                            "Charges for H atoms do not match.")

    # Note that the whole system gets the sysID 'SYSTEM' for the PQR file (when
    # read with a PSF it is 's4AKE')
    def test_ArgCACharges(self):
        assert_almost_equal(
            self.universe.SYSTEM.ARG.CA.charges, self.ref_charmm_ArgCAcharges,
            3, "Charges for CA atoms in Arg residues do not match.")

    def test_ProNCharges(self):
        assert_almost_equal(
            self.universe.SYSTEM.PRO.N.charges, self.ref_charmm_ProNcharges, 3,
            "Charges for N atoms in Pro residues do not match.")


class TestPQRWriter(TestCase, RefAdKSmall):
    def setUp(self):
        self.universe = mda.Universe(PQR)
        self.prec = 3
        ext = ".pqr"
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/pqr-writer-test' + ext

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe
        del self.tmpdir

    def test_writer_noChainID(self):
        assert_equal(self.universe.segments.segids[0], 'SYSTEM')
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_equal(u.segments.segids[0], 'SYSTEM')
        assert_almost_equal(u.atoms.positions,
                            self.universe.atoms.positions, self.prec,
                            err_msg="Writing PQR file with PQRWriter does "
                            "not reproduce original coordinates")
        assert_almost_equal(u.atoms.charges, self.universe.atoms.charges,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original charges")
        assert_almost_equal(u.atoms.radii, self.universe.atoms.radii,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original radii")

    def test_write_withChainID(self):
        self.universe.atoms.set_segids('A')
        assert_equal(self.universe.segments.segids[0], 'A')  # sanity check
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_equal(u.segments.segids[0], 'A')
        assert_almost_equal(u.atoms.positions,
                            self.universe.atoms.positions, self.prec,
                            err_msg="Writing PQR file with PQRWriter does "
                            "not reproduce original coordinates")
        assert_almost_equal(u.atoms.charges, self.universe.atoms.charges,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original charges")
        assert_almost_equal(u.atoms.radii, self.universe.atoms.radii,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original radii")

    def test_timestep_not_modified_by_writer(self):
        ts = self.universe.trajectory.ts
        x = ts._pos.copy()
        self.universe.atoms.write(self.outfile)
        assert_equal(ts._pos,
                     x,
                     err_msg="Positions in Timestep were modified by writer.")

    def test_total_charge(self):
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_almost_equal(
            u.atoms.total_charge(), self.ref_charmm_totalcharge, 3,
            "Total charge (in CHARMM) does not match expected value.")
