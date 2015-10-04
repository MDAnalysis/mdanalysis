import MDAnalysis as mda
import numpy as np
import os

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_almost_equal, dec,
                           assert_array_almost_equal, assert_raises)
from unittest import TestCase
import tempdir

from MDAnalysisTests.datafiles import (GRO)
from MDAnalysisTests.coordinates.reference import RefAdK


class TestGROReader(TestCase, RefAdK):
    def setUp(self):
        self.universe = mda.Universe(GRO)
        self.ts = self.universe.trajectory.ts
        # lower prec in gro!! (3 decimals nm -> 2 decimals in Angstroem)
        self.prec = 2

    def tearDown(self):
        del self.universe
        del self.ts

    def test_flag_convert_lengths(self):
        assert_equal(mda.core.flags['convert_lengths'], True,
                     "MDAnalysis.core.flags['convert_lengths'] should be True "
                     "by default")

    def test_load_gro(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from small GRO")
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1,
                     "wrong number of frames")

    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0,
                     "wrong time of the frame")

    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0,
                     "wrong frame number (should be 0 for 0-based ts.frame)")

    def test_dt(self):
        """testing that accessing universe.trajectory.dt gives
        default of 1.0"""
        assert_equal(self.universe.trajectory.dt, 1.0)

    def test_coordinates(self):
        A10CA = self.universe.SYSTEM.CA[10]
        assert_almost_equal(A10CA.pos,
                            self.ref_coordinates['A10CA'],
                            self.prec,
                            err_msg="wrong coordinates for A10:CA")

    def test_distances(self):
        # NOTe that the prec is only 1 decimal: subtracting two low precision
        #      coordinates low prec: 9.3455122920041109; high prec (from pdb):
        #      9.3513174
        NTERM = self.universe.SYSTEM.N[0]
        CTERM = self.universe.SYSTEM.C[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec - 1,
                            err_msg="distance between M1:N and G214:C")

    def test_selection(self):
        na = self.universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size,
                     "Atom selection of last atoms in file")

    def test_unitcell(self):
        assert_array_almost_equal(
            self.ts.dimensions,
            self.ref_unitcell,
            self.prec,
            err_msg="unit cell dimensions (rhombic dodecahedron)")

    def test_volume(self):
        # test_volume: reduce precision for Gromacs comparison to 0 decimals
        # (A**3 <--> nm**3!)
        assert_almost_equal(
            self.ts.volume,
            self.ref_volume,
            0,
            err_msg="wrong volume for unitcell (rhombic dodecahedron)")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))


class TestGROReaderNoConversion(TestCase, RefAdK):
    def setUp(self):
        self.universe = mda.Universe(GRO, convert_units=False)
        self.ts = self.universe.trajectory.ts
        self.prec = 3

    def tearDown(self):
        del self.universe
        del self.ts

    def test_coordinates(self):
        # note: these are the native coordinates in nm; for the test to succeed
        # we loaded with convert_units=False
        A10CA = self.universe.SYSTEM.CA[10]
        # coordinates in nm
        assert_almost_equal(A10CA.pos, RefAdK.ref_coordinates['A10CA'] / 10.0,
                            self.prec, err_msg="wrong native coordinates "
                            "(in nm) for A10:CA")

    def test_distances(self):
        # 3 decimals on nm in gro but we compare to the distance
        # computed from the pdb file, so the effective precision is 2 again.
        # (Otherwise the distance test fails:
        #  Arrays are not almost equal distance between M1:N and G214:C
        #    ACTUAL: 0.93455122920041123
        #    DESIRED: 0.93513173999999988
        NTERM = self.universe.SYSTEM.N[0]
        CTERM = self.universe.SYSTEM.C[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        # coordinates in nm
        assert_almost_equal(d, RefAdK.ref_distances['endtoend'] / 10.0,
                            self.prec - 1, err_msg="distance between M1:N "
                            "and G214:C")

    def test_unitcell(self):
        # lengths in A : convert to nm
        assert_array_almost_equal(
            self.ts.dimensions[:3],
            self.ref_unitcell[:3] / 10.0,
            self.prec,
            err_msg="unit cell A,B,C (rhombic dodecahedron)")
        # angles should not have changed
        assert_array_almost_equal(
            self.ts.dimensions[3:],
            self.ref_unitcell[3:],
            self.prec,
            err_msg="unit cell alpha,beta,gamma (rhombic dodecahedron)")

    def test_volume(self):
        # ref lengths in A (which was originally converted from nm)
        assert_almost_equal(
            self.ts.volume,
            self.ref_volume / 1000.,
            3,
            err_msg="wrong volume for unitcell (rhombic dodecahedron)")


class TestGROWriter(TestCase, tempdir.TempDir):
    def setUp(self):
        self.universe = mda.Universe(GRO)
        self.prec = 2  # 3 decimals in file in nm but MDAnalysis is in A
        ext = ".gro"
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/gro-writer' + ext
        self.outfile2 = self.tmpdir.name + '/gro-writer2' + ext

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        try:
            os.unlink(self.outfile2)
        except OSError:
            pass
        del self.universe
        del self.tmpdir

    @dec.slow
    def test_writer(self):
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_almost_equal(u.atoms.coordinates(),
                            self.universe.atoms.coordinates(), self.prec,
                            err_msg="Writing GRO file with GROWriter does "
                            "not reproduce original coordinates")

    @dec.slow
    def test_timestep_not_modified_by_writer(self):
        ts = self.universe.trajectory.ts
        x = ts._pos.copy()
        self.universe.atoms.write(self.outfile)
        assert_equal(ts._pos,
                     x,
                     err_msg="Positions in Timestep were modified by writer.")

    @dec.slow
    @attr('issue')
    def test_check_coordinate_limits_min(self):
        """Test that illegal GRO coordinates (x <= -999.9995 nm) are caught
        with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(GRO)
        u.atoms[2000].pos[1] = -999.9995 * 10  # nm -> A
        assert_raises(ValueError, u.atoms.write, self.outfile2)
        del u

    @dec.slow
    @attr('issue')
    def test_check_coordinate_limits_max(self):
        """Test that illegal GRO coordinates (x > 9999.9995 nm) are caught
        with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(GRO)
        # nm -> A  ; [ob] 9999.9996 not caught
        u.atoms[1000].pos[1] = 9999.9999 * 10
        assert_raises(ValueError, u.atoms.write, self.outfile2)
        del u

    @dec.slow
    def test_check_coordinate_limits_max_noconversion(self):
        """Test that illegal GRO coordinates (x > 9999.9995 nm) also
        raises exception for convert_units=False"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(GRO, convert_units=False)
        u.atoms[1000].pos[1] = 9999.9999
        assert_raises(ValueError, u.atoms.write, self.outfile2,
                      convert_units=False)
        del u
