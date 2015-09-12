# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import MDAnalysis
import MDAnalysis as mda
import MDAnalysis.coordinates
import MDAnalysis.coordinates.core
from MDAnalysis import NoDataError

import numpy as np
from numpy.testing import *
from nose.plugins.attrib import attr

from MDAnalysisTests.datafiles import (
    PSF, DCD, DCD_empty, PDB_small, XPDB_small, PDB_closed, PDB_multiframe,
    PDB_full, PDB, CRD, XTC, TRR, GRO, DMS, CONECT, INC_PDB, XYZ, XYZ_bz2,
    XYZ_psf, PRM, TRJ, TRJ_bz2, PRMpbc, TRJpbc_bz2, PRMncdf, NCDF, PQR,
    PDB_sub_dry, TRR_sub_sol, PDB_sub_sol, TRZ, TRZ_psf, LAMMPSdata,
    LAMMPSdata_mini, PSF_TRICLINIC, DCD_TRICLINIC, PSF_NAMD_TRICLINIC,
    DCD_NAMD_TRICLINIC, GMS_ASYMOPT, GMS_SYMOPT, GMS_ASYMSURF, XYZ_mini,
    PFncdf_Top, PFncdf_Trj, INPCRD, XYZ_five, two_water_gro, DLP_CONFIG,
    DLP_CONFIG_order, DLP_CONFIG_minimal, DLP_HISTORY, DLP_HISTORY_order,
    DLP_HISTORY_minimal)

from MDAnalysisTests.plugins.knownfailure import knownfailure

import os
import errno
import tempfile
import itertools


def atom_distance(a, b):
    """Calculate the distance between two atoms a and b."""
    r = a.position - b.position
    return np.sqrt(np.sum(r ** 2))


class RefAdKSmall(object):
    """Mixin class to provide comparison numbers.

    Based on small PDB with AdK (:data:`PDB_small`).

    .. Note::

       All distances must be in ANGSTROEM as this is the MDAnalysis
       default unit. All readers must return Angstroem by default.
    """
    filename = PDB_small
    ref_coordinates = {
        'A10CA': np.array([-1.198, 7.937, 22.654]),  # G11:CA, copied frm adk_open.pdb
    }
    ref_distances = {'endtoend': 11.016959}
    ref_E151HA2_index = 2314
    ref_n_atoms = 3341
    ref_charmm_totalcharge = -4.0
    ref_charmm_Hcharges = [0.33] + 203 * [0.31]
    ref_charmm_ArgCAcharges = 13 * [0.07]
    ref_charmm_ProNcharges = 10 * [-0.29]
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.],
                            dtype=np.float32)
    ref_volume = 0.0


class RefAdK(object):
    """Mixin class to provide comparison numbers.

    Based on PDB/GRO with AdK in water + Na+ (:data:`PDB`).

    .. Note::

       All distances must be in ANGSTROEM as this is the MDAnalysis
       default unit. All readers must return Angstroem by default.
    """
    filename = PDB
    ref_coordinates = {
        'A10CA': np.array([62.97600174, 62.08800125, 20.2329998]),  # Angstroem as MDAnalysis unit!!
    }
    ref_distances = {'endtoend': 9.3513174}
    ref_E151HA2_index = 2314
    ref_n_atoms = 47681
    ref_Na_sel_size = 4
    # CRYST1 80.017   80.017   80.017  60.00  60.00  90.00
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.], dtype=np.float32)
    ref_volume = 362270.0  # computed with Gromacs


class Ref2r9r(object):
    """Mixin class to provide comparison numbers.

    Based on S6 helices of chimeric Kv channel

    .. Note::

       All distances must be in ANGSTROEM as this is the MDAnalysis
       default unit. All readers must return Angstroem by default.
    """
    ref_n_atoms = 1284
    ref_sum_centre_of_geometry = -98.24146
    ref_n_frames = 10

class Ref4e43(object):
    """Mixin class for a clean Protein Databank PDB file"""
    filename = PDB_full
    header = "HYDROLASE                               11-MAR-12   4E43"
    title = ["HIV PROTEASE (PR) DIMER WITH ACETATE IN EXO SITE AND PEPTIDE IN ACTIVE",
             "2 SITE"]
    compnd = ["MOL_ID: 1;",
              "2 MOLECULE: PROTEASE;",
              "3 CHAIN: A, B;",
              "4 ENGINEERED: YES;",
              "5 MUTATION: YES;",
              "6 MOL_ID: 2;",
              "7 MOLECULE: RANDOM PEPTIDE;",
              "8 CHAIN: C;",
              "9 ENGINEERED: YES;",
              "10 OTHER_DETAILS: UNKNOWN IMPURITY",
        ]
    num_remarks = 333
    # only first 5 remarks for comparison
    nmax_remarks = 5
    remarks = [
        "2",
        "2 RESOLUTION.    1.54 ANGSTROMS.",
        "3",
        "3 REFINEMENT.",
        "3   PROGRAM     : REFMAC 5.5.0110",
    ]

class TestXYZReader(TestCase, Ref2r9r):
    def setUp(self):
        self.universe = mda.Universe(XYZ_psf, XYZ)
        self.prec = 3  # 4 decimals in xyz file

    def tearDown(self):
        del self.universe

    def test_load_xyz(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms, "load Universe from PSF and XYZ")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms, "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, self.ref_n_frames, "wrong number of frames in xyz")

    def test_sum_centres_of_geometry(self):
        centreOfGeometry = 0

        for i in self.universe.trajectory:
            sel = self.universe.select_atoms("all")
            centreOfGeometry += sum(sel.center_of_geometry())

        assert_almost_equal(centreOfGeometry, self.ref_sum_centre_of_geometry, self.prec,
                            err_msg="sum of centers of geometry over the trajectory do not match")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_slice(self):
        frames = [ts.frame for ts in self.universe.trajectory[::2]]
        assert_equal(frames, np.arange(len(self.universe.trajectory))[::2])


class TestXYZReaderAsTopology(object):
    """Test that an XYZ file can act as its own topology"""
    def setUp(self):
        self.universe = mda.Universe(XYZ_mini)

    def tearDown(self):
        del self.universe

    def test_coords(self):
        ref = np.array([[0.0, 0.0, 0.0],
                        [1.0, 1.0, 1.0],
                        [2.0, 2.0, 2.0]],
                       dtype=np.float32)
        assert_array_almost_equal(self.universe.atoms.positions, ref)

    def test_rewind(self):
        self.universe.trajectory.rewind()
        assert_equal(self.universe.trajectory.ts.frame, 0, "rewinding to frame 0")

    def test_dt(self):
        assert_almost_equal(self.universe.trajectory.dt, 1.0, 4,
                            err_msg="wrong timestep dt")

class TestCompressedXYZReader(TestCase, Ref2r9r):
    def setUp(self):
        self.universe = mda.Universe(XYZ_psf, XYZ_bz2)
        self.prec = 3  # 4 decimals in xyz file

    def tearDown(self):
        del self.universe

    def test_load_xyz(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms, "load Universe from PSF and XYZ")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms, "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, self.ref_n_frames, "wrong number of frames in xyz")

    def test_sum_centres_of_geometry(self):
        centreOfGeometry = 0

        for i in self.universe.trajectory:
            sel = self.universe.select_atoms("all")
            centreOfGeometry += sum(sel.center_of_geometry())

        assert_almost_equal(centreOfGeometry, self.ref_sum_centre_of_geometry, self.prec,
                            err_msg="sum of centers of geometry over the trajectory do not match")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_slice(self):
        frames = [ts.frame for ts in self.universe.trajectory[::2]]
        assert_equal(frames, np.arange(len(self.universe.trajectory))[::2])

    def test_rewind(self):
        self.universe.trajectory.rewind()
        assert_equal(self.universe.trajectory.ts.frame, 0, "rewinding to frame 0")

    def test_next(self):
        self.universe.trajectory.rewind()
        self.universe.trajectory.next()
        assert_equal(self.universe.trajectory.ts.frame, 1, "loading frame 1")

    def test_dt(self):
        assert_almost_equal(self.universe.trajectory.dt, 1.0, 4,
                            err_msg="wrong timestep dt")


class TestXYZWriter(TestCase, Ref2r9r):
    def setUp(self):
        self.universe = mda.Universe(XYZ_psf, XYZ_bz2)
        self.prec = 3  # 4 decimals in xyz file
        ext = ".xyz"
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)
        self.Writer = MDAnalysis.coordinates.XYZ.XYZWriter

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe

    def test_write_trajectory_timestep(self):
        W = self.Writer(self.outfile)
        self._copy_traj(W)

    def test_write_trajectory_atomgroup(self):
        W = self.Writer(self.outfile)
        for ts in self.universe.trajectory:
            W.write(self.universe.atoms)
        W.close()
        self._check_copy()

    def test_ReaderWriter(self):
        t = self.universe.trajectory
        W = t.Writer(self.outfile)
        self._copy_traj(W)

    def _copy_traj(self, writer):
        for ts in self.universe.trajectory:
            writer.write_next_timestep(ts)
        writer.close()
        self._check_copy()

    def _check_copy(self):
        uw = mda.Universe(XYZ_psf, self.outfile)
        assert_equal(self.universe.trajectory.n_frames, uw.trajectory.n_frames)
        # check that the trajectories are identical for each time step
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between original and written trajectory at frame "
                                              "%d (orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))


class RefACHE(object):
    """Mixin class to provide comparison numbers.

    ACHE peptide

    # COM check in VMD::

        set p [atomselect top "not water"]
        set total {0 0 0};
        for {set i 0} {$i < 11} {incr i} {
           $p frame $i; set total [vecadd $total [measure center $p]]}

        puts [vecsum $total]
        # 472.2592159509659

    """
    ref_n_atoms = 252
    ref_proteinatoms = ref_n_atoms
    ref_sum_centre_of_geometry = 472.2592159509659  # 430.44807815551758
    ref_n_frames = 11
    ref_periodic = False


class RefCappedAla(object):
    """Mixin class to provide comparison numbers.

    Capped Ala in water

    # COM check in VMD (load trajectory as *AMBER with periodic box*!)::

        set p [atomselect top "not water"]
        set total {0 0 0};
        for {set i 0} {$i < 11} {incr i} {
           $p frame $i; set total [vecadd $total [measure center $p]]}

        puts [vecsum $total]
        # 686.276834487915

    """
    ref_n_atoms = 5071
    ref_proteinatoms = 22
    ref_sum_centre_of_geometry = 686.276834487915
    ref_n_frames = 11
    ref_periodic = True


class RefVGV(object):
    """Mixin class to provide comparison numbers.

    Computed from bala.trj::

      w = MDAnalysis.Universe(PRMncdf, TRJncdf)
      ref_n_atoms = len(w.atoms)
      ref_proteinatoms = len(w.select_atoms("protein"))
      ref_sum_centre_of_geometry = np.sum([protein.center_of_geometry() for ts in w.trajectory])
    """
    ref_n_atoms = 2661
    ref_proteinatoms = 50
    ref_sum_centre_of_geometry = 1552.9125
    ref_n_frames = 30
    ref_periodic = True


class _TRJReaderTest(TestCase):
    # use as a base class (override setUp()) and mixin a reference
    def tearDown(self):
        del self.universe

    def test_load_prm(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms, "load Universe from PRM and TRJ")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms, "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, self.ref_n_frames, "wrong number of frames in xyz")

    def test_periodic(self):
        assert_equal(self.universe.trajectory.periodic, self.ref_periodic)

    def test_amber_proteinselection(self):
        protein = self.universe.select_atoms('protein')
        assert_equal(protein.n_atoms, self.ref_proteinatoms, "error in protein selection (HIS or termini?)")

    def test_sum_centres_of_geometry(self):
        protein = self.universe.select_atoms('protein')
        total = np.sum([protein.center_of_geometry() for ts in self.universe.trajectory])
        assert_almost_equal(total, self.ref_sum_centre_of_geometry, self.prec,
                            err_msg="sum of centers of geometry over the trajectory do not match")

    def test_initial_frame_is_0(self):
        assert_equal(self.universe.trajectory.ts.frame, 0,
                     "initial frame is not 0 but {0}".format(self.universe.trajectory.ts.frame))

    def test_starts_with_first_frame(self):
        """Test that coordinate arrays are filled as soon as the trajectory has been opened."""
        assert_(np.any(self.universe.atoms.coordinates() > 0),
                "Reader does not populate coordinates() right away.")

    def test_rewind(self):
        trj = self.universe.trajectory
        trj.next()
        trj.next()  # for readers that do not support indexing
        assert_equal(trj.ts.frame, 2, "failed to forward to frame 2 (frameindex 2)")
        trj.rewind()
        assert_equal(trj.ts.frame, 0, "failed to rewind to first frame")
        assert_(np.any(self.universe.atoms.coordinates() > 0),
                "Reader does not populate coordinates() after rewinding.")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))


class TestTRJReader(_TRJReaderTest, RefACHE):
    def setUp(self):
        self.universe = mda.Universe(PRM, TRJ)
        self.prec = 3

    def test_slice_raises_TypeError(self):
        def trj_iter():
            return list(self.universe.trajectory[::2])

        assert_raises(TypeError, trj_iter)


class TestBzippedTRJReader(TestTRJReader):
    def setUp(self):
        self.universe = mda.Universe(PRM, TRJ_bz2)
        self.prec = 3


class TestBzippedTRJReaderPBC(_TRJReaderTest, RefCappedAla):
    def setUp(self):
        self.universe = mda.Universe(PRMpbc, TRJpbc_bz2)
        self.prec = 3

    def test_slice_raises_TypeError(self):
        def trj_iter():
            return list(self.universe.trajectory[::2])

        assert_raises(TypeError, trj_iter)


class TestNCDFReader(_TRJReaderTest, RefVGV):
    def setUp(self):
        self.universe = mda.Universe(PRMncdf, NCDF)
        self.prec = 3

    def test_slice_iteration(self):
        frames = [ts.frame for ts in self.universe.trajectory[4:-2:4]]
        assert_equal(frames,
                     np.arange(self.universe.trajectory.n_frames)[4:-2:4],
                     err_msg="slicing did not produce the expected frames")

    def test_metadata(self):
        data = self.universe.trajectory.trjfile
        assert_equal(data.Conventions, 'AMBER')
        assert_equal(data.ConventionVersion, '1.0')

class TestNCDFReader2(TestCase):
    """NCDF Trajectory with positions and forces.

    Contributed by Albert Solernou
    """
    def setUp(self):
        self.u = mda.Universe(PFncdf_Top, PFncdf_Trj)
        self.prec = 3

    def tearDown(self):
        self.u.trajectory.close()
        del self.u

    def test_positions_1(self):
        """Check positions on first frame"""
        self.u.trajectory[0]
        ref_1 = np.array([[ -0.11980818,  18.70524979,  11.6477766 ],
                          [ -0.44717646,  18.61727142,  12.59919548],
                          [ -0.60952115,  19.47885513,  11.22137547]], dtype=np.float32)
        assert_array_almost_equal(ref_1, self.u.atoms.positions[:3], self.prec)

    def test_positions_2(self):
        """Check positions on second frame"""
        self.u.trajectory[1]
        ref_2= np.array([[ -0.13042036,  18.6671524 ,  11.69647026],
                         [ -0.46643803,  18.60186768,  12.646698  ],
                         [ -0.46567637,  19.49173927,  11.21922874]], dtype=np.float32)
        assert_array_almost_equal(ref_2, self.u.atoms.positions[:3], self.prec)

    def test_forces_1(self):
        """Check forces on first frame"""
        self.u.trajectory[0]
        ref_1 = np.array([[ 49.23017883, -97.05565643, -86.09863281],
                          [  2.97547197,  29.84169388,  11.12069607],
                          [-15.93093777,  14.43616867,  30.25889015]], dtype=np.float32)
        assert_array_almost_equal(ref_1, self.u.atoms.forces[:3], self.prec)

    def test_forces_2(self):
        """Check forces on second frame"""
        self.u.trajectory[1]
        ref_2 = np.array([[ 116.39096832, -145.44448853, -151.3155365 ],
                          [ -18.90058327,   27.20145798,    1.95245135],
                          [ -31.08556366,   14.95863628,   41.10367966]], dtype=np.float32)
        assert_array_almost_equal(ref_2, self.u.atoms.forces[:3], self.prec)

    def test_time_1(self):
        """Check time on first frame"""
        ref = 35.02
        assert_almost_equal(ref, self.u.trajectory[0].time, self.prec)

    def test_time_2(self):
        """Check time on second frame"""
        ref = 35.04
        assert_almost_equal(ref, self.u.trajectory[1].time, self.prec)


class TestNCDFWriter(TestCase, RefVGV):
    def setUp(self):
        self.universe = mda.Universe(PRMncdf, NCDF)
        self.prec = 6
        ext = ".ncdf"
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)
        fd, self.outtop = tempfile.mkstemp(suffix=".pdb")
        os.close(fd)
        self.Writer = MDAnalysis.coordinates.TRJ.NCDFWriter

    def tearDown(self):
        for f in self.outfile, self.outtop:
            try:
                os.unlink(f)
            except OSError:
                pass
        del self.universe
        del self.Writer

    def test_write_trajectory(self):
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms, dt=t.dt)
        self._copy_traj(W)

    def test_OtherWriter(self):
        t = self.universe.trajectory
        W = t.OtherWriter(self.outfile)
        self._copy_traj(W)

    def _copy_traj(self, writer):
        for ts in self.universe.trajectory:
            writer.write_next_timestep(ts)
        writer.close()

        uw = mda.Universe(PRMncdf, self.outfile)

        # check that the trajectories are identical for each time step
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between original and written trajectory at frame "
                                              "%d (orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))
            # not a good test because in the example trajectory all times are 0
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the same.".format(orig_ts.frame))
            assert_array_almost_equal(written_ts.dimensions, orig_ts.dimensions, self.prec,
                                      err_msg="unitcells are not identical")

    @attr('slow')
    def test_TRR2NCDF(self):
        trr = MDAnalysis.Universe(GRO, TRR)
        W = self.Writer(self.outfile, trr.trajectory.n_atoms, velocities=True)
        for ts in trr.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = MDAnalysis.Universe(GRO, self.outfile)

        for orig_ts, written_ts in itertools.izip(trr.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between original and written trajectory at frame "
                                              "%d (orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))
            assert_array_almost_equal(written_ts._velocities, orig_ts._velocities, self.prec,
                                      err_msg="velocity mismatch between original and written trajectory at frame %d "
                                              "(orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the same.".format(orig_ts.frame))
            assert_array_almost_equal(written_ts.dimensions, orig_ts.dimensions, self.prec,
                                      err_msg="unitcells are not identical")
        del trr

    @attr('issue')
    def test_write_AtomGroup(self):
        """test to write NCDF from AtomGroup (Issue 116)"""
        p = self.universe.select_atoms("not resname WAT")
        p.write(self.outtop)
        W = self.Writer(self.outfile, n_atoms=p.n_atoms)
        for ts in self.universe.trajectory:
            W.write(p)
        W.close()

        uw = MDAnalysis.Universe(self.outtop, self.outfile)
        pw = uw.atoms

        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(p.positions, pw.positions, self.prec,
                                      err_msg="coordinate mismatch between original and written trajectory at frame "
                                              "%d (orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the same.".format(orig_ts.frame))
            assert_array_almost_equal(written_ts.dimensions, orig_ts.dimensions, self.prec,
                                      err_msg="unitcells are not identical")


class TestINPCRDReader(TestCase):
    """Test reading Amber restart coordinate files"""
    def _check_ts(self, ts):
        # Check a ts has the right values in
        ref_pos = np.array([[6.6528795, 6.6711416, -8.5963255],
                        [7.3133773, 5.8359736, -8.8294175],
                        [8.3254058, 6.2227613, -8.7098593],
                        [7.0833200, 5.5038197, -9.8417650],
                        [7.1129439, 4.6170351, -7.9729560]])
        for ref, val in itertools.izip(ref_pos, ts._pos):
            assert_allclose(ref, val)

    def test_reader(self):
        from MDAnalysis.coordinates.INPCRD import INPReader

        r = INPReader(INPCRD)

        assert_equal(r.n_atoms, 5)
        self._check_ts(r.ts)

    def test_universe_inpcrd(self):
        u = MDAnalysis.Universe(XYZ_five, INPCRD)

        self._check_ts(u.trajectory.ts)

    def test_universe_restrt(self):
        u = MDAnalysis.Universe(XYZ_five, INPCRD, format='RESTRT')
        self._check_ts(u.trajectory.ts)


class TestNCDFWriterVelsForces(TestCase):
    """Test writing NCDF trajectories with a mixture of options"""
    def setUp(self):
        fd, self.outfile = tempfile.mkstemp(suffix='.ncdf')
        os.close(fd)
        self.prec = 3
        self.top = XYZ_mini
        self.n_atoms = 3

        self.ts1 = MDAnalysis.coordinates.TRJ.Timestep(self.n_atoms, velocities=True, forces=True)
        self.ts1._pos[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms, 3)
        self.ts1._velocities[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms, 3) + 100
        self.ts1._forces[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms, 3) + 200

        self.ts2 = MDAnalysis.coordinates.TRJ.Timestep(self.n_atoms, velocities=True, forces=True)
        self.ts2._pos[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms, 3) + 300
        self.ts2._velocities[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms, 3) + 400
        self.ts2._forces[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms, 3) + 500

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass

        del self.n_atoms
        del self.ts1
        del self.ts2

    def _write_ts(self, pos, vel, force):
        """Write the two reference timesteps, then open them up and check values

        pos vel and force are bools which define whether these properties should
        be in TS
        """
        with MDAnalysis.Writer(self.outfile, n_atoms=self.n_atoms,
                               velocities=vel, forces=force) as w:
            w.write(self.ts1)
            w.write(self.ts2)

        u = MDAnalysis.Universe(self.top, self.outfile)
        for ts, ref_ts in itertools.izip(u.trajectory, [self.ts1, self.ts2]):
            if pos:
                assert_almost_equal(ts._pos, ref_ts._pos, self.prec)
            else:
                assert_raises(NoDataError, getattr, ts, 'positions')
            if vel:
                assert_almost_equal(ts._velocities, ref_ts._velocities, self.prec)
            else:
                assert_raises(NoDataError, getattr, ts, 'velocities')
            if force:
                assert_almost_equal(ts._forces, ref_ts._forces, self.prec)
            else:
                assert_raises(NoDataError, getattr, ts, 'forces')

        u.trajectory.close()

    def test_pos(self):
        self._write_ts(True, False, False)

    def test_pos_vel(self):
        self._write_ts(True, True, False)

    def test_pos_force(self):
        self._write_ts(True, False, True)

    def test_pos_vel_force(self):
        self._write_ts(True, True, True)


class _SingleFrameReader(TestCase, RefAdKSmall):
    # see TestPDBReader how to set up!

    def tearDown(self):
        del self.universe

    def test_flag_permissive_pdb_reader(self):
        """test_flag_permissive_pdb_reader: permissive_pdb_reader==True enables primitive PDB reader"""
        assert_equal(mda.core.flags['permissive_pdb_reader'], True,
                     "'permissive_pdb_reader' flag should be True as MDAnalysis default")

    def test_load_file(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms, "load Universe from file %s" % U.trajectory.filename)
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms, "wrong number of atoms")

    def test_numres(self):
        assert_equal(self.universe.atoms.n_residues, 214, "wrong number of residues")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1, "wrong number of frames in pdb")

    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0, "wrong time of the frame")

    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0,
                     "wrong frame number (0-based, should be 0 for single frame readers)")

    def test_frame_index_0(self):
        self.universe.trajectory[0]
        assert_equal(self.universe.trajectory.ts.frame, 0, "frame number for frame index 0 should be 0")

    def test_frame_index_1_raises_IndexError(self):
        def go_to_2(traj=self.universe.trajectory):
            traj[1]

        assert_raises(IndexError, go_to_2)

    def test_dt(self):
        """testing that accessing universe.trajectory.dt gives 1.0 (the default)"""
        assert_equal(1.0, self.universe.trajectory.dt)

    def test_coordinates(self):
        A10CA = self.universe.atoms.CA[10]
        # restrict accuracy to maximum in PDB files (3 decimals)
        assert_almost_equal(A10CA.pos, self.ref_coordinates['A10CA'], 3,
                            err_msg="wrong coordinates for A10:CA")

    def test_distances(self):
        NTERM = self.universe.atoms.N[0]
        CTERM = self.universe.atoms.C[-1]
        d = atom_distance(NTERM, CTERM)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec,
                            err_msg="distance between M1:N and G214:C")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_last_slice(self):
        trj_iter = self.universe.trajectory[-1:]  # should be same as above: only 1 frame!
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))


class TestPDBReader(_SingleFrameReader):
    def setUp(self):
        ##mda.core.flags['permissive_pdb_reader'] = False # enable Bio.PDB reader!!
        # use permissive=False instead of changing the global flag as this
        # can lead to race conditions when testing in parallel
        self.universe = mda.Universe(RefAdKSmall.filename, permissive=False)
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

    def test_uses_Biopython(self):
        from MDAnalysis.coordinates.PDB import PDBReader

        assert_(isinstance(self.universe.trajectory, PDBReader), "failed to choose Biopython PDBReader")

    @knownfailure("Biopython PDB reader does not parse CRYST1", AssertionError)
    def test_dimensions(self):
        assert_almost_equal(self.universe.trajectory.ts.dimensions, RefAdKSmall.ref_unitcell,
                            self.prec,
                            "Biopython reader failed to get unitcell dimensions from CRYST1")

class _PDBMetadata(TestCase, Ref4e43):
    permissive = True

    def setUp(self):
        self.universe = mda.Universe(self.filename, permissive=self.permissive)

    def tearDown(self):
        del self.universe

    def test_HEADER(self):
        assert_equal(self.universe.trajectory.header, self.header,
                     err_msg="HEADER record not correctly parsed")

    def test_TITLE(self):
        try:
            title = self.universe.trajectory.title
        except AttributeError:
            raise AssertionError("Reader does not have a 'title' attribute.")
        assert_equal(len(title), len(self.title),
                     err_msg="TITLE does not contain same number of lines")
        for lineno, (parsed, reference) in enumerate(zip(title, self.title), start=1):
            assert_equal(parsed, reference,
                         err_msg="TITLE line {0} do not match".format(lineno))

    def test_COMPND(self):
        try:
            compound = self.universe.trajectory.compound
        except AttributeError:
            raise AssertionError("Reader does not have a 'compound' attribute.")
        assert_equal(len(compound), len(self.compnd),
                     err_msg="COMPND does not contain same number of lines")
        for lineno, (parsed, reference) in enumerate(zip(compound, self.compnd), start=1):
            assert_equal(parsed, reference,
                         err_msg="COMPND line {0} do not match".format(lineno))

    def test_REMARK(self):
        try:
            remarks = self.universe.trajectory.remarks
        except AttributeError:
            raise AssertionError("Reader does not have a 'remarks' attribute.")
        assert_equal(len(remarks), self.num_remarks,
                     err_msg="REMARK does not contain same number of lines")
        # only look at the first 5 entries
        for lineno, (parsed, reference) in enumerate(
                zip(remarks[:self.nmax_remarks], self.remarks[:self.nmax_remarks]), start=1):
            assert_equal(parsed, reference,
                         err_msg="REMARK line {0} do not match".format(lineno))

class TestPrimitivePDBReader_Metadata(_PDBMetadata):
    permissive = True

### Does not implement Reader.remarks, Reader.header, Reader.title, Reader.compounds
### because the PDB header data in trajectory.metadata are already parsed; should perhaps
### update the PrimitivePDBReader to do the same. [orbeckst]
#class TestPDBReader_Metadata(_PDBMetadata):
#    permissive = False

class TestPSF_CRDReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(PSF, CRD)
        self.prec = 5  # precision in CRD (at least we are writing %9.5f)


class TestPSF_PDBReader(TestPDBReader):
    def setUp(self):
        # mda.core.flags['permissive_pdb_reader'] = False
        self.universe = mda.Universe(PSF, PDB_small, permissive=False)
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

    def test_uses_Biopython(self):
        from MDAnalysis.coordinates.PDB import PDBReader

        assert_(isinstance(self.universe.trajectory, PDBReader), "failed to choose Biopython PDBReader")


class TestPrimitivePDBReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(PDB_small, permissive=True)
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

    def test_missing_natoms(self):
        from MDAnalysis.coordinates.PDB import PrimitivePDBReader

        assert_raises(ValueError, PrimitivePDBReader, 'something.pdb')

    def test_wrong_natoms(self):
        from MDAnalysis.coordinates.PDB import PrimitivePDBReader

        assert_raises(ValueError, PrimitivePDBReader, PDB_small, n_atoms=4000)


class TestExtendedPDBReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(PDB_small, topology_format="XPDB", format="XPDB")
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

    def test_long_resSeq(self):
        #it checks that it can read a 5-digit resid
        self.universe = mda.Universe(XPDB_small, topology_format="XPDB")
        u = self.universe.select_atoms('resid 1 or resid 10 or resid 100 or resid 1000 or resid 10000')
        assert_equal(u[4].resid, 10000, "can't read a five digit resid")


class TestPSF_PrimitivePDBReader(TestPrimitivePDBReader):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small, permissive=True)
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

    def test_dimensions(self):
        assert_almost_equal(self.universe.trajectory.ts.dimensions, RefAdKSmall.ref_unitcell,
                            self.prec,
                            "Primitive PDB reader failed to get unitcell dimensions from CRYST1")


class TestPrimitivePDBWriter(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small, permissive=True)
        self.universe2 = mda.Universe(PSF, DCD, permissive=True)
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        ext = ".pdb"
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe, self.universe2

    def test_writer(self):
        "Test writing from a single frame PDB file to a PDB file."""
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(PSF, self.outfile, permissive=True)
        assert_almost_equal(u.atoms.coordinates(), self.universe.atoms.coordinates(), self.prec,
                            err_msg="Writing PDB file with PrimitivePDBWriter does not reproduce original coordinates")

    @attr('issue')
    def test_write_single_frame_Writer(self):
        """Test writing a single frame from a DCD trajectory to a PDB using MDAnalysis.Writer (Issue 105)"""
        u = self.universe2
        W = mda.Writer(self.outfile)
        u.trajectory[50]
        W.write(u.select_atoms('all'))
        W.close()
        u2 = mda.Universe(self.outfile)
        assert_equal(u2.trajectory.n_frames, 1, err_msg="The number of frames should be 1.")

    @attr('issue')
    def test_write_single_frame_AtomGroup(self):
        """Test writing a single frame from a DCD trajectory to a PDB using AtomGroup.write() (Issue 105)"""
        u = self.universe2
        u.trajectory[50]
        u.atoms.write(self.outfile)
        u2 = MDAnalysis.Universe(PSF, self.outfile)
        assert_equal(u2.trajectory.n_frames, 1, err_msg="Output PDB should only contain a single frame")
        assert_almost_equal(u2.atoms.coordinates(), u.atoms.coordinates(), self.prec,
                            err_msg="Written coordinates do not agree with original coordinates from frame %d" %
                                    u.trajectory.frame)

    @attr('issue')
    def test_check_coordinate_limits_min(self):
        """Test that illegal PDB coordinates (x <= -999.9995 A) are caught with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up parallel tests
        u = mda.Universe(PSF, PDB_small, permissive=True)
        u.atoms[2000].pos[1] = -999.9995
        assert_raises(ValueError, u.atoms.write, self.outfile)
        del u

    @attr('issue')
    def test_check_coordinate_limits_max(self):
        """Test that illegal PDB coordinates (x > 9999.9995 A) are caught with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up parallel tests
        u = mda.Universe(PSF, PDB_small, permissive=True)
        u.atoms[1000].pos[1] = 9999.9996  # OB: 9999.99951 is not caught by '<=' ?!?
        assert_raises(ValueError, u.atoms.write, self.outfile)
        del u


class TestGMSReader(TestCase):
    ''' Test cases for GAMESS output log-files '''

    def setUp(self):
        # optimize no-symmetry
        self.u_aso = mda.Universe(GMS_ASYMOPT, GMS_ASYMOPT, format='GMS',
                topology_format='GMS')
        self.u_so =  mda.Universe(GMS_SYMOPT,  GMS_SYMOPT)
        self.u_ass = mda.Universe(GMS_ASYMSURF, GMS_ASYMSURF)

    def test_n_frames(self):
        desired = [21,8,10]
        assert_equal(self.u_aso.trajectory.n_frames, desired[0],
                err_msg="Wrong number of frames read from GAMESS C1 optimization")
        assert_equal(self.u_so.trajectory.n_frames, desired[1],
                err_msg="Wrong number of frames read from GAMESS D4H optimization")
        assert_equal(self.u_ass.trajectory.n_frames, desired[2],
                err_msg="Wrong number of frames read from GAMESS C1 surface")

    def test_step5distances_asymopt(self):
        '''TestGMSReader: C1 optimization:
            distance between 1st and 4th atoms changes after 5 steps '''
        desired = -0.0484664
        assert_almost_equal(self.__calcFD(self.u_aso), desired, decimal=5,
                err_msg="Wrong 1-4 atom distance change after 5 steps for GAMESS C1 optimization")

    def test_step5distances_symopt(self):
        '''TestGMSReader: Symmetry-input optimization:
            distance between 1st and 4th atoms changes after 5 steps '''
        desired = 0.227637
        assert_almost_equal(self.__calcFD(self.u_so), desired, decimal=5,
                err_msg="Wrong 1-4 atom distance change after 5 steps for GAMESS D4H optimization")

    def test_step5distances_asymsurf(self):
        '''TestGMSReader: Symmetry-input potential-energy surface:
            distance between 1st and 4th atoms changes after 5 steps '''
        desired = -0.499996
        assert_almost_equal(self.__calcFD(self.u_ass), desired, decimal=5,
                err_msg="Wrong 1-4 atom distance change after 5 steps for GAMESS C1 surface")

    def __calcFD(self, u):
        u.trajectory.rewind()
        pp = (u.trajectory.ts._pos[0] - u.trajectory.ts._pos[3])
        z1 = np.sqrt(sum(pp**2))
        for i in range(5):
            u.trajectory.next()
        pp = (u.trajectory.ts._pos[0] - u.trajectory.ts._pos[3])
        z2 = np.sqrt(sum(pp**2))
        return z1-z2

    def test_rewind(self):
        self.u_aso.trajectory.rewind()
        assert_equal(self.u_aso.trajectory.ts.frame, 0, "rewinding to frame 0")

    def test_next(self):
        self.u_aso.trajectory.rewind()
        self.u_aso.trajectory.next()
        assert_equal(self.u_aso.trajectory.ts.frame, 1, "loading frame 1")

    def test_dt(self):
        assert_almost_equal(self.u_aso.trajectory.dt, 1.0, 4,
                            err_msg="wrong timestep dt")

    def tearDown(self):
        del self.u_aso
        del self.u_so
        del self.u_ass


class TestMultiPDBReader(TestCase):
    def setUp(self):
        self.multiverse = mda.Universe(PDB_multiframe, permissive=True, guess_bonds=True)
        self.multiverse.build_topology()
        self.conect = mda.Universe(CONECT, guess_bonds=True)
        self.conect.build_topology()

    def tearDown(self):
        del self.multiverse
        del self.conect

    @attr('slow')
    def test_n_frames(self):
        assert_equal(self.multiverse.trajectory.n_frames, 24,
                     "Wrong number of frames read from PDB muliple model file")

    @attr('slow')
    def test_n_atoms_frame(self):
        u = self.multiverse
        desired = 392
        for frame in u.trajectory:
            assert_equal(len(u.atoms), desired,
                         err_msg=("The number of atoms in the Universe (%d) does not"
                                  " match the number of atoms in the test case (%d) at frame %d"
                                  % (len(u.atoms), desired, u.trajectory.frame)))

    @attr('slow')
    def test_rewind(self):
        u = self.multiverse
        u.trajectory[11]
        assert_equal(u.trajectory.ts.frame, 11, "Failed to forward to 11th frame (frame index 11)")
        u.trajectory.rewind()
        assert_equal(u.trajectory.ts.frame, 0, "Failed to rewind to 0th frame (frame index 0)")

    @attr('slow')
    def test_iteration(self):
        u = self.multiverse
        frames = []
        for frame in u.trajectory:
            pass
        # should rewind after previous test
        # problem was: the iterator is NoneType and next() cannot be called
        for ts in u.trajectory:
            frames.append(ts)
        assert_equal(len(frames), u.trajectory.n_frames,
                     "iterated number of frames %d is not the expected number %d; "
                     "trajectory iterator fails to rewind" % (len(frames), u.trajectory.n_frames))

    @attr('slow')
    def test_slice_iteration(self):
        u = self.multiverse
        frames = []
        for ts in u.trajectory[4:-2:4]:
            frames.append(ts.frame)
        assert_equal(np.array(frames),
                     np.arange(u.trajectory.n_frames)[4:-2:4],
                     err_msg="slicing did not produce the expected frames")

    @attr('slow')
    def test_conect_bonds_conect(self):
        conect = self.conect
        assert_equal(len(conect.atoms), 1890)
        assert_equal(len(conect.bonds), 1922)

        try:
            fd, outfile = tempfile.mkstemp(suffix=".pdb")
            os.close(fd)
            self.conect.atoms.write(outfile, bonds="conect")
            u1 = mda.Universe(outfile, guess_bonds=True)
        finally:
            os.unlink(outfile)
        assert_equal(len(u1.atoms), 1890)
        assert_equal(len(u1.bonds), 1922)

    @attr('slow')
    def test_conect_bonds_all(self):
        conect = self.conect
        assert_equal(len(conect.atoms), 1890)
        assert_equal(len(conect.bonds), 1922)

        try:
            fd, outfile = tempfile.mkstemp(suffix=".pdb")
            os.close(fd)
            self.conect.atoms.write(outfile, bonds="all")
            u2 = mda.Universe(outfile, guess_bonds=True)
        finally:
            os.unlink(outfile)
        assert_equal(len(u2.atoms), 1890)
        assert_equal(len([b for b in u2.bonds if not b.is_guessed]), 1922)

        #assert_equal(len([b for b in conect.bonds if not b.is_guessed]), 1922)

    @attr('slow')
    def test_numconnections(self):
        u = self.multiverse

        # the bond list is sorted - so swaps in input pdb sequence should not be a problem
        desired = [
            [48, 365],
            [99, 166],
            [166, 99],
            [249, 387],
            [313, 331],
            [331, 313, 332, 340],
            [332, 331, 333, 338, 341],
            [333, 332, 334, 342, 343],
            [334, 333, 335, 344, 345],
            [335, 334, 336, 337],
            [336, 335],
            [337, 335, 346, 347, 348],
            [338, 332, 339, 349],
            [339, 338],
            [340, 331],
            [341, 332],
            [342, 333],
            [343, 333],
            [344, 334],
            [345, 334],
            [346, 337],
            [347, 337],
            [348, 337],
            [349, 338],
            [365, 48],
            [387, 249]]

        def helper(atoms, bonds):
            """
            Convert a bunch of atoms and bonds into a list of CONECT records
            """
            con = {}

            for bond in bonds:
                a1, a2 = bond[0].index, bond[1].index
                if not a1 in con:
                    con[a1] = []
                if not a2 in con:
                    con[a2] = []
                con[a2].append(a1)
                con[a1].append(a2)

            #print con
            atoms = sorted([a.index for a in atoms])

            conect = [([a, ] + sorted(con[a])) for a in atoms if a in con]
            conect = [[a + 1 for a in c] for c in conect]

            return conect

        conect = helper(self.multiverse.atoms, [b for b in u.bonds if not b.is_guessed])
        #for r in conect:
            #print r
        assert_equal(conect, desired,
                     err_msg="The bond list does not match the test reference; len(actual) is %d, len(desired) is %d "
                             "" % (
                     len(u._topology['bonds']), len(desired)))


class TestMultiPDBWriter(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small, permissive=True)
        self.multiverse = mda.Universe(PDB_multiframe, permissive=True)
        self.universe2 = mda.Universe(PSF, DCD, permissive=True)
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        ext = ".pdb"
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)
        fd, self.outfile2 = tempfile.mkstemp(suffix=ext)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        try:
            os.unlink(self.outfile2)
        except OSError:
            pass
        del self.universe, self.multiverse, self.universe2

    @attr('slow')
    def test_write_atomselection(self):
        """Test if multiframe writer can write selected frames for an atomselection."""
        u = self.multiverse
        group = u.select_atoms('name CA', 'name C')
        desired_group = 56
        desired_frames = 6
        pdb = MDAnalysis.Writer(self.outfile, multiframe=True, start=12, step=2)
        for ts in u.trajectory[-6:]:
            pdb.write(group)
        pdb.close()
        u2 = mda.Universe(self.outfile)
        assert_equal(len(u2.atoms), desired_group,
                     err_msg="MultiPDBWriter trajectory written for an AtomGroup contains %d atoms, it should contain "
                             "%d" % (
                     len(u2.atoms), desired_group))

        assert_equal(len(u2.trajectory), desired_frames,
                     err_msg="MultiPDBWriter trajectory written for an AtomGroup contains %d frames, it should have "
                             "%d" % (
                     len(u.trajectory), desired_frames))

    @attr('slow')
    def test_write_all_timesteps(self):
        """
        Test write_all_timesteps() of the  multiframe writer (selected frames for an atomselection)
        """
        u = self.multiverse
        group = u.select_atoms('name CA', 'name C')
        desired_group = 56
        desired_frames = 6

        pdb = MDAnalysis.Writer(self.outfile, multiframe=True, start=12, step=2)
        pdb.write_all_timesteps(group)
        u2 = mda.Universe(self.outfile)
        assert_equal(len(u2.atoms), desired_group,
                     err_msg="MultiPDBWriter trajectory written for an AtomGroup contains %d atoms, it should contain "
                             "%d" % (
                     len(u2.atoms), desired_group))

        assert_equal(len(u2.trajectory), desired_frames,
                     err_msg="MultiPDBWriter trajectory written for an AtomGroup contains %d frames, it should have "
                             "%d" % (
                     len(u.trajectory), desired_frames))

    @attr('slow')
    def test_write_atoms(self):
        u = self.universe2
        W = mda.Writer(self.outfile, multiframe=True)
        # 2 frames expceted
        for ts in u.trajectory[-2:]:
            W.write(u.atoms)
        W.close()
        u0 = mda.Universe(self.outfile)
        assert_equal(u0.trajectory.n_frames, 2, err_msg="The number of frames should be 3.")


class TestPQRReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(PQR)
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

    def test_total_charge(self):
        assert_almost_equal(self.universe.atoms.total_charge(), self.ref_charmm_totalcharge, 3,
                            "Total charge (in CHARMM) does not match expected value.")

    def test_hydrogenCharges(self):
        assert_almost_equal(self.universe.atoms.H.charges, self.ref_charmm_Hcharges, 3,
                            "Charges for H atoms do not match.")

    # Note that the whole system gets the sysID 'SYSTEM' for the PQR file (when read with
    # a PSF it is 's4AKE')
    def test_ArgCACharges(self):
        assert_almost_equal(self.universe.SYSTEM.ARG.CA.charges, self.ref_charmm_ArgCAcharges, 3,
                            "Charges for CA atoms in Arg residues do not match.")

    def test_ProNCharges(self):
        assert_almost_equal(self.universe.SYSTEM.PRO.N.charges, self.ref_charmm_ProNcharges, 3,
                            "Charges for N atoms in Pro residues do not match.")


class TestPQRWriter(TestCase, RefAdKSmall):
    def setUp(self):
        self.universe = mda.Universe(PQR)
        self.prec = 3
        ext = ".pqr"
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe

    def test_writer_noChainID(self):
        assert_equal(self.universe.segments.segids[0], 'SYSTEM')  # sanity check
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_equal(u.segments.segids[0], 'SYSTEM')
        assert_almost_equal(u.atoms.coordinates(), self.universe.atoms.coordinates(), self.prec,
                            err_msg="Writing PQR file with PQRWriter does not reproduce original coordinates")
        assert_almost_equal(u.atoms.charges, self.universe.atoms.charges, self.prec,
                            err_msg="Writing PQR file with PQRWriter does not reproduce original charges")
        assert_almost_equal(u.atoms.radii, self.universe.atoms.radii, self.prec,
                            err_msg="Writing PQR file with PQRWriter does not reproduce original radii")

    def test_write_withChainID(self):
        self.universe.atoms.set_segids('A')
        assert_equal(self.universe.segments.segids[0], 'A')  # sanity check
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_equal(u.segments.segids[0], 'A')
        assert_almost_equal(u.atoms.coordinates(), self.universe.atoms.coordinates(), self.prec,
                            err_msg="Writing PQR file with PQRWriter does not reproduce original coordinates")
        assert_almost_equal(u.atoms.charges, self.universe.atoms.charges, self.prec,
                            err_msg="Writing PQR file with PQRWriter does not reproduce original charges")
        assert_almost_equal(u.atoms.radii, self.universe.atoms.radii, self.prec,
                            err_msg="Writing PQR file with PQRWriter does not reproduce original radii")

    def test_timestep_not_modified_by_writer(self):
        ts = self.universe.trajectory.ts
        x = ts._pos.copy()
        self.universe.atoms.write(self.outfile)
        assert_equal(ts._pos, x, err_msg="Positions in Timestep were modified by writer.")

    def test_total_charge(self):
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_almost_equal(u.atoms.total_charge(), self.ref_charmm_totalcharge, 3,
                            "Total charge (in CHARMM) does not match expected value.")


class TestGROReader(TestCase, RefAdK):
    def setUp(self):
        self.universe = mda.Universe(GRO)
        self.ts = self.universe.trajectory.ts
        self.prec = 2  # lower prec in gro!! (3 decimals nm -> 2 decimals in Angstroem)

    def tearDown(self):
        del self.universe
        del self.ts

    def test_flag_convert_lengths(self):
        assert_equal(mda.core.flags['convert_lengths'], True,
                     "MDAnalysis.core.flags['convert_lengths'] should be True by default")

    def test_load_gro(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms, "load Universe from small GRO")
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms, "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1, "wrong number of frames")

    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0, "wrong time of the frame")

    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0, "wrong frame number (should be 0 for 0-based ts.frame)")

    def test_dt(self):
        """testing that accessing universe.trajectory.dt gives default of 1.0"""
        assert_equal(self.universe.trajectory.dt, 1.0)

    def test_coordinates(self):
        A10CA = self.universe.SYSTEM.CA[10]
        assert_almost_equal(A10CA.pos, self.ref_coordinates['A10CA'], self.prec,
                            err_msg="wrong coordinates for A10:CA")

    def test_distances(self):
        # NOTE that the prec is only 1 decimal: subtracting two low precision coordinates
        #      low prec: 9.3455122920041109; high prec (from pdb): 9.3513174
        NTERM = self.universe.SYSTEM.N[0]
        CTERM = self.universe.SYSTEM.C[-1]
        d = atom_distance(NTERM, CTERM)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec - 1,  # note low prec!!
                            err_msg="distance between M1:N and G214:C")

    def test_selection(self):
        na = self.universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size, "Atom selection of last atoms in file")

    def test_unitcell(self):
        assert_array_almost_equal(self.ts.dimensions, self.ref_unitcell, self.prec,
                                  err_msg="unit cell dimensions (rhombic dodecahedron)")

    def test_volume(self):
        # test_volume: reduce precision for Gromacs comparison to 0 decimals (A**3 <--> nm**3!)
        assert_almost_equal(self.ts.volume, self.ref_volume, 0,
                            err_msg="wrong volume for unitcell (rhombic dodecahedron)")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))


class TestDMSReader(TestCase):
    def setUp(self):
        self.universe = mda.Universe(DMS)
        self.ts = self.universe.trajectory.ts

    def tearDown(self):
        del self.universe
        del self.ts

    def test_global_cell(self):
        assert_equal(self.ts.dimensions, [0., 0., 0., 0., 0., 0.])

    def test_velocities(self):
        assert_equal(hasattr(self.ts, "_velocities"), False)

    def test_number_of_coords(self):
        # Desired value taken from VMD
        #      Info)    Atoms: 3341
        assert_equal(len(self.universe.atoms), 3341)

    def test_coords_atom_0(self):
        # Desired coordinates taken directly from the SQLite file. Check unit conversion
        coords_0 = np.array([-11.0530004501343, 26.6800003051758, 12.7419996261597, ], dtype=np.float32)
        assert_array_equal(self.universe.atoms[0].pos, coords_0)

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1, "wrong number of frames in pdb")

    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0, "wrong time of the frame")

    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0,
                     "wrong frame number (0-based, should be 0 for single frame readers)")

    def test_frame_index_0(self):
        self.universe.trajectory[0]
        assert_equal(self.universe.trajectory.ts.frame, 0, "frame number for frame index 0 should be 0")

    def test_frame_index_1_raises_IndexError(self):
        def go_to_2(traj=self.universe.trajectory):
            traj[1]

        assert_raises(IndexError, go_to_2)

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
        assert_almost_equal(A10CA.pos, RefAdK.ref_coordinates['A10CA'] / 10.0,  # coordinates in nm
                            self.prec,
                            err_msg="wrong native coordinates (in nm) for A10:CA")

    def test_distances(self):
        # 3 decimals on nm in gro but we compare to the distance
        # computed from the pdb file, so the effective precision is 2 again.
        # (Otherwise the distance test fails:
        #  Arrays are not almost equal distance between M1:N and G214:C
        #    ACTUAL: 0.93455122920041123
        #    DESIRED: 0.93513173999999988
        NTERM = self.universe.SYSTEM.N[0]
        CTERM = self.universe.SYSTEM.C[-1]
        d = atom_distance(NTERM, CTERM)
        assert_almost_equal(d, RefAdK.ref_distances['endtoend'] / 10.0,  # coordinates in nm
                            self.prec - 1,
                            err_msg="distance between M1:N and G214:C")

    def test_unitcell(self):
        # lengths in A : convert to nm
        assert_array_almost_equal(self.ts.dimensions[:3], self.ref_unitcell[:3] / 10.0, self.prec,
                                  err_msg="unit cell A,B,C (rhombic dodecahedron)")
        # angles should not have changed
        assert_array_almost_equal(self.ts.dimensions[3:], self.ref_unitcell[3:], self.prec,
                                  err_msg="unit cell alpha,beta,gamma (rhombic dodecahedron)")

    def test_volume(self):
        # ref lengths in A (which was originally converted from nm)
        assert_almost_equal(self.ts.volume, self.ref_volume / 1000., 3,
                            err_msg="wrong volume for unitcell (rhombic dodecahedron)")


class TestGROWriter(TestCase):
    def setUp(self):
        self.universe = mda.Universe(GRO)
        self.prec = 2  # 3 decimals in file in nm but MDAnalysis is in A
        ext = ".gro"
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)
        fd, self.outfile2 = tempfile.mkstemp(suffix=ext)
        os.close(fd)

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

    @dec.slow
    def test_writer(self):
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_almost_equal(u.atoms.coordinates(), self.universe.atoms.coordinates(), self.prec,
                            err_msg="Writing GRO file with GROWriter does not reproduce original coordinates")

    @dec.slow
    def test_timestep_not_modified_by_writer(self):
        ts = self.universe.trajectory.ts
        x = ts._pos.copy()
        self.universe.atoms.write(self.outfile)
        assert_equal(ts._pos, x, err_msg="Positions in Timestep were modified by writer.")

    @dec.slow
    @attr('issue')
    def test_check_coordinate_limits_min(self):
        """Test that illegal GRO coordinates (x <= -999.9995 nm) are caught with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up parallel tests
        u = mda.Universe(GRO)
        u.atoms[2000].pos[1] = -999.9995 * 10  # nm -> A
        assert_raises(ValueError, u.atoms.write, self.outfile2)
        del u

    @dec.slow
    @attr('issue')
    def test_check_coordinate_limits_max(self):
        """Test that illegal GRO coordinates (x > 9999.9995 nm) are caught with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up parallel tests
        u = mda.Universe(GRO)
        u.atoms[1000].pos[1] = 9999.9999 * 10  # nm -> A  ; [ob] 9999.9996 not caught
        assert_raises(ValueError, u.atoms.write, self.outfile2)
        del u

    @dec.slow
    def test_check_coordinate_limits_max_noconversion(self):
        """Test that illegal GRO coordinates (x > 9999.9995 nm) also raises exception for convert_units=False"""
        # modify coordinates so we need our own copy or we could mess up parallel tests
        u = mda.Universe(GRO, convert_units=False)
        u.atoms[1000].pos[1] = 9999.9999
        assert_raises(ValueError, u.atoms.write, self.outfile2, convert_units=False)
        del u


class TestPDBReaderBig(TestCase, RefAdK):
    def setUp(self):
        self.universe = mda.Universe(PDB)
        self.prec = 6

    def tearDown(self):
        del self.universe

    @dec.slow
    def test_load_pdb(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms, "load Universe from big PDB")
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    @dec.slow
    def test_selection(self):
        na = self.universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size, "Atom selection of last atoms in file")

    @dec.slow
    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms, "wrong number of atoms")

    @dec.slow
    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1, "wrong number of frames")

    @dec.slow
    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0, "wrong time of the frame")

    @dec.slow
    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0, "wrong frame number")

    @dec.slow
    def test_dt(self):
        """testing that accessing universe.trajectory.dt returns the default of 1.0 ps"""
        assert_equal(self.universe.trajectory.dt, 1.0)

    @dec.slow
    def test_coordinates(self):
        A10CA = self.universe.SYSTEM.CA[10]
        assert_almost_equal(A10CA.pos, self.ref_coordinates['A10CA'], self.prec,
                            err_msg="wrong coordinates for A10:CA")

    @dec.slow
    def test_distances(self):
        NTERM = self.universe.SYSTEM.N[0]
        CTERM = self.universe.SYSTEM.C[-1]
        d = atom_distance(NTERM, CTERM)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec,
                            err_msg="wrong distance between M1:N and G214:C")

    @dec.slow
    def test_selection(self):
        na = self.universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size, "Atom selection of last atoms in file")

    @dec.slow
    @attr('issue')
    def test_unitcell(self):
        assert_array_almost_equal(self.universe.coord.dimensions, self.ref_unitcell, self.prec,
                                  err_msg="unit cell dimensions (rhombic dodecahedron), issue 60")

    @dec.slow
    def test_volume(self):
        assert_almost_equal(self.universe.coord.volume, self.ref_volume, 0,
                            err_msg="wrong volume for unitcell (rhombic dodecahedron)")


@attr('issue')
def TestDCD_Issue32():
    """Test for Issue 32: 0-size dcds lead to a segfault: now caught with IOError"""
    assert_raises(IOError, mda.Universe, PSF, DCD_empty)


class _TestDCD(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.dcd = self.universe.trajectory
        self.ts = self.universe.coord

    def tearDown(self):
        del self.universe
        del self.dcd
        del self.ts


class TestDCDReaderClass(TestCase):
    def test_with_statement(self):
        from MDAnalysis.coordinates.DCD import DCDReader

        try:
            with DCDReader(DCD) as trj:
                N = trj.n_frames
                frames = [ts.frame for ts in trj]
        except:
            raise AssertionError("with_statement not working for DCDReader")
        assert_equal(N, 98, err_msg="with_statement: DCDReader reads wrong number of frames")
        assert_array_equal(frames, np.arange(0, N), err_msg="with_statement: DCDReader does not read all frames")


class TestDCDReader(_TestDCD):
    def test_rewind_dcd(self):
        self.dcd.rewind()
        assert_equal(self.ts.frame, 0, "rewinding to frame 0")

    def test_next_dcd(self):
        self.dcd.rewind()
        self.dcd.next()
        assert_equal(self.ts.frame, 1, "loading frame 1")

    def test_jump_dcd(self):
        self.dcd[15]  # index is 0-based and frames are 0-based
        assert_equal(self.ts.frame, 15, "jumping to frame 15")

    def test_jump_lastframe_dcd(self):
        self.dcd[-1]
        assert_equal(self.ts.frame, 97, "indexing last frame with dcd[-1]")

    def test_slice_dcd(self):
        frames = [ts.frame for ts in self.dcd[5:17:3]]
        assert_equal(frames, [5, 8, 11, 14], "slicing dcd [5:17:3]")

    def test_reverse_dcd(self):
        frames = [ts.frame for ts in self.dcd[20:5:-1]]
        assert_equal(frames, range(20, 5, -1), "reversing dcd [20:5:-1]")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, 3341, "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 98, "wrong number of frames in dcd")

    def test_dt(self):
        assert_almost_equal(self.universe.trajectory.dt, 1.0, 4,
                            err_msg="wrong timestep dt")

    def test_totaltime(self):
        # test_totaltime(): need to reduce precision because dt is only precise
        # to ~4 decimals and accumulating the inaccuracy leads to even lower
        # precision in the totaltime (consequence of fixing Issue 64)
        assert_almost_equal(self.universe.trajectory.totaltime, 98.0, 3,
                            err_msg="wrong total length of AdK trajectory")

    def test_frame(self):
        self.dcd[15]  # index is 0-based and frames are 0-based
        assert_equal(self.universe.trajectory.frame, 15, "wrong frame number")

    def test_time(self):
        self.dcd[15]  # index is 0-based and frames are 0-based
        assert_almost_equal(self.universe.trajectory.time, 15.0, 5,
                            err_msg="wrong time of frame")

    def test_volume(self):
        assert_almost_equal(self.ts.volume, 0.0, 3,
                            err_msg="wrong volume for unitcell (no unitcell in DCD so this should be 0)")


class TestDCDWriter(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        ext = ".dcd"
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)
        self.Writer = MDAnalysis.coordinates.DCD.DCDWriter

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe
        del self.Writer

    @attr('issue')
    def test_write_trajectory(self):
        """Test writing DCD trajectories (Issue 50)"""
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms, dt=t.dt, step=t.skip_timestep)
        for ts in self.universe.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(PSF, self.outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, 3,
                                      err_msg="coordinate mismatch between original and written trajectory at frame "
                                              "%d (orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))

    def test_dt(self):
        DT = 5.0
        t = self.universe.trajectory
        with self.Writer(self.outfile, t.n_atoms, dt=DT) as W:  # set time step to 5 ps
            for ts in self.universe.trajectory:
                W.write_next_timestep(ts)

        uw = mda.Universe(PSF, self.outfile)
        assert_almost_equal(uw.trajectory.totaltime, uw.trajectory.n_frames * DT, 5)
        times = np.array([uw.trajectory.time for ts in uw.trajectory])
        frames = np.array([ts.frame for ts in uw.trajectory])
        assert_array_almost_equal(times, frames * DT, 5)

    def test_OtherWriter(self):
        t = self.universe.trajectory
        W = t.OtherWriter(self.outfile)
        for ts in self.universe.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(PSF, self.outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, 3,
                                      err_msg="coordinate mismatch between original and written trajectory at frame "
                                              "%d (orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))

    def test_single_frame(self):
        u = MDAnalysis.Universe(PSF, CRD)
        W = MDAnalysis.Writer(self.outfile, u.atoms.n_atoms)
        W.write(u.atoms)
        W.close()
        w = MDAnalysis.Universe(PSF, self.outfile)
        assert_equal(w.trajectory.n_frames, 1, "single frame trajectory has wrong number of frames")
        assert_almost_equal(w.atoms.coordinates(), u.atoms.coordinates(), 3,
                            err_msg="coordinates do not match")

    def test_with_statement(self):
        u = MDAnalysis.Universe(PSF, CRD)
        try:
            with MDAnalysis.Writer(self.outfile, u.atoms.n_atoms) as W:
                W.write(u.atoms)
        except AttributeError:  # misses __exit__
            raise AssertionError("DCDWriter: does not support with statement")
        w = MDAnalysis.Universe(PSF, self.outfile)
        assert_equal(w.trajectory.n_frames, 1, "with_statement: single frame trajectory has wrong number of frames")
        assert_almost_equal(w.atoms.coordinates(), u.atoms.coordinates(), 3,
                            err_msg="with_statement: coordinates do not match")


class TestDCDWriter_Issue59(TestCase):
    def setUp(self):
        """Generate input xtc."""
        self.u = MDAnalysis.Universe(PSF, DCD)
        fd, self.xtc = tempfile.mkstemp(suffix='.xtc')
        os.close(fd)
        wXTC = MDAnalysis.Writer(self.xtc, self.u.atoms.n_atoms)
        for ts in self.u.trajectory:
            wXTC.write(ts)
        wXTC.close()

    def tearDown(self):
        try:
            os.unlink(self.xtc)
        except OSError:
            pass
        try:
            os.unlink(self.dcd)
        except (AttributeError, OSError):
            pass
        del self.u

    @attr('issue')
    def test_issue59(self):
        """Test writing of XTC to DCD (Issue 59)"""
        xtc = MDAnalysis.Universe(PSF, self.xtc)
        fd, self.dcd = tempfile.mkstemp(suffix='.dcd')
        os.close(fd)
        wDCD = MDAnalysis.Writer(self.dcd, xtc.atoms.n_atoms)
        for ts in xtc.trajectory:
            wDCD.write(ts)
        wDCD.close()

        dcd = MDAnalysis.Universe(PSF, self.dcd)

        xtc.trajectory.rewind()
        dcd.trajectory.rewind()

        assert_array_almost_equal(xtc.atoms.coordinates(), dcd.atoms.coordinates(), 3,
                                  err_msg="XTC -> DCD: DCD coordinates are messed up (Issue 59)")

    def test_OtherWriter(self):
        dcd = self.u
        wXTC = dcd.trajectory.OtherWriter(self.xtc)
        for ts in dcd.trajectory:
            wXTC.write(ts)
        wXTC.close()

        xtc = MDAnalysis.Universe(PSF, self.xtc)
        xtc.trajectory.rewind()
        dcd.trajectory.rewind()

        assert_array_almost_equal(dcd.atoms.coordinates(), xtc.atoms.coordinates(), 2,
                                  err_msg="DCD -> XTC: coordinates are messed up (frame %d)" % dcd.trajectory.frame)
        xtc.trajectory[3]
        dcd.trajectory[3]
        assert_array_almost_equal(dcd.atoms.coordinates(), xtc.atoms.coordinates(), 2,
                                  err_msg="DCD -> XTC: coordinates are messed up (frame %d)" % dcd.trajectory.frame)

class RefCHARMMtriclinicDCD(object):
    topology = PSF_TRICLINIC
    trajectory = DCD_TRICLINIC
    # time(ps) A B C alpha beta gamma (length in Angstrome, angles in degrees)
    ref_dimensions = np.array([
            # [  0.     ,  35.     ,  35.     ,  35.     ,  90.     ,  60.     ,         45.     ], # dcd starts at t=1ps
            [  1.     ,  35.44604,  35.06156,  34.1585 ,  91.32802,  61.73521,         44.40703],
            [  2.     ,  34.65957,  34.22689,  33.09897,  90.56206,  61.79192,         44.14549],
            [  3.     ,  34.52772,  34.66422,  33.53881,  90.55859,  63.11228,         40.14044],
            [  4.     ,  34.43749,  33.38432,  34.02133,  88.82457,  64.98057,         36.77397],
            [  5.     ,  33.73129,  32.47752,  34.18961,  89.88102,  65.89032,         36.10921],
            [  6.     ,  33.78703,  31.90317,  34.98833,  90.03092,  66.12877,         35.07141],
            [  7.     ,  33.24708,  31.18271,  34.9654 ,  93.11122,  68.17743,         35.73643],
            [  8.     ,  32.92599,  30.31393,  34.99197,  93.89051,  69.3799 ,         33.48945],
            [  9.     ,  32.15295,  30.43056,  34.96157,  96.01416,  71.50115,         32.56111],
            [ 10.     ,  31.99748,  30.21518,  35.24292,  95.85821,  71.08429,         31.85939]])

class RefNAMDtriclinicDCD(object):
    topology = PSF_NAMD_TRICLINIC
    trajectory = DCD_NAMD_TRICLINIC
    # vmd topology trajectory
    # molinfo 0 get {a b c alpha beta gamma}
    # time(ps) A B C alpha beta gamma (length in Angstrome, angles in degrees)
    ref_dimensions = np.array([
            [1., 38.426594, 38.393101, 44.759800, 90.000000, 90.000000, 60.028915],
            ])

class _TestDCDReader_TriclinicUnitcell(TestCase):
    def setUp(self):
        self.u = MDAnalysis.Universe(self.topology, self.trajectory)
        fd, self.dcd = tempfile.mkstemp(suffix='.dcd')
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.dcd)
        except (AttributeError, OSError):
            pass
        del self.u

    @attr('issue')
    def test_read_triclinic(self):
        """test reading of triclinic unitcell (Issue 187) for NAMD or new CHARMM format (at least since c36b2)"""
        for ts, box in itertools.izip(self.u.trajectory, self.ref_dimensions[:, 1:]):
            assert_array_almost_equal(ts.dimensions, box, 4,
                                      err_msg="box dimensions A,B,C,alpha,beta,gamma not identical at frame {0}".format(ts.frame))
    @attr('issue')
    def test_write_triclinic(self):
        """test writing of triclinic unitcell (Issue 187) for NAMD or new CHARMM format (at least since c36b2)"""
        with self.u.trajectory.OtherWriter(self.dcd) as w:
            for ts in self.u.trajectory:
                w.write(ts)
        w = MDAnalysis.Universe(self.topology, self.dcd)
        for ts_orig, ts_copy in itertools.izip(self.u.trajectory, w.trajectory):
            assert_almost_equal(ts_orig.dimensions, ts_copy.dimensions, 4,
                                err_msg="DCD->DCD: unit cell dimensions wrong at frame {0}".format(ts_orig.frame))
        del w

class TestDCDReader_CHARMM_Unitcell(_TestDCDReader_TriclinicUnitcell, RefCHARMMtriclinicDCD):
    pass

class TestDCDReader_NAMD_Unitcell(_TestDCDReader_TriclinicUnitcell, RefNAMDtriclinicDCD):
    pass


class TestNCDF2DCD(TestCase):
    def setUp(self):
        self.u = MDAnalysis.Universe(PRMncdf, NCDF)
        # create the DCD
        fd, self.dcd = tempfile.mkstemp(suffix='.dcd')
        os.close(fd)
        DCD = MDAnalysis.Writer(self.dcd, n_atoms=self.u.atoms.n_atoms)
        for ts in self.u.trajectory:
            DCD.write(ts)
        DCD.close()
        self.w = MDAnalysis.Universe(PRMncdf, self.dcd)

    def tearDown(self):
        try:
            os.unlink(self.dcd)
        except (AttributeError, OSError):
            pass
        del self.u
        del self.w

    @attr('issue')
    def test_unitcell(self):
        """NCDFReader: Test that DCDWriter correctly writes the CHARMM unit cell"""
        for ts_orig, ts_copy in itertools.izip(self.u.trajectory, self.w.trajectory):
            assert_almost_equal(ts_orig.dimensions, ts_copy.dimensions, 3,
                                err_msg="NCDF->DCD: unit cell dimensions wrong at frame %d" % ts_orig.frame)

    def test_coordinates(self):
        for ts_orig, ts_copy in itertools.izip(self.u.trajectory, self.w.trajectory):
            assert_almost_equal(self.u.atoms.positions, self.w.atoms.positions, 3,
                                err_msg="NCDF->DCD: coordinates wrong at frame %d" % ts_orig.frame)


class TestDCDCorrel(_TestDCD):
    def setUp(self):
        # Note: setUp is executed for *every* test !
        super(TestDCDCorrel, self).setUp()
        import MDAnalysis.core.Timeseries as TS

        self.collection = TS.TimeseriesCollection()
        C = self.collection
        all = self.universe.atoms
        ca = self.universe.s4AKE.CA
        ca_termini = mda.core.AtomGroup.AtomGroup([ca[0], ca[-1]])
        # note that this is not quite phi... HN should be C of prec. residue
        phi151 = self.universe.select_atoms('resid 151').select_atoms('name HN', 'name N', 'name CA', 'name CB')
        C.addTimeseries(TS.Atom('v', ca_termini))  # 0
        C.addTimeseries(TS.Bond(ca_termini))  # 1
        C.addTimeseries(TS.Bond([ca[0], ca[-1]]))  # 2
        C.addTimeseries(TS.Angle(phi151[1:4]))  # 3
        C.addTimeseries(TS.Dihedral(phi151))  # 4
        C.addTimeseries(TS.Distance('r', ca_termini))  # 5
        C.addTimeseries(TS.CenterOfMass(ca))  # 6
        C.addTimeseries(TS.CenterOfGeometry(ca))  # 7
        C.addTimeseries(TS.CenterOfMass(all))  # 8
        C.addTimeseries(TS.CenterOfGeometry(all))  # 9
        # cannot test WaterDipole because there's no water in the test dcd
        C.compute(self.dcd)

    def tearDown(self):
        del self.collection
        super(TestDCDCorrel, self).tearDown()

    def test_correl(self):
        assert_equal(len(self.collection), 10, "Correl: len(collection)")

    def test_Atom(self):
        assert_equal(self.collection[0].shape, (2, 3, 98),
                     "Correl: Atom positions")

    def test_Bonds(self):
        C = self.collection
        assert_array_equal(C[1].__data__, C[2].__data__,
                           "Correl: Bonds with lists and AtomGroup")

    def test_Angle(self):
        C = self.collection
        avg_angle = 1.9111695972912988
        assert_almost_equal(C[3].__data__.mean(), avg_angle,
                            err_msg="Correl: average Angle")

    def test_Dihedral(self):
        C = self.collection
        avg_phi151 = 0.0088003870749735619
        assert_almost_equal(C[4].__data__.mean(), avg_phi151,
                            err_msg="Correl: average Dihedral")

    def test_scalarDistance(self):
        C = self.collection
        avg_dist = 9.7960210987736236
        assert_almost_equal(C[5].__data__.mean(), avg_dist,
                            err_msg="Correl: average scalar Distance")

    def test_CenterOfMass(self):
        C = self.collection
        avg_com_ca = np.array([0.0043688, -0.27812258, 0.0284051])
        avg_com_all = np.array([-0.10086529, -0.16357276, 0.12724672])
        assert_array_almost_equal(C[6].__data__.mean(axis=1), avg_com_ca,
                                  err_msg="Correl: average CA CenterOfMass")
        assert_almost_equal(C[8].__data__.mean(axis=1), avg_com_all,
                            err_msg="Correl: average all CenterOfMass")

    def test_CenterOfGeometry(self):
        C = self.collection
        avg_cog_all = np.array([-0.13554797, -0.20521885, 0.2118998])
        assert_almost_equal(C[9].__data__.mean(axis=1), avg_cog_all,
                            err_msg="Correl: average all CenterOfGeometry")

    def test_CA_COMeqCOG(self):
        C = self.collection
        assert_array_almost_equal(C[6].__data__, C[7].__data__,
                                  err_msg="Correl: CA CentreOfMass == CenterOfGeometry")

    def test_clear(self):
        C = self.collection
        C.clear()
        assert_equal(len(C), 0, "Correl: clear()")


# notes:
def compute_correl_references():
    universe = MDAnalysis.Universe(PSF, DCD)

    all = universe.atoms
    ca = universe.s4AKE.CA
    ca_termini = mda.core.AtomGroup.AtomGroup([ca[0], ca[-1]])
    phi151 = universe.select_atoms('resid 151').select_atoms('name HN', 'name N', 'name CA', 'name CB')

    C = MDAnalysis.collection
    C.clear()

    C.addTimeseries(TS.Atom('v', ca_termini))  # 0
    C.addTimeseries(TS.Bond(ca_termini))  # 1
    C.addTimeseries(TS.Bond([ca[0], ca[-1]]))  # 2
    C.addTimeseries(TS.Angle(phi151[1:4]))  # 3
    C.addTimeseries(TS.Dihedral(phi151))  # 4
    C.addTimeseries(TS.Distance('r', ca_termini))  # 5
    C.addTimeseries(TS.CenterOfMass(ca))  # 6
    C.addTimeseries(TS.CenterOfGeometry(ca))  # 7
    C.addTimeseries(TS.CenterOfMass(all))  # 8
    C.addTimeseries(TS.CenterOfGeometry(all))  # 9

    C.compute(universe.dcd)

    results = {
        "avg_angle": C[3].__data__.mean(),
        "avg_phi151": C[4].__data__.mean(),
        "avg_dist": C[5].__data__.mean(),
        "avg_com_ca": C[6].__data__.mean(axis=1),
        "avg_com_all": C[8].__data__.mean(axis=1),
        "avg_cog_all": C[9].__data__.mean(axis=1),
    }
    C.clear()
    return results


class TestChainReader(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, [DCD, CRD, DCD, CRD, DCD, CRD, CRD])
        self.trajectory = self.universe.trajectory
        self.prec = 3
        # dummy output DCD file
        fd, self.outfile = tempfile.mkstemp(suffix=".dcd")
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe

    def test_next_trajectory(self):
        self.trajectory.rewind()
        self.trajectory.next()
        assert_equal(self.trajectory.ts.frame, 1, "loading frame 2")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, 3341, "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 3 * 98 + 4, "wrong number of frames in chained dcd")

    def test_iteration(self):
        for ts in self.trajectory:
            pass  # just forward to last frame
        assert_equal(self.trajectory.n_frames - 1, ts.frame,
                     "iteration yielded wrong number of frames (%d), should be %d"
                     % (ts.frame, self.trajectory.n_frames))

    def test_jump_lastframe_trajectory(self):
        self.trajectory[-1]
        #print self.trajectory.ts, self.trajectory.ts.frame
        assert_equal(self.trajectory.ts.frame + 1, self.trajectory.n_frames, "indexing last frame with trajectory[-1]")

    def test_slice_trajectory(self):
        frames = [ts.frame for ts in self.trajectory[5:17:3]]
        assert_equal(frames, [5, 8, 11, 14], "slicing dcd [5:17:3]")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_frame_numbering(self):
        self.trajectory[98]  # index is 0-based and frames are 0-based
        assert_equal(self.universe.trajectory.frame, 98, "wrong frame number")

    def test_frame(self):
        self.trajectory[0]
        coord0 = self.universe.atoms.coordinates().copy()
        # forward to frame where we repeat original dcd again:
        # dcd:0..97 crd:98 dcd:99..196
        self.trajectory[99]
        assert_array_equal(self.universe.atoms.coordinates(), coord0,
                           "coordinates at frame 1 and 100 should be the same!")

    @knownfailure("time attribute not implemented for chained reader", ValueError)
    def test_time(self):
        self.trajectory[30]  # index is 0-based but frames are 1-based
        assert_almost_equal(self.universe.trajectory.time, 31.0, 5,
                            err_msg="wrong time of frame")

    @dec.slow
    def test_write_dcd(self):
        """test that ChainReader written dcd (containing crds) is correct (Issue 81)"""
        W = MDAnalysis.Writer(self.outfile, self.universe.atoms.n_atoms)
        for ts in self.universe.trajectory:
            W.write(self.universe)
        W.close()
        self.universe.trajectory.rewind()
        u = MDAnalysis.Universe(PSF, self.outfile)
        for (ts_orig, ts_new) in itertools.izip(self.universe.trajectory, u.trajectory):
            assert_almost_equal(ts_orig._pos, ts_new._pos, self.prec,
                                err_msg="Coordinates disagree at frame %d" % ts_orig.frame)


class TestChainReaderFormats(TestCase):
    """Test of ChainReader with explicit formats (Issue 76)."""

    @attr('issue')
    def test_set_all_format_tuples(self):
        universe = MDAnalysis.Universe(GRO, [(PDB, 'pdb'), (XTC, 'xtc'), (TRR, 'trr')])
        assert_equal(universe.trajectory.n_frames, 21)

    @attr('issue')
    def test_set_one_format_tuple(self):
        universe = MDAnalysis.Universe(PSF, [(PDB_small, 'pdb'), DCD])
        assert_equal(universe.trajectory.n_frames, 99)

    @attr('issue')
    def test_set_all_formats(self):
        universe = MDAnalysis.Universe(PSF, [PDB_small, PDB_closed], format='pdb')
        assert_equal(universe.trajectory.n_frames, 2)


class TestTRRReader_Sub(TestCase):
    def setUp(self):
        """
        grab values from selected atoms from full solvated traj,
        later compare to using 'sub'
        """
        usol = mda.Universe(PDB_sub_sol, TRR_sub_sol)
        atoms = usol.select_atoms("not resname SOL")
        self.pos = atoms.positions
        self.vel = atoms.velocities
        self.force = atoms.forces
        self.sub = atoms.indices
        # universe from un-solvated protein
        self.udry = mda.Universe(PDB_sub_dry)

    def test_load_new_raises_ValueError(self):
        # should fail if we load universe with a trajectory with different
        # number of atoms when NOT using sub, same as before.
        def load_new_without_sub():
            self.udry.load_new(TRR_sub_sol)

        assert_raises(ValueError, load_new_without_sub)

    def test_sub_coordinates(self):
        """
        load solvated trajectory into universe with unsolvated protein.
        """
        self.udry.load_new(TRR_sub_sol, sub=self.sub)
        assert_array_almost_equal(self.pos, self.udry.atoms.positions,
                                  err_msg="positions differ")
        assert_array_almost_equal(self.vel, self.udry.atoms.velocities,
                                  err_msg="positions differ")
        assert_array_almost_equal(self.force, self.udry.atoms.forces,
                                  err_msg="positions differ")


class _GromacsReader(TestCase):
    # This base class assumes same lengths and dt for XTC and TRR test cases!
    filename = None
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.], dtype=np.float32)
    ref_volume = 362270.0  # computed with Gromacs: 362.26999999999998 nm**3 * 1000 A**3/nm**3

    def setUp(self):
        # loading from GRO is 4x faster than the PDB reader
        self.universe = mda.Universe(GRO, self.filename, convert_units=True)
        self.trajectory = self.universe.trajectory
        self.prec = 3
        self.ts = self.universe.coord
        # dummy output file
        ext = os.path.splitext(self.filename)[1]
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe

    @dec.slow
    def test_flag_convert_lengths(self):
        assert_equal(mda.core.flags['convert_lengths'], True,
                     "MDAnalysis.core.flags['convert_lengths'] should be True by default")

    @dec.slow
    def test_rewind_xdrtrj(self):
        self.trajectory.rewind()
        assert_equal(self.ts.frame, 0, "rewinding to frame 1")

    @dec.slow
    def test_next_xdrtrj(self):
        self.trajectory.rewind()
        self.trajectory.next()
        assert_equal(self.ts.frame, 1, "loading frame 1")

    @dec.slow
    def test_jump_xdrtrj(self):
        self.trajectory[4]  # index is 0-based and frames are 0-based
        assert_equal(self.ts.frame, 4, "jumping to frame 4")

    @dec.slow
    def test_jump_lastframe_xdrtrj(self):
        self.trajectory[-1]
        assert_equal(self.ts.frame, 9, "indexing last frame with trajectory[-1]")

    @dec.slow
    def test_slice_xdrtrj(self):
        frames = [ts.frame for ts in self.trajectory[2:9:3]]
        assert_equal(frames, [2, 5, 8], "slicing xdrtrj [2:9:3]")

    @dec.slow
    def test_reverse_xdrtrj(self):
        frames = [ts.frame for ts in self.trajectory[::-1]]
        assert_equal(frames, range(9, -1, -1), "slicing xdrtrj [::-1]")

    @dec.slow
    def test_coordinates(self):
        ca_nm = np.array([[6.043369675, 7.385184479, 1.381425762]], dtype=np.float32)
        # coordinates in the base unit (needed for True)
        ca_Angstrom = ca_nm * 10.0
        U = self.universe
        T = U.trajectory
        T.rewind()
        T.next()
        T.next()
        assert_equal(self.ts.frame, 2, "failed to step to frame 3")
        ca = U.select_atoms('name CA and resid 122')
        # low precision match (2 decimals in A, 3 in nm) because the above are the trr coords
        assert_array_almost_equal(ca.coordinates(), ca_Angstrom, 2,
                                  err_msg="coords of Ca of resid 122 do not match for frame 3")

    @dec.slow
    @attr('issue')
    def test_unitcell(self):
        """Test that xtc/trr unitcell is read correctly (Issue 34)"""
        self.universe.trajectory.rewind()
        uc = self.ts.dimensions
        assert_array_almost_equal(uc, self.ref_unitcell, self.prec,
                                  err_msg="unit cell dimensions (rhombic dodecahedron)")

    @dec.slow
    def test_volume(self):
        # need to reduce precision for test (nm**3 <--> A**3)
        self.universe.trajectory.rewind()
        vol = self.ts.volume
        assert_array_almost_equal(vol, self.ref_volume, 0, err_msg="unit cell volume (rhombic dodecahedron)")

    @dec.slow
    def test_dt(self):
        assert_almost_equal(self.universe.trajectory.dt, 100.0, 4,
                            err_msg="wrong timestep dt")

    @dec.slow
    def test_totaltime(self):
        # test_totaltime(): need to reduce precision because dt is only precise
        # to ~4 decimals and accumulating the inaccuracy leads to even lower
        # precision in the totaltime (consequence of fixing Issue 64)
        assert_almost_equal(self.universe.trajectory.totaltime, 1000.0, 3,
                            err_msg="wrong total length of trajectory")

    @dec.slow
    def test_frame(self):
        self.trajectory[4]  # index is 0-based and frames are 0-based
        assert_equal(self.universe.trajectory.frame, 4, "wrong frame number")

    @dec.slow
    def test_time(self):
        self.trajectory[4]
        assert_almost_equal(self.universe.trajectory.time, 400.0, 3,
                            err_msg="wrong time of frame")

    @dec.slow
    def test_get_Writer(self):
        W = self.universe.trajectory.Writer(self.outfile)
        assert_equal(self.universe.trajectory.format, W.format)
        assert_equal(self.universe.atoms.n_atoms, W.n_atoms)
        W.close()

    @dec.slow
    def test_Writer(self):
        W = self.universe.trajectory.Writer(self.outfile)
        W.write(self.universe.atoms)
        self.universe.trajectory.next()
        W.write(self.universe.atoms)
        W.close()
        self.universe.trajectory.rewind()
        u = MDAnalysis.Universe(GRO, self.outfile)
        assert_equal(u.trajectory.n_frames, 2)
        # prec = 6: TRR test fails; here I am generous and take self.prec = 3...
        assert_almost_equal(u.atoms.coordinates(), self.universe.atoms.coordinates(), self.prec)

    @dec.slow
    def test_EOFraisesIOErrorEIO(self):
        def go_beyond_EOF():
            self.universe.trajectory[-1]
            self.universe.trajectory.next()

        assert_raises(IOError, go_beyond_EOF)
        try:
            go_beyond_EOF()
        except IOError as err:
            assert_equal(err.errno, errno.EIO, "IOError produces wrong error code")


class TestXTCReader(_GromacsReader):
    filename = XTC


class TestXTCReaderClass(TestCase):
    def test_with_statement(self):
        from MDAnalysis.coordinates.XTC import XTCReader

        try:
            with XTCReader(XTC) as trj:
                N = trj.n_frames
                frames = [ts.frame for ts in trj]
        except:
            raise AssertionError("with_statement not working for XTCReader")
        assert_equal(N, 10, err_msg="with_statement: XTCReader reads wrong number of frames")
        assert_array_equal(frames, np.arange(0, N), err_msg="with_statement: XTCReader does not read all frames")


class TestTRRReader(_GromacsReader):
    filename = TRR

    @dec.slow
    def test_velocities(self):
        # frame 0, v in nm/ps
        # from gmxdump -f MDAnalysisTests/data/adk_oplsaa.trr
        #      v[47675]={-7.86469e-01,  1.57479e+00,  2.79722e-01}
        #      v[47676]={ 2.70593e-08,  1.08052e-06,  6.97028e-07}
        v_native = np.array(
            [
                [-7.86469e-01, 1.57479e+00, 2.79722e-01],
                [2.70593e-08, 1.08052e-06, 6.97028e-07]
            ], dtype=np.float32)

        # velocities in the MDA base unit A/ps (needed for True)
        v_base = v_native * 10.0
        self.universe.trajectory.rewind()
        assert_equal(self.ts.frame, 0, "failed to read frame 1")

        assert_array_almost_equal(self.universe.trajectory.ts._velocities[[47675, 47676]], v_base, self.prec,
                                  err_msg="ts._velocities for indices 47675,47676 do not match known values")

        assert_array_almost_equal(self.universe.atoms.velocities[[47675, 47676]], v_base, self.prec,
                                  err_msg="velocities for indices 47675,47676 do not match known values")

        for index, v_known in zip([47675, 47676], v_base):
            assert_array_almost_equal(self.universe.atoms[index].velocity, v_known, self.prec,
                                      err_msg="atom[%d].velocity does not match known values" % index)

class _XDRNoConversion(TestCase):
    filename = None

    def setUp(self):
        self.universe = mda.Universe(PDB, self.filename, convert_units=False)
        self.ts = self.universe.trajectory.ts

    def tearDown(self):
        del self.universe
        del self.ts

    @dec.slow
    def test_coordinates(self):
        # note: these are the native coordinates in nm
        ca_nm = np.array([[6.043369675, 7.385184479, 1.381425762]], dtype=np.float32)
        U = self.universe
        T = U.trajectory
        T.rewind()
        T.next()
        T.next()
        assert_equal(self.ts.frame, 2, "failed to step to frame 3")
        ca = U.select_atoms('name CA and resid 122')
        # low precision match because we also look at the trr: only 3 decimals in nm in xtc!
        assert_array_almost_equal(ca.coordinates(), ca_nm, 3,
                                  err_msg="native coords of Ca of resid 122 do not match for frame 3 "
                                          "with convert_units=False")


class TestXTCNoConversion(_XDRNoConversion):
    filename = XTC


class TestTRRNoConversion(_XDRNoConversion):
    filename = TRR


class _GromacsWriter(TestCase):
    infilename = None  # XTC or TRR
    Writers = {
        '.trr': MDAnalysis.coordinates.TRR.TRRWriter,
        '.xtc': MDAnalysis.coordinates.XTC.XTCWriter,
    }

    def setUp(self):
        self.universe = mda.Universe(GRO, self.infilename)
        ext = os.path.splitext(self.infilename)[1]
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)
        self.Writer = self.Writers[ext]

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe
        del self.Writer

    @dec.slow
    @attr('issue')
    def test_write_trajectory(self):
        """Test writing Gromacs trajectories (Issue 38)"""
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms, dt=t.dt)
        for ts in self.universe.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(GRO, self.outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, 3,
                                      err_msg="coordinate mismatch between original and written trajectory at frame "
                                              "%d (orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))

    @dec.slow
    def test_timestep_not_modified_by_writer(self):
        trj = self.universe.trajectory
        ts = trj.ts

        trj[-1]  # last timestep (so that time != 0)
        x = ts._pos.copy()
        time = ts.time

        W = self.Writer(self.outfile, trj.n_atoms, dt=trj.dt)
        trj[-1]  # last timestep (so that time != 0) (say it again, just in case...)
        W.write_next_timestep(ts)
        W.close()

        assert_equal(ts._pos, x, err_msg="Positions in Timestep were modified by writer.")
        assert_equal(ts.time, time, err_msg="Time in Timestep was modified by writer.")


class TestXTCWriter(_GromacsWriter):
    infilename = XTC


class TestTRRWriter(_GromacsWriter):
    infilename = TRR

    def test_velocities(self):
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms, dt=t.dt)
        for ts in self.universe.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(GRO, self.outfile)

        # check that the velocities are identical for each time step
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._velocities, orig_ts._velocities, 3,
                                      err_msg="velocities mismatch between original and written trajectory at frame "
                                              "%d (orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))

    def test_gaps(self):
        """Tests the writing and reading back of TRRs with gaps in any of the coordinates/velocities properties."""
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms, dt=t.dt)
        for ts in self.universe.trajectory:
            # Inset some gaps in the properties: coords every 4 steps, vels every 2.
            if not ts.frame % 4:
                ts.has_positions = False
            if not ts.frame % 2:
                ts.has_velocities = False
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(GRO, self.outfile)

        # check that the velocities are identical for each time step, except for the gaps
        # (that we must make sure to raise exceptions on).
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            if ts.frame % 4:
                assert_array_almost_equal(written_ts.positions, orig_ts.positions, 3,
                                          err_msg="coordinates mismatch between original and written trajectory at "
                                                  "frame {0} (orig) vs {1} (written)".format(
                                                      orig_ts.frame, written_ts.frame))
            else:
                assert_raises(NoDataError, getattr, written_ts, 'positions')

            if ts.frame % 2:
                assert_array_almost_equal(written_ts.velocities, orig_ts.velocities, 3,
                                          err_msg="velocities mismatch between original and written trajectory at "
                                                  "frame {0} (orig) vs {1} (written)".format(
                                                      orig_ts.frame, written_ts.frame))
            else:
                assert_raises(NoDataError, getattr, written_ts, 'velocities')


class _GromacsWriterIssue101(TestCase):
    Writers = {
        '.trr': MDAnalysis.coordinates.TRR.TRRWriter,
        '.xtc': MDAnalysis.coordinates.XTC.XTCWriter,
    }
    ext = None  # set to '.xtc' or '.trr'
    prec = 3

    def setUp(self):
        fd, self.outfile = tempfile.mkstemp(suffix=self.ext)
        os.close(fd)
        self.Writer = self.Writers[self.ext]

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.Writer

    @dec.slow
    @attr('issue')
    def test_single_frame_GRO(self):
        self._single_frame(GRO)

    @dec.slow
    @attr('issue')
    def test_single_frame_PDB(self):
        self._single_frame(PDB)

    @attr('issue')
    def test_single_frame_CRD(self):
        self._single_frame(CRD)

    def _single_frame(self, filename):
        u = MDAnalysis.Universe(filename)
        with self.Writer(self.outfile, u.atoms.n_atoms) as W:
            W.write(u.atoms)
        w = MDAnalysis.Universe(filename, self.outfile)
        assert_equal(w.trajectory.n_frames, 1, "single frame trajectory has wrong number of frames")
        assert_almost_equal(w.atoms.coordinates(), u.atoms.coordinates(), self.prec,
                            err_msg="coordinates do not match for %r" % filename)


class TestXTCWriterSingleFrame(_GromacsWriterIssue101):
    ext = ".xtc"
    prec = 2


class TestTRRWriterSingleFrame(_GromacsWriterIssue101):
    ext = ".trr"


class _GromacsWriterIssue117(TestCase):
    """Issue 117: Cannot write XTC or TRR from AMBER NCDF"""
    ext = None
    prec = 5

    def setUp(self):
        self.universe = mda.Universe(PRMncdf, NCDF)
        fd, self.outfile = tempfile.mkstemp(suffix=self.ext)
        os.close(fd)
        self.Writer = MDAnalysis.Writer(self.outfile, n_atoms=self.universe.atoms.n_atoms)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe
        del self.Writer

    @attr('issue')
    def test_write_trajectory(self):
        """Test writing Gromacs trajectories from AMBER NCDF (Issue 117)"""
        t = self.universe.trajectory
        for ts in self.universe.trajectory:
            self.Writer.write_next_timestep(ts)
        self.Writer.close()

        uw = MDAnalysis.Universe(PRMncdf, self.outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between original and written trajectory at frame "
                                              "%d (orig) vs %d (written)" % (
                                      orig_ts.frame, written_ts.frame))


class TestXTCWriterIssue117(_GromacsWriterIssue117):
    ext = ".xtc"
    prec = 2


class TestTRRWriterIssue117(_GromacsWriterIssue117):
    ext = ".trr"


@attr('issue')
def test_triclinic_box():
    """Test coordinates.core.triclinic_box() (Issue 61)"""
    unitcell = np.array([80.017, 55, 100.11, 60.00, 30.50, 90.00])
    box = MDAnalysis.coordinates.core.triclinic_vectors(unitcell)
    new_unitcell = MDAnalysis.coordinates.core.triclinic_box(box[0], box[1], box[2])
    assert_array_almost_equal(new_unitcell, unitcell, 3,
                              err_msg="unitcell round-trip connversion failed (Issue 61)")


class RefTRZ(object):
    #    ref_coordinates = {}
    #    ref_distances = {'endtoend': }
    ref_n_atoms = 8184
    ref_dimensions = np.array([55.422830581665039, 55.422830581665039, 55.422830581665039, 90., 90., 90.],
                              dtype=np.float32)
    ref_volume = 170241.762765
    ref_n_frames = 6
    ref_coordinates = np.array([72.3163681, -130.31130981, 19.97969055], dtype=np.float32)
    ref_velocities = np.array([[14.83297443, 18.02611542, 6.07733774]], dtype=np.float32)
    ref_delta = 0.001
    ref_time = 0.01


class TestTRZReader(TestCase, RefTRZ):
    def setUp(self):
        self.universe = mda.Universe(TRZ_psf, TRZ)
        self.trz = self.universe.trajectory
        self.ts = self.universe.trajectory.ts
        self.prec = 3

    def tearDown(self):
        del self.universe
        del self.trz
        del self.ts

    def test_load_trz(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms, "load Universe from PSF and TRZ")

    def test_next_trz(self):
        assert_equal(self.ts.frame, 0, "starts at first frame")
        self.trz.next()
        assert_equal(self.ts.frame, 1, "next returns frame index 1")

    def test_rewind_trz(self):
        # move to different frame and rewind to get first frame back
        self.trz[2]
        self.trz.rewind()
        assert_equal(self.ts.frame, 0, "rewinding to frame 1")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, self.ref_n_frames, "wrong number of frames in trz")

    def test_seeking(self):
        self.universe.trajectory[3]
        assert_equal(self.ts.frame, 3, "loading frame 3")

        orig = self.universe.atoms[0:3].positions.copy()

        self.universe.trajectory[4]
        assert_equal(self.ts.frame, 4, "loading frame 4")
        self.universe.trajectory[3]

        assert_almost_equal(self.universe.atoms[0:3].positions, orig, self.prec)

        self.universe.trajectory[0]
        assert_equal(self.ts.frame, 0, "loading frame 0")
        self.universe.trajectory[3]

        assert_almost_equal(self.universe.atoms[0:3].positions, orig, self.prec)

    def test_volume(self):
        assert_almost_equal(self.ts.volume, self.ref_volume, 1,
                            "wrong volume for trz")  # Lower precision here because errors seem to accumulate and
                            # throw this off (is rounded value**3)

    def test_unitcell(self):
        assert_almost_equal(self.ts.dimensions, self.ref_dimensions, self.prec, "wrong dimensions for trz")

    def test_coordinates(self):
        fortytwo = self.universe.atoms[41]  # 41 because is 0 based
        assert_almost_equal(fortytwo.pos, self.ref_coordinates, self.prec, "wrong coordinates in trz")

    def test_velocities(self):
        fortytwo = self.universe.select_atoms('bynum 42')
        assert_almost_equal(fortytwo.velocities, self.ref_velocities, self.prec, "wrong velocities in trz")

    def test_delta(self):
        assert_almost_equal(self.trz.delta, self.ref_delta, self.prec, "wrong time delta in trz")

    def test_time(self):
        assert_almost_equal(self.trz.time, self.ref_time, self.prec, "wrong time value in trz")

    def test_get_writer(self):
        fd, self.outfile = tempfile.mkstemp(suffix='.trz')
        os.close(fd)
        W = self.trz.Writer(self.outfile)
        assert_equal(isinstance(W, MDAnalysis.coordinates.TRZ.TRZWriter), True)
        assert_equal(W.n_atoms, self.trz.n_atoms)
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

    def test_get_writer_2(self):
        fd, self.outfile = tempfile.mkstemp(suffix='.trz')
        os.close(fd)
        W = self.trz.Writer(self.outfile, n_atoms=100)
        assert_equal(isinstance(W, MDAnalysis.coordinates.TRZ.TRZWriter), True)
        assert_equal(W.n_atoms, 100)
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

class TestTRZWriter(TestCase, RefTRZ):
    def setUp(self):
        self.universe = mda.Universe(TRZ_psf, TRZ)
        self.prec = 3
        fd, self.outfile = tempfile.mkstemp(suffix='.trz')
        os.close(fd)
        self.Writer = MDAnalysis.coordinates.TRZ.TRZWriter

    def tearDown(self):
        del self.universe
        del self.prec
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.Writer

    def test_write_trajectory(self):
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms)
        self._copy_traj(W)

    def _copy_traj(self, writer):
        for ts in self.universe.trajectory:
            writer.write_next_timestep(ts)
        writer.close()

        uw = mda.Universe(TRZ_psf, self.outfile)

        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(orig_ts._pos, written_ts._pos, self.prec,
                                      err_msg="Coordinate mismatch between orig and written at frame %d" %
                                              orig_ts.frame)
            assert_array_almost_equal(orig_ts._velocities, written_ts._velocities, self.prec,
                                      err_msg="Coordinate mismatch between orig and written at frame %d" %
                                              orig_ts.frame)
            assert_array_almost_equal(orig_ts._unitcell, written_ts._unitcell, self.prec,
                                      err_msg="Unitcell mismatch between orig and written at frame %d" % orig_ts.frame)
            for att in orig_ts.data:
                assert_array_almost_equal(orig_ts.data[att], written_ts.data[att], self.prec,
                                          err_msg="TS equal failed for %s" % att)

class TestTRZWriter2(object):
    def setUp(self):
        self.u = mda.Universe(two_water_gro)

    def tearDown(self):
        del self.u
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

    def test_writer_trz_from_other(self):
        fd, self.outfile = tempfile.mkstemp(suffix='.trz')
        os.close(fd)
        W = MDAnalysis.coordinates.TRZ.TRZWriter(self.outfile, n_atoms=len(self.u.atoms))

        W.write(self.u.trajectory.ts)
        W.close()

        u2 = mda.Universe(two_water_gro, self.outfile)

        assert_array_almost_equal(self.u.atoms.positions, u2.atoms.positions, 3)


class TestWrite_Partial_Timestep(TestCase):
    """Test writing a partial timestep made by passing only an atomgroup to Writer. (Issue 163)

    The contents of the AtomGroup.ts are checked in test_atomgroup, this test just checks that Writer
    is receiving this information properly.
    """

    def setUp(self):
        self.universe = mda.Universe(TRZ_psf, TRZ)
        self.ag = self.universe.select_atoms('name N')
        self.prec = 3
        fd, self.outfile = tempfile.mkstemp(suffix='.pdb')
        os.close(fd)
        self.Writer = MDAnalysis.Writer(self.outfile, n_atoms=len(self.ag))

    def tearDown(self):
        del self.universe
        del self.ag
        del self.prec
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.Writer

    def test_write_trajectory(self):
        self.Writer.write(self.ag)
        self.Writer.close()

        u_ag = mda.Universe(self.outfile)

        assert_array_almost_equal(self.ag.coordinates(), u_ag.atoms.coordinates(), self.prec,
                                  err_msg="Writing AtomGroup timestep failed.")


class RefLAMMPSData(object):
    filename = LAMMPSdata
    n_atoms = 18360
    pos_atom1 = np.array([11.89985657, 48.4455719, 19.09719849], dtype=np.float32)
    vel_atom1 = np.array([-5.667593, 7.91380978, -3.00779533], dtype=np.float32)
    dimensions = np.array([55.42282867, 55.42282867, 55.42282867, 90., 90., 90.],
                          dtype=np.float32)


class RefLAMMPSDataMini(object):
    filename = LAMMPSdata_mini
    n_atoms = 1
    pos_atom1 = np.array([11.89985657, 48.4455719, 19.09719849], dtype=np.float32)
    vel_atom1 = np.array([-5.667593, 7.91380978, -3.00779533], dtype=np.float32)
    dimensions = np.array([60., 50., 30., 90., 90., 90.], dtype=np.float32)


def test_datareader_VE():
    from MDAnalysis.coordinates.LAMMPS import DATAReader
    assert_raises(ValueError, DATAReader, 'filename')


class _TestLammpsData_Coords(TestCase):
    """Tests using a .data file for loading single frame.

    All topology loading from MDAnalysisTests.data is done in test_topology
    """

    def setUp(self):
        self.u = MDAnalysis.Universe(self.filename)

    def tearDown(self):
        del self.u

    def test_n_atoms(self):
        assert_equal(self.u.atoms.n_atoms, self.n_atoms)

    def test_coords(self):
        assert_equal(self.u.atoms[0].pos, self.pos_atom1)

    def test_velos(self):
        assert_equal(self.u.atoms[0].velocity, self.vel_atom1)

    def test_dimensions(self):
        assert_equal(self.u.dimensions, self.dimensions)

    def test_singleframe(self):
        assert_raises(IOError, self.u.trajectory.next)

    def test_seek(self):
        assert_raises(IndexError, self.u.trajectory.__getitem__, 1)

    def test_seek_2(self):
        ts = self.u.trajectory[0]
        assert_equal(type(ts), MDAnalysis.coordinates.base.Timestep)

    def test_iter(self):
        # Check that iterating works, but only gives a single frame
        assert_equal(len(list(iter(self.u.trajectory))), 1)


class TestLammpsData_Coords(_TestLammpsData_Coords, RefLAMMPSData):
    pass


class TestLammpsDataMini_Coords(_TestLammpsData_Coords, RefLAMMPSDataMini):
    pass


class _DLPConfig(object):
    def setUp(self):
        self.r = MDAnalysis.coordinates.DLPoly.ConfigReader
        rd = self.rd = self.r(self.f)
        self.ts = rd.ts

    def tearDown(self):
        del self.r
        del self.rd
        del self.ts

    def test_read_unitcell(self):
        ref = np.array([[18.6960000000, 0.0000000000, 0.0000000000],
                        [0.0000000000, 18.6960000000, 0.0000000000],
                        [0.0000000000, 0.0000000000, 18.6960000000]])
        assert_allclose(self.ts._unitcell, ref)

    def test_positions(self):
        ref = np.array([-7.608595309, -7.897790000, -7.892053559])
        assert_allclose(self.ts._pos[0], ref)

    def test_velocities(self):
        ref = np.array([1.056610291, -1.218664448, 3.345828610])
        assert_allclose(self.ts._velocities[0], ref)

    def test_forces(self):
        ref = np.array([-1979.558687, 739.7961625, 1027.996603])
        assert_allclose(self.ts._forces[0], ref)


class TestConfigReader(_DLPConfig):
    f = DLP_CONFIG

    def test_read(self):
        assert self.rd.title == "DL_POLY: Potassium Chloride Test Case"


class TestConfigOrder(_DLPConfig):
    f = DLP_CONFIG_order

class TestConfigMinimal(_DLPConfig):
    f = DLP_CONFIG_minimal

    def test_read_unitcell(self):
        pass

    def test_velocities(self):
        assert_raises(AttributeError, getattr, self.ts, "_velocities")

    def test_forces(self):
        assert_raises(AttributeError, getattr, self.ts, "_forces")


class _DLPConfig2(object):
    def setUp(self):
        self.u = mda.Universe(self.f, format='CONFIG')

    def tearDown(self):
        del self.u

    def test_names(self):
        ref = ['C', 'B', 'A']
        assert_equal([a.name for a in self.u.atoms], ref)

    def test_pos(self):
        ref = np.array([-7.821414265, -4.635443539, -4.732164540])
        assert_allclose(self.u.atoms[2].pos, ref)

    def test_vel(self):
        ref = np.array([2.637614561, 0.5778767520E-01, -1.704765568])
        assert_allclose(self.u.atoms[2].velocity, ref)

    def test_for(self):
        ref = np.array([150.3309776, -812.6932914, 1429.413120])
        assert_allclose(self.u.atoms[2].force, ref)

    def test_number(self):
        ref = [0, 1, 2]
        assert_equal([a.index for a in self.u.atoms], ref)


class TestConfigReader2(_DLPConfig2):
    f = DLP_CONFIG_order

class TestConfigReaderMinimal2(_DLPConfig2):
    f = DLP_CONFIG_minimal

    def test_vel(self):
        pass

    def test_for(self):
        pass


class _DLHistory(object):
    def setUp(self):
        self.u = mda.Universe(self.f, format='HISTORY')

    def tearDown(self):
        self.u.trajectory.close()
        del self.u

    def test_len(self):
        assert_equal(len(self.u.trajectory), 3)
        assert_equal([ts.frame for ts in self.u.trajectory], [1, 2, 3])

    def test_getting(self):
        ts = self.u.trajectory[1]
        assert_equal(ts.frame, 2)

    def test_slicing(self):
        nums = [ts.frame for ts in self.u.trajectory[::2]]
        assert_equal(nums, [1, 3])

    def test_slicing_2(self):
        nums = [ts.frame for ts in self.u.trajectory[1::-2]]
        assert_equal(nums, [2])

    def test_position(self):
        ref = np.array([[-7.595541651, -7.898808509, -7.861763110],
                        [-7.019565641, -7.264933320, -7.045213551],
                        [-6.787470785, -6.912685099, -6.922156843]])
        for ts, r in itertools.izip(self.u.trajectory, ref):
            assert_allclose(self.u.atoms[0].pos, r)

    def test_velocity(self):
        ref = np.array([[1.109901682, -1.500264697, 4.752251711],
                        [-1.398479696, 2.091141311, 1.957430003],
                        [0.2570827995, -0.7146878577, -3.547444215]])
        for ts, r in itertools.izip(self.u.trajectory, ref):
            assert_allclose(self.u.atoms[0].velocity, r)

    def test_force(self):
        ref = np.array([[-2621.386432, 1579.334443, 1041.103241],
                        [-1472.262341, 2450.379615, -8149.916193],
                        [2471.802059, -3828.467296, 3596.679326]])
        for ts, r in itertools.izip(self.u.trajectory, ref):
            assert_allclose(self.u.atoms[0].force, r)

    def test_unitcell(self):
        ref1 = np.array([[18.6796195135, 0.0000058913, -0.0000139999],
                        [0.0000058913, 18.6794658887, -0.0000016255],
                        [-0.0000139999, -0.0000016255, 18.6797229304]])
        ref2 = np.array([[17.2277221163, -0.0044216126, -0.0003229237],
                         [-0.0044205826, 17.2124253987, 0.0019439244],
                         [-0.0003226531, 0.0019445826, 17.2416976104]])
        ref3 = np.array([[16.5435673205, -0.0108424742, 0.0014935464],
                         [-0.0108333201, 16.5270298891, 0.0011094612],
                         [0.0014948739, 0.0011058349, 16.5725517831]])
        for ts, r in itertools.izip(self.u.trajectory, [ref1, ref2, ref3]):
            assert_allclose(ts._unitcell, r)


class TestDLPolyHistory(_DLHistory):
    f = DLP_HISTORY

class TestDLPolyHistoryOrder(_DLHistory):
    f = DLP_HISTORY_order

class TestDLPolyHistoryMinimal(_DLHistory):
    f = DLP_HISTORY_minimal

    def test_velocity(self):
        assert_raises(NoDataError, getattr, self.u.atoms[0], 'velocity')

    def test_force(self):
        assert_raises(NoDataError, getattr, self.u.atoms[0], 'force')

    def test_unitcell(self):
        pass


class TestIncompletePDB(object):
    """Tests for Issue #396

    Reads an incomplete (but still intelligible) PDB file
    """
    def setUp(self):
        self.u = MDAnalysis.Universe(INC_PDB)

    def tearDown(self):
        del self.u

    def test_natoms(self):
        assert len(self.u.atoms) == 3

    def test_coords(self):
        assert_array_almost_equal(self.u.atoms.positions,
                                  np.array([[111.2519989, 98.3730011, 98.18699646],
                                            [111.20300293, 101.74199677, 96.43000031],
                                            [107.60700226, 102.96800232, 96.31600189]],
                                           dtype=np.float32))

    def test_dims(self):
        assert_array_almost_equal(self.u.dimensions,
                                  np.array([ 216.48899841, 216.48899841, 216.48899841,
                                             90., 90., 90.], dtype=np.float32))

    def test_names(self):
        assert all(self.u.atoms.names == 'CA')

    def test_residues(self):
        assert len(self.u.residues) == 3

    def test_resnames(self):
        assert len(self.u.atoms.resnames) == 3
        assert 'VAL' in self.u.atoms.resnames
        assert 'LYS' in self.u.atoms.resnames
        assert 'PHE' in self.u.atoms.resnames

    def test_reading_trajectory(self):
        for ts in self.u.trajectory:
            pass

    def test_occupancy(self):
        occupancies = self.u.atoms.occupancies
        assert_array_almost_equal(occupancies,
                                  np.ones(len(occupancies)))

    def test_set_occupancy(self):
        for atom in self.u.atoms:
            atom.occupancy = 0
        assert_almost_equal(self.u.atoms.occupancies,
                            np.zeros(self.u.atoms.n_atoms))

    def test_set_occupancies(self):
        self.u.atoms.occupancies = 0.0
        assert_almost_equal(self.u.atoms.occupancies,
                            np.zeros(self.u.atoms.n_atoms))
