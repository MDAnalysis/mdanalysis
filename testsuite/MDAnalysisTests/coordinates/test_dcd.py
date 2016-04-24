import MDAnalysis as mda
import numpy as np
import os
from six.moves import zip, range

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_array_equal, assert_raises,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, dec)
from unittest import TestCase

from MDAnalysisTests.datafiles import (DCD, PSF, DCD_empty, CRD, PRMncdf, NCDF)
from MDAnalysisTests.coordinates.reference import (RefCHARMMtriclinicDCD,
                                                   RefNAMDtriclinicDCD)
from MDAnalysisTests.coordinates.base import BaseTimestepTest
from MDAnalysisTests import module_not_found, tempdir


@attr('issue')
def TestDCD_Issue32():
    """Test for Issue 32: 0-size dcds lead to a segfault: now caught with
    IOError"""
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
        assert_equal(
            N,
            98,
            err_msg="with_statement: DCDReader reads wrong number of frames")
        assert_array_equal(
            frames,
            np.arange(0, N),
            err_msg="with_statement: DCDReader does not read all frames")


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

    def test_list_trajectory(self):
        frames = [ts.frame for ts in self.dcd[[0, 3, 4, 5]]]
        assert_equal(frames, [0, 3, 4, 5])

    def test_array_trajectory(self):
        frames = [ts.frame for ts in self.dcd[np.array([0, 3, 4, 5])]]
        assert_equal(frames, [0, 3, 4, 5])

    def test_list_reverse_trajectory(self):
        frames = [ts.frame for ts in self.dcd[[0, 4, 2, 3, 0, 1]]]
        assert_equal(frames, [0, 4, 2, 3, 0, 1])

    def test_list_repeated_trajectory(self):
        frames = [ts.frame for ts in self.dcd[[0, 0, 1, 1, 2, 1, 1]]]
        assert_equal(frames, [0, 0, 1, 1, 2, 1, 1])

    def test_reverse_dcd(self):
        frames = [ts.frame for ts in self.dcd[20:5:-1]]
        assert_equal(frames, list(range(20, 5, -1)),
                     "reversing dcd [20:5:-1]")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, 3341,
                     "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 98,
                     "wrong number of frames in dcd")

    def test_dt(self):
        assert_almost_equal(self.universe.trajectory.dt,
                            1.0,
                            4,
                            err_msg="wrong timestep dt")

    def test_totaltime(self):
        # test_totaltime(): need to reduce precision because dt is only precise
        # to ~4 decimals and accumulating the inaccuracy leads to even lower
        # precision in the totaltime (consequence of fixing Issue 64)
        assert_almost_equal(self.universe.trajectory.totaltime,
                            98.0,
                            3,
                            err_msg="wrong total length of AdK trajectory")

    def test_frame(self):
        self.dcd[15]  # index is 0-based and frames are 0-based
        assert_equal(self.universe.trajectory.frame, 15, "wrong frame number")

    def test_time(self):
        self.dcd[15]  # index is 0-based and frames are 0-based
        assert_almost_equal(self.universe.trajectory.time,
                            15.0,
                            5,
                            err_msg="wrong time of frame")

    def test_volume(self):
        assert_almost_equal(self.ts.volume, 0.0, 3,
                            err_msg="wrong volume for unitcell (no unitcell "
                            "in DCD so this should be 0)")

def test_DCDReader_set_dt(dt=100., frame=3):
    u = mda.Universe(PSF, DCD, dt=dt)
    assert_almost_equal(u.trajectory[frame].time, frame*dt,
                        err_msg="setting time step dt={0} failed: "
                        "actually used dt={1}".format(
            dt, u.trajectory._ts_kwargs['dt']))
    assert_almost_equal(u.trajectory.dt, dt,
                        err_msg="trajectory.dt does not match set dt")

class TestDCDWriter(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        ext = ".dcd"
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/dcd-writer' + ext
        self.Writer = mda.coordinates.DCD.DCDWriter

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe
        del self.Writer
        del self.tmpdir

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
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, 3,
                                      err_msg="coordinate mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                          orig_ts.frame, written_ts.frame))

    def test_dt(self):
        DT = 5.0
        t = self.universe.trajectory
        with self.Writer(self.outfile,
                         t.n_atoms,
                         dt=DT) as W:  # set time step to 5 ps
            for ts in self.universe.trajectory:
                W.write_next_timestep(ts)

        uw = mda.Universe(PSF, self.outfile)
        assert_almost_equal(uw.trajectory.totaltime,
                            uw.trajectory.n_frames * DT, 5)
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
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, 3,
                                      err_msg="coordinate mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                          orig_ts.frame, written_ts.frame))

    def test_single_frame(self):
        u = mda.Universe(PSF, CRD)
        W = mda.Writer(self.outfile, u.atoms.n_atoms)
        W.write(u.atoms)
        W.close()
        w = mda.Universe(PSF, self.outfile)
        assert_equal(w.trajectory.n_frames, 1,
                     "single frame trajectory has wrong number of frames")
        assert_almost_equal(w.atoms.positions,
                            u.atoms.positions,
                            3,
                            err_msg="coordinates do not match")

    def test_with_statement(self):
        u = mda.Universe(PSF, CRD)
        try:
            with mda.Writer(self.outfile, u.atoms.n_atoms) as W:
                W.write(u.atoms)
        except AttributeError:  # misses __exit__
            raise AssertionError("DCDWriter: does not support with statement")
        w = mda.Universe(PSF, self.outfile)
        assert_equal(w.trajectory.n_frames, 1,
                     "with_statement: single frame trajectory has wrong "
                     "number of frames")
        assert_almost_equal(w.atoms.positions,
                            u.atoms.positions,
                            3,
                            err_msg="with_statement: coordinates do not match")


class TestDCDWriter_Issue59(TestCase):
    def setUp(self):
        """Generate input xtc."""
        self.u = mda.Universe(PSF, DCD)
        self.tmpdir = tempdir.TempDir()
        self.xtc = self.tmpdir.name + '/dcd-writer-issue59-test.xtc'
        wXTC = mda.Writer(self.xtc, self.u.atoms.n_atoms)
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
        del self.tmpdir

    @attr('issue')
    def test_issue59(self):
        """Test writing of XTC to DCD (Issue 59)"""
        xtc = mda.Universe(PSF, self.xtc)
        self.dcd = self.tmpdir.name + '/dcd-writer-issue59-test.dcd'
        wDCD = mda.Writer(self.dcd, xtc.atoms.n_atoms)
        for ts in xtc.trajectory:
            wDCD.write(ts)
        wDCD.close()

        dcd = mda.Universe(PSF, self.dcd)

        xtc.trajectory.rewind()
        dcd.trajectory.rewind()

        assert_array_almost_equal(
            xtc.atoms.positions,
            dcd.atoms.positions,
            3,
            err_msg="XTC -> DCD: DCD coordinates are messed up (Issue 59)")

    def test_OtherWriter(self):
        dcd = self.u
        wXTC = dcd.trajectory.OtherWriter(self.xtc)
        for ts in dcd.trajectory:
            wXTC.write(ts)
        wXTC.close()

        xtc = mda.Universe(PSF, self.xtc)
        xtc.trajectory.rewind()
        dcd.trajectory.rewind()

        assert_array_almost_equal(
            dcd.atoms.positions,
            xtc.atoms.positions,
            2,
            err_msg="DCD -> XTC: coordinates are messed up (frame {0:d})".format(
            dcd.trajectory.frame))
        xtc.trajectory[3]
        dcd.trajectory[3]
        assert_array_almost_equal(
            dcd.atoms.positions,
            xtc.atoms.positions,
            2,
            err_msg="DCD -> XTC: coordinates are messed up (frame {0:d})".format(
            dcd.trajectory.frame))


class _TestDCDReader_TriclinicUnitcell(TestCase):
    def setUp(self):
        self.u = mda.Universe(self.topology, self.trajectory)
        self.tempdir = tempdir.TempDir()
        self.dcd = self.tempdir.name + '/dcd-reader-triclinic.dcd'

    def tearDown(self):
        try:
            os.unlink(self.dcd)
        except (AttributeError, OSError):
            pass
        del self.u
        del self.tempdir

    @attr('issue')
    def test_read_triclinic(self):
        """test reading of triclinic unitcell (Issue 187) for NAMD or new
        CHARMM format (at least since c36b2)"""
        for ts, box in zip(self.u.trajectory,
                           self.ref_dimensions[:, 1:]):
            assert_array_almost_equal(ts.dimensions, box, 4,
                                      err_msg="box dimensions A,B,C,alpha,"
                                      "beta,gamma not identical at frame "
                                      "{}".format(ts.frame))

    @attr('issue')
    def test_write_triclinic(self):
        """test writing of triclinic unitcell (Issue 187) for NAMD or new
        CHARMM format (at least since c36b2)"""
        with self.u.trajectory.OtherWriter(self.dcd) as w:
            for ts in self.u.trajectory:
                w.write(ts)
        w = mda.Universe(self.topology, self.dcd)
        for ts_orig, ts_copy in zip(self.u.trajectory,
                                    w.trajectory):
            assert_almost_equal(ts_orig.dimensions, ts_copy.dimensions, 4,
                                err_msg="DCD->DCD: unit cell dimensions wrong "
                                "at frame {0}".format(ts_orig.frame))
        del w


class TestDCDReader_CHARMM_Unitcell(_TestDCDReader_TriclinicUnitcell,
                                    RefCHARMMtriclinicDCD):
    pass


class TestDCDReader_NAMD_Unitcell(_TestDCDReader_TriclinicUnitcell,
                                  RefNAMDtriclinicDCD):
    pass


class TestNCDF2DCD(TestCase):
    @dec.skipif(module_not_found("netCDF4"),
                "Test skipped because netCDF is not available.")
    def setUp(self):
        self.u = mda.Universe(PRMncdf, NCDF)
        # create the DCD
        self.tmpdir = tempdir.TempDir()
        self.dcd = self.tmpdir.name + '/ncdf-2-dcd.dcd'
        DCD = mda.Writer(self.dcd, n_atoms=self.u.atoms.n_atoms)
        for ts in self.u.trajectory:
            DCD.write(ts)
        DCD.close()
        self.w = mda.Universe(PRMncdf, self.dcd)

    def tearDown(self):
        try:
            os.unlink(self.dcd)
        except (AttributeError, OSError):
            pass
        del self.u
        del self.w
        del self.tmpdir

    @attr('issue')
    def test_unitcell(self):
        """NCDFReader: Test that DCDWriter correctly writes the CHARMM
        unit cell"""
        for ts_orig, ts_copy in zip(self.u.trajectory,
                                    self.w.trajectory):
            assert_almost_equal(
                ts_orig.dimensions,
                ts_copy.dimensions,
                3,
                err_msg="NCDF->DCD: unit cell dimensions wrong at frame {0:d}".format(
                ts_orig.frame))

    def test_coordinates(self):
        for ts_orig, ts_copy in zip(self.u.trajectory,
                                    self.w.trajectory):
            assert_almost_equal(
                self.u.atoms.positions,
                self.w.atoms.positions,
                3,
                err_msg="NCDF->DCD: coordinates wrong at frame {0:d}".format(
                ts_orig.frame))


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
        phi151 = self.universe.select_atoms('resid 151').select_atoms(
            'name HN', 'name N', 'name CA', 'name CB')
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
        assert_almost_equal(C[3].__data__.mean(),
                            avg_angle,
                            err_msg="Correl: average Angle")

    def test_Dihedral(self):
        C = self.collection
        avg_phi151 = 0.0088003870749735619
        assert_almost_equal(C[4].__data__.mean(),
                            avg_phi151,
                            err_msg="Correl: average Dihedral")

    def test_scalarDistance(self):
        C = self.collection
        avg_dist = 9.7960210987736236
        assert_almost_equal(C[5].__data__.mean(),
                            avg_dist,
                            err_msg="Correl: average scalar Distance")

    def test_CenterOfMass(self):
        C = self.collection
        avg_com_ca = np.array([0.0043688, -0.27812258, 0.0284051])
        avg_com_all = np.array([-0.10086529, -0.16357276, 0.12724672])
        assert_array_almost_equal(C[6].__data__.mean(axis=1),
                                  avg_com_ca,
                                  err_msg="Correl: average CA CenterOfMass")
        assert_almost_equal(C[8].__data__.mean(axis=1),
                            avg_com_all,
                            err_msg="Correl: average all CenterOfMass")

    def test_CenterOfGeometry(self):
        C = self.collection
        avg_cog_all = np.array([-0.13554797, -0.20521885, 0.2118998])
        assert_almost_equal(C[9].__data__.mean(axis=1),
                            avg_cog_all,
                            err_msg="Correl: average all CenterOfGeometry")

    def test_CA_COMeqCOG(self):
        C = self.collection
        assert_array_almost_equal(
            C[6].__data__,
            C[7].__data__,
            err_msg="Correl: CA CentreOfMass == CenterOfGeometry")

    def test_clear(self):
        C = self.collection
        C.clear()
        assert_equal(len(C), 0, "Correl: clear()")


# notes:
def compute_correl_references():
    universe = mda.Universe(PSF, DCD)

    all = universe.atoms
    ca = universe.s4AKE.CA
    ca_termini = mda.core.AtomGroup.AtomGroup([ca[0], ca[-1]])
    phi151 = universe.select_atoms('resid 151').select_atoms(
        'name HN', 'name N', 'name CA', 'name CB')

    C = mda.collection
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


class TestDCDTimestep(BaseTimestepTest):
    Timestep = mda.coordinates.DCD.Timestep
    name = "DCD"
    has_box = True
    set_box = True
    unitcell = np.array([10., 90., 11., 90., 90., 12.])
    uni_args = (PSF, DCD)

    def test_ts_order_define(self):
        """Check that users can hack in a custom unitcell order"""
        old = self.Timestep._ts_order
        self.ts._ts_order = [0, 2, 5, 1, 3, 4]
        self.ts.dimensions = np.array([10, 11, 12, 80, 85, 90])
        assert_allclose(self.ts._unitcell, np.array([10, 80, 11, 85, 90, 12]))
        self.ts._ts_order = old
        self.ts.dimensions = np.zeros(6)
