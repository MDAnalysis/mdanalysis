import MDAnalysis as mda
import numpy as np
import os
from six.moves import zip

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_array_almost_equal,
                           assert_almost_equal, assert_raises)
import tempdir
from unittest import TestCase

from MDAnalysisTests.datafiles import (PRMncdf, NCDF, PFncdf_Top, PFncdf_Trj,
                                       GRO, TRR, XYZ_mini)
from MDAnalysisTests.coordinates.test_trj import _TRJReaderTest
from MDAnalysisTests.coordinates.reference import (RefVGV,)


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
        ref_1 = np.array([[-0.11980818, 18.70524979, 11.6477766
                           ], [-0.44717646, 18.61727142, 12.59919548],
                          [-0.60952115, 19.47885513, 11.22137547]],
                         dtype=np.float32)
        assert_array_almost_equal(ref_1, self.u.atoms.positions[:3], self.prec)

    def test_positions_2(self):
        """Check positions on second frame"""
        self.u.trajectory[1]
        ref_2 = np.array([[-0.13042036, 18.6671524, 11.69647026
                           ], [-0.46643803, 18.60186768, 12.646698],
                          [-0.46567637, 19.49173927, 11.21922874]],
                         dtype=np.float32)
        assert_array_almost_equal(ref_2, self.u.atoms.positions[:3], self.prec)

    def test_forces_1(self):
        """Check forces on first frame"""
        self.u.trajectory[0]
        ref_1 = np.array([[49.23017883, -97.05565643, -86.09863281
                           ], [2.97547197, 29.84169388, 11.12069607],
                          [-15.93093777, 14.43616867, 30.25889015]],
                         dtype=np.float32)
        assert_array_almost_equal(ref_1, self.u.atoms.forces[:3], self.prec)

    def test_forces_2(self):
        """Check forces on second frame"""
        self.u.trajectory[1]
        ref_2 = np.array([[116.39096832, -145.44448853, -151.3155365
                           ], [-18.90058327, 27.20145798, 1.95245135],
                          [-31.08556366, 14.95863628, 41.10367966]],
                         dtype=np.float32)
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
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/ncdf-writer-1' + ext
        self.outtop = self.tmpdir.name + '/ncdf-writer-top.pdb'
        self.Writer = mda.coordinates.TRJ.NCDFWriter

    def tearDown(self):
        for f in self.outfile, self.outtop:
            try:
                os.unlink(f)
            except OSError:
                pass
        del self.universe
        del self.Writer
        del self.tmpdir

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
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                          orig_ts.frame, written_ts.frame))
            # not a good test because in the example trajectory all times are 0
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the "
                                "same.".format(orig_ts.frame))
            assert_array_almost_equal(written_ts.dimensions,
                                      orig_ts.dimensions,
                                      self.prec,
                                      err_msg="unitcells are not identical")

    @attr('slow')
    def test_TRR2NCDF(self):
        trr = mda.Universe(GRO, TRR)
        W = self.Writer(self.outfile, trr.trajectory.n_atoms, velocities=True)
        for ts in trr.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(GRO, self.outfile)

        for orig_ts, written_ts in zip(trr.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                          orig_ts.frame, written_ts.frame))
            assert_array_almost_equal(written_ts._velocities,
                                      orig_ts._velocities, self.prec,
                                      err_msg="velocity mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                          orig_ts.frame, written_ts.frame))
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the "
                                "same.".format(orig_ts.frame))
            assert_array_almost_equal(written_ts.dimensions,
                                      orig_ts.dimensions,
                                      self.prec,
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

        uw = mda.Universe(self.outtop, self.outfile)
        pw = uw.atoms

        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(p.positions, pw.positions, self.prec,
                                      err_msg="coordinate mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                          orig_ts.frame, written_ts.frame))
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the "
                                "same.".format(orig_ts.frame))
            assert_array_almost_equal(written_ts.dimensions,
                                      orig_ts.dimensions,
                                      self.prec,
                                      err_msg="unitcells are not identical")


class TestNCDFWriterVelsForces(TestCase):
    """Test writing NCDF trajectories with a mixture of options"""

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/ncdf-write-vels-force.ncdf'
        self.prec = 3
        self.top = XYZ_mini
        self.n_atoms = 3

        self.ts1 = mda.coordinates.TRJ.Timestep(self.n_atoms,
                                                velocities=True,
                                                forces=True)
        self.ts1._pos[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms, 3)
        self.ts1._velocities[:] = np.arange(self.n_atoms * 3).reshape(
            self.n_atoms, 3) + 100
        self.ts1._forces[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms,
                                                                  3) + 200

        self.ts2 = mda.coordinates.TRJ.Timestep(self.n_atoms,
                                                velocities=True,
                                                forces=True)
        self.ts2._pos[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms,
                                                               3) + 300
        self.ts2._velocities[:] = np.arange(self.n_atoms * 3).reshape(
            self.n_atoms, 3) + 400
        self.ts2._forces[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms,
                                                                  3) + 500

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass

        del self.n_atoms
        del self.ts1
        del self.ts2
        del self.tmpdir

    def _write_ts(self, pos, vel, force):
        """Write the two reference timesteps, then open them up and check values

        pos vel and force are bools which define whether these properties
        should be in TS

        """
        with mda.Writer(self.outfile,
                        n_atoms=self.n_atoms,
                        velocities=vel,
                        forces=force) as w:
            w.write(self.ts1)
            w.write(self.ts2)

        u = mda.Universe(self.top, self.outfile)
        for ts, ref_ts in zip(u.trajectory, [self.ts1, self.ts2]):
            if pos:
                assert_almost_equal(ts._pos, ref_ts._pos, self.prec)
            else:
                assert_raises(mda.NoDataError, getattr, ts, 'positions')
            if vel:
                assert_almost_equal(ts._velocities, ref_ts._velocities,
                                    self.prec)
            else:
                assert_raises(mda.NoDataError, getattr, ts, 'velocities')
            if force:
                assert_almost_equal(ts._forces, ref_ts._forces, self.prec)
            else:
                assert_raises(mda.NoDataError, getattr, ts, 'forces')

        u.trajectory.close()

    def test_pos(self):
        self._write_ts(True, False, False)

    def test_pos_vel(self):
        self._write_ts(True, True, False)

    def test_pos_force(self):
        self._write_ts(True, False, True)

    def test_pos_vel_force(self):
        self._write_ts(True, True, True)
