import MDAnalysis as mda
import numpy as np
import os
from six.moves import zip

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_array_almost_equal,
                           assert_array_equal,
                           assert_almost_equal, assert_raises, dec)
from unittest import TestCase
from MDAnalysisTests import module_not_found

from MDAnalysisTests.datafiles import (PRMncdf, NCDF, PFncdf_Top, PFncdf_Trj,
                                       GRO, TRR, XYZ_mini)
from MDAnalysisTests.coordinates.test_trj import _TRJReaderTest
from MDAnalysisTests.coordinates.reference import (RefVGV, RefTZ2)
from MDAnalysisTests import tempdir


class _NCDFReaderTest(_TRJReaderTest):
    @dec.skipif(module_not_found("netCDF4"), "Test skipped because netCDF is not available.")
    def setUp(self):
        self.universe = mda.Universe(self.topology, self.filename)
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

    def test_dt(self):
        ref = 0.0
        assert_almost_equal(ref, self.universe.trajectory.dt, self.prec)
        assert_almost_equal(ref, self.universe.trajectory.ts.dt, self.prec)


class TestNCDFReader(_NCDFReaderTest, RefVGV):
    pass

class TestNCDFReaderTZ2(_NCDFReaderTest, RefTZ2):
    pass


class TestNCDFReader2(TestCase):
    """NCDF Trajectory with positions and forces.

    Contributed by Albert Solernou
    """

    @dec.skipif(module_not_found("netCDF4"), "Test skipped because netCDF is not available.")
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

    def test_dt(self):
        ref = 0.02
        assert_almost_equal(ref, self.u.trajectory.dt, self.prec)
        assert_almost_equal(ref, self.u.trajectory.ts.dt, self.prec)


class _NCDFWriterTest(TestCase):
    @dec.skipif(module_not_found("netCDF4"), "Test skipped because netCDF is not available.")
    def setUp(self):
        self.universe = mda.Universe(self.topology, self.filename)
        self.prec = 5
        ext = ".ncdf"
        self.tmpdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tmpdir.name, 'ncdf-writer-1' + ext)
        self.outtop = os.path.join(self.tmpdir.name, 'ncdf-writer-top.pdb')
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
        with self.Writer(self.outfile, t.n_atoms, dt=t.dt) as W:
            self._copy_traj(W)
        self._check_new_traj()
        import netCDF4
        #for issue #518 -- preserve float32 data in ncdf output
        dataset = netCDF4.Dataset(self.outfile, 'r', format='NETCDF3')
        coords = dataset.variables['coordinates']
        time = dataset.variables['time']
        assert_equal(coords.dtype, np.float32,
                     err_msg='ncdf coord output not float32')
        assert_equal(time.dtype, np.float32,
                     err_msg='ncdf time output not float32')

    def test_OtherWriter(self):
        t = self.universe.trajectory
        with t.OtherWriter(self.outfile) as W:
            self._copy_traj(W)
        self._check_new_traj()

    def _copy_traj(self, writer):
        for ts in self.universe.trajectory:
            writer.write_next_timestep(ts)

    def _check_new_traj(self):
        uw = mda.Universe(self.topology, self.outfile)

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
        # check that the NCDF data structures are the same
        nc_orig = self.universe.trajectory.trjfile
        nc_copy = uw.trajectory.trjfile

        for k, dim in nc_orig.dimensions.items():
            try:
                dim_new = nc_copy.dimensions[k]
            except KeyError:
                raise AssertionError("NCDFWriter did not write "
                                     "dimension '{}'".format(k))
            else:
                assert_equal(len(dim), len(dim_new),
                             err_msg="Dimension '{0}' size mismatch".format(k))


        for k, v in nc_orig.variables.items():
            try:
                v_new = nc_copy.variables[k]
            except KeyError:
                raise AssertionError("NCDFWriter did not write "
                                     "variable '{}'".format(k))
            else:
                try:
                    assert_array_almost_equal(v[:], v_new[:], self.prec,
                                              err_msg="Variable '{}' not "
                                              "written correctly".format(k))
                except TypeError:
                    assert_array_equal(v[:], v_new[:],
                                              err_msg="Variable {} not written "
                                    "correctly".format(k))

    @attr('slow')
    def test_TRR2NCDF(self):
        trr = mda.Universe(GRO, TRR)
        with self.Writer(self.outfile, trr.trajectory.n_atoms,
                         velocities=True) as W:
            for ts in trr.trajectory:
                W.write_next_timestep(ts)

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
        with self.Writer(self.outfile, n_atoms=p.n_atoms) as W:
            for ts in self.universe.trajectory:
                W.write(p)

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

class TestNCDFWriter(_NCDFWriterTest, RefVGV):
    pass

class TestNCDFWriterTZ2(_NCDFWriterTest, RefTZ2):
    pass

class TestNCDFWriterVelsForces(TestCase):
    """Test writing NCDF trajectories with a mixture of options"""

    @dec.skipif(module_not_found("netCDF4"), "Test skipped because netCDF is not available.")
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
        except OSError:
            pass

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
