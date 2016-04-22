import os
import numpy as np

import MDAnalysis as mda

from numpy.testing import (assert_equal, assert_almost_equal, assert_raises,
                           assert_)
from unittest import TestCase

from MDAnalysisTests.coordinates.reference import (RefLAMMPSData,
                                                   RefLAMMPSDataMini,
                                                   RefLAMMPSDataDCD)
from MDAnalysis.tests.datafiles import LAMMPScnt
from MDAnalysisTests import tempdir


def test_datareader_ValueError():
    from MDAnalysis.coordinates.LAMMPS import DATAReader
    assert_raises(ValueError, DATAReader, 'filename')


class _TestLammpsData_Coords(TestCase):
    """Tests using a .data file for loading single frame.

    All topology loading from MDAnalysisTests.data is done in test_topology
    """

    def setUp(self):
        self.u = mda.Universe(self.filename)

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
        assert_equal(type(ts), mda.coordinates.base.Timestep)

    def test_iter(self):
        # Check that iterating works, but only gives a single frame
        assert_equal(len(list(iter(self.u.trajectory))), 1)


class TestLammpsData_Coords(_TestLammpsData_Coords, RefLAMMPSData):
    pass


class TestLammpsDataMini_Coords(_TestLammpsData_Coords, RefLAMMPSDataMini):
    pass

# need more tests of the LAMMPS DCDReader

class TestLAMMPSDCDReader(TestCase, RefLAMMPSDataDCD):
    flavor = 'LAMMPS'

    def setUp(self):
        self.u = mda.Universe(self.topology, self.trajectory,
                              format=self.format)

    def tearDown(self):
        del self.u

    def test_Reader_is_LAMMPS(self):
        assert_(self.u.trajectory.flavor, self.flavor)

    def get_frame_from_end(self, offset):
        iframe = self.u.trajectory.n_frames - 1 - offset
        iframe = iframe if iframe > 0 else 0
        return iframe

    def test_n_atoms(self):
        assert_equal(self.u.atoms.n_atoms, self.n_atoms)

    def test_n_frames(self):
        assert_equal(self.u.trajectory.n_frames, self.n_frames)

    def test_dimensions(self):
        mean_dimensions = np.mean([ts.dimensions for ts in self.u.trajectory],
                                  axis=0)
        assert_almost_equal(mean_dimensions, self.mean_dimensions)

    def test_dt(self):
        assert_almost_equal(self.u.trajectory.dt, self.dt,
                            err_msg="Time between frames dt is wrong.")

    def test_Timestep_time(self):
        iframe = self.get_frame_from_end(1)
        assert_almost_equal(self.u.trajectory[iframe].time,
                            iframe * self.dt,
                            err_msg="Time for frame {0} (dt={1}) is wrong.".format(
                iframe, self.dt))

    def test_LAMMPSDCDReader_set_dt(self, dt=1500.):
        u = mda.Universe(self.topology, self.trajectory, format=self.format,
                         dt=dt)
        iframe = self.get_frame_from_end(1)
        assert_almost_equal(u.trajectory[iframe].time, iframe*dt,
                            err_msg="setting time step dt={0} failed: "
                            "actually used dt={1}".format(
                dt, u.trajectory._ts_kwargs['dt']))

    def test_wrong_time_unit(self):
        def wrong_load(unit="nm"):
            return mda.Universe(self.topology, self.trajectory, format=self.format,
                                timeunit=unit)
        assert_raises(TypeError, wrong_load)

    def test_wrong_unit(self):
        def wrong_load(unit="GARBAGE"):
            return mda.Universe(self.topology, self.trajectory, format=self.format,
                                timeunit=unit)
        assert_raises(ValueError, wrong_load)


class TestLAMMPSDCDWriter(TestCase, RefLAMMPSDataDCD):
    flavor = 'LAMMPS'

    def setUp(self):
        self.u = mda.Universe(self.topology, self.trajectory,
                              format=self.format)
        # dummy output file
        ext = os.path.splitext(self.trajectory)[1]
        self.tmpdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tmpdir.name,  'lammps-writer-test' + ext)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.u
        del self.tmpdir

    def test_Writer_is_LAMMPS(self):
        with mda.Writer(self.outfile, n_atoms=self.u.atoms.n_atoms,
                       format=self.format) as W:
            assert_(W.flavor, self.flavor)

    def test_Writer(self, n_frames=3):
        W = mda.Writer(self.outfile, n_atoms=self.u.atoms.n_atoms,
                       format=self.format)
        with self.u.trajectory.OtherWriter(self.outfile) as w:
            for ts in self.u.trajectory[:n_frames]:
                w.write(ts)
        short = mda.Universe(self.topology, self.outfile)
        assert_equal(short.trajectory.n_frames, n_frames,
                     err_msg="number of frames mismatch")
        assert_almost_equal(short.trajectory[n_frames - 1].positions,
                            self.u.trajectory[n_frames - 1].positions,
                            6, err_msg="coordinate mismatch between corresponding frames")

    def test_OtherWriter_is_LAMMPS(self):
        with self.u.trajectory.OtherWriter(self.outfile) as W:
            assert_(W.flavor, self.flavor)

    def test_OtherWriter(self):
        times = []
        with self.u.trajectory.OtherWriter(self.outfile) as w:
            for ts in self.u.trajectory[::-1]:
                times.append(ts.time)
                w.write(ts)
        # note: the reversed trajectory records times in increasing
        #       steps, and NOT reversed, i.e. the time markers are not
        #       attached to their frames. This could be considered a bug
        #       but DCD has no way to store timestamps. Right now, we'll simply
        #       test that this is the case and pass.
        reversed = mda.Universe(self.topology, self.outfile)
        assert_equal(reversed.trajectory.n_frames, self.u.trajectory.n_frames,
                     err_msg="number of frames mismatch")
        rev_times = [ts.time for ts in reversed.trajectory]
        assert_almost_equal(rev_times, times[::-1], 6,
                            err_msg="time steps of written DCD mismatch")
        assert_almost_equal(reversed.trajectory[-1].positions,
                            self.u.trajectory[0].positions,
                            6, err_msg="coordinate mismatch between corresponding frames")

class TestLAMMPSDCDWriterClass(TestCase):
    flavor = 'LAMMPS'

    def setUp(self):
        # dummy output file
        ext = ".dcd"
        self.tmpdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tmpdir.name,  'lammps-writer-test' + ext)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.tmpdir

    def test_Writer_is_LAMMPS(self):
        with mda.coordinates.LAMMPS.DCDWriter(self.outfile, n_atoms=10) as W:
            assert_(W.flavor, self.flavor)

    def test_open(self):
        def open_dcd():
            try:
                with mda.coordinates.LAMMPS.DCDWriter(self.outfile, n_atoms=10):
                    pass
            except Exception:
                return False
            else:
                return True
        assert_(open_dcd(), True)

    def test_wrong_time_unit(self):
        def wrong_load(unit="nm"):
                with mda.coordinates.LAMMPS.DCDWriter(self.outfile, n_atoms=10,
                                                      timeunit=unit):
                    pass
        assert_raises(TypeError, wrong_load)

    def test_wrong_unit(self):
        def wrong_load(unit="GARBAGE"):
                with mda.coordinates.LAMMPS.DCDWriter(self.outfile, n_atoms=10,
                                                      timeunit=unit):
                    pass
        assert_raises(ValueError, wrong_load)


class TestLammpsDataTriclinic(object):
    def setUp(self):
        self.u = mda.Universe(LAMMPScnt)

    def tearDown(self):
        del self.u

    def test_triclinicness(self):
        assert_(self.u.dimensions[3] == 90.)
        assert_(self.u.dimensions[4] == 90.)
        assert_(self.u.dimensions[5] == 120.)
