import numpy as np

import MDAnalysis as mda

from numpy.testing import (assert_equal, assert_almost_equal, assert_raises)
from unittest import TestCase

from MDAnalysisTests.coordinates.reference import (RefLAMMPSData,
                                                   RefLAMMPSDataMini,
                                                   RefLAMMPSDataDCD)


def test_datareader_VE():
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

class TestLAMPPSDCDReader(TestCase, RefLAMMPSDataDCD):
    def setUp(self):
        self.u = mda.Universe(self.topology, self.trajectory,
                              format=self.format)

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
