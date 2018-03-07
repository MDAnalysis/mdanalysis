# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
import MDAnalysis as mda
import numpy as np
import sys
from six.moves import zip

from scipy.io import netcdf

import pytest
from numpy.testing import (
    assert_equal,
    assert_almost_equal
)

from MDAnalysisTests.datafiles import (PFncdf_Top, PFncdf_Trj,
                                       GRO, TRR, XYZ_mini)
from MDAnalysisTests.coordinates.test_trj import _TRJReaderTest
from MDAnalysisTests.coordinates.reference import (RefVGV, RefTZ2)
from MDAnalysisTests import make_Universe


class _NCDFReaderTest(_TRJReaderTest):
    prec = 3

    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename)

    def test_slice_iteration(self, universe):
        frames = [ts.frame for ts in universe.trajectory[4:-2:4]]
        assert_equal(frames,
                     np.arange(universe.trajectory.n_frames)[4:-2:4],
                     err_msg="slicing did not produce the expected frames")

    def test_metadata(self, universe):
        data = universe.trajectory.trjfile
        assert_equal(data.Conventions.decode('utf-8'), 'AMBER')
        assert_equal(data.ConventionVersion.decode('utf-8'), '1.0')

    def test_dt(self, universe):
        ref = 0.0
        assert_almost_equal(ref, universe.trajectory.dt, self.prec)
        assert_almost_equal(ref, universe.trajectory.ts.dt, self.prec)

    def test_get_writer(self, universe):
        with universe.trajectory.Writer('out.ncdf') as w:
            assert w.n_atoms == len(universe.atoms)
            assert w.remarks.startswith('AMBER NetCDF format')

    def test_get_writer_custom_n_atoms(self, universe):
        with universe.trajectory.Writer('out.ncdf', n_atoms=42,
                                        remarks='Hi!') as w:
            assert w.n_atoms == 42
            assert w.remarks == 'Hi!'

    def test_wrong_natoms(self):
        with pytest.raises(ValueError):
            mda.coordinates.TRJ.NCDFReader(self.filename, n_atoms=2)

    def test_read_on_closed(self, universe):
        universe.trajectory.close()

        with pytest.raises(IOError):
            universe.trajectory.__getitem__(2)

    def test_mmap_kwarg(self, universe):
        # default is None
        assert universe.trajectory._mmap == None


# Ugly way to create the tests for mmap

class _NCDFReaderTest_mmap_None(_NCDFReaderTest):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename, mmap=None)

class _NCDFReaderTest_mmap_True(_NCDFReaderTest):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename, mmap=True)

    def test_mmap_kwarg(self, universe):
        # default is None
        assert universe.trajectory._mmap == True

class _NCDFReaderTest_mmap_False(_NCDFReaderTest):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename, mmap=False)

    def test_mmap_kwarg(self, universe):
        assert universe.trajectory._mmap == False


class TestNCDFReader(_NCDFReaderTest, RefVGV):
    pass

class TestNCDFReader_mmap_None(_NCDFReaderTest_mmap_None, RefVGV):
    pass

class TestNCDFReader_mmap_True(_NCDFReaderTest_mmap_True, RefVGV):
    pass

class TestNCDFReader_mmap_False(_NCDFReaderTest_mmap_False, RefVGV):
    pass



class TestNCDFReaderTZ2(_NCDFReaderTest, RefTZ2):
    pass


class TestNCDFReader2(object):
    """NCDF Trajectory with positions and forces.

    Contributed by Albert Solernou
    """
    prec = 3

    @pytest.fixture(scope='class')
    def u(self):
        return mda.Universe(PFncdf_Top, PFncdf_Trj)

    def test_positions_1(self, u):
        """Check positions on first frame"""
        u.trajectory[0]
        ref_1 = np.array([[-0.11980818, 18.70524979, 11.6477766
                           ], [-0.44717646, 18.61727142, 12.59919548],
                          [-0.60952115, 19.47885513, 11.22137547]],
                         dtype=np.float32)
        assert_almost_equal(ref_1, u.atoms.positions[:3], self.prec)

    def test_positions_2(self, u):
        """Check positions on second frame"""
        u.trajectory[1]
        ref_2 = np.array([[-0.13042036, 18.6671524, 11.69647026
                           ], [-0.46643803, 18.60186768, 12.646698],
                          [-0.46567637, 19.49173927, 11.21922874]],
                         dtype=np.float32)
        assert_almost_equal(ref_2, u.atoms.positions[:3], self.prec)

    def test_forces_1(self, u):
        """Check forces on first frame"""
        u.trajectory[0]
        ref_1 = np.array([[49.23017883, -97.05565643, -86.09863281
                           ], [2.97547197, 29.84169388, 11.12069607],
                          [-15.93093777, 14.43616867, 30.25889015]],
                         dtype=np.float32)
        assert_almost_equal(ref_1, u.atoms.forces[:3], self.prec)

    def test_forces_2(self, u):
        """Check forces on second frame"""
        u.trajectory[1]
        ref_2 = np.array([[116.39096832, -145.44448853, -151.3155365
                           ], [-18.90058327, 27.20145798, 1.95245135],
                          [-31.08556366, 14.95863628, 41.10367966]],
                         dtype=np.float32)
        assert_almost_equal(ref_2, u.atoms.forces[:3], self.prec)

    def test_time_1(self, u):
        """Check time on first frame"""
        ref = 35.02
        assert_almost_equal(ref, u.trajectory[0].time, self.prec)

    def test_time_2(self, u):
        """Check time on second frame"""
        ref = 35.04
        assert_almost_equal(ref, u.trajectory[1].time, self.prec)

    def test_dt(self, u):
        ref = 0.02
        assert_almost_equal(ref, u.trajectory.dt, self.prec)
        assert_almost_equal(ref, u.trajectory.ts.dt, self.prec)


class _NCDFWriterTest(object):
    prec = 5

    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename)

    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + 'ncdf-writer-1.ncdf'

    @pytest.fixture()
    def outtop(self, tmpdir):
        return str(tmpdir) + 'ncdf-writer-top.pdb'

    def _test_write_trajectory(self, universe, outfile):
        # explicit import so that we can artifically remove netCDF4
        # before calling
        from MDAnalysis.coordinates import TRJ

        t = universe.trajectory
        with TRJ.NCDFWriter(outfile, t.n_atoms, dt=t.dt) as W:
            self._copy_traj(W, universe)
        self._check_new_traj(universe, outfile)
        # for issue #518 -- preserve float32 data in ncdf output
        # NOTE: This originally failed with the dtype('>f4') instead
        #       of dtype('<f4') == dtype('f') == np.float32, i.e. then
        #       endianness is different. The current hack-ish solution
        #       ignores endianness by comparing the name of the types,
        #       which should be "float32".
        #       See http://docs.scipy.org/doc/numpy-1.10.0/reference/arrays.dtypes.html
        #       and https://github.com/MDAnalysis/mdanalysis/pull/503
        dataset = netcdf.netcdf_file(outfile, 'r')
        coords = dataset.variables['coordinates']
        time = dataset.variables['time']
        assert_equal(coords[:].dtype.name, np.dtype(np.float32).name,
                     err_msg='ncdf coord output not float32 '
                             'but {}'.format(coords[:].dtype))
        assert_equal(time[:].dtype.name, np.dtype(np.float32).name,
                     err_msg='ncdf time output not float32 '
                             'but {}'.format(time[:].dtype))

    def test_write_trajectory_netCDF4(self, universe, outfile):
        pytest.importorskip("netCDF4")
        return self._test_write_trajectory(universe, outfile)

    def test_write_trajectory_netcdf(self, universe, outfile):
        import MDAnalysis.coordinates.TRJ
        loaded_netCDF4 = sys.modules['MDAnalysis.coordinates.TRJ'].netCDF4
        try:
            # cannot use @block_import('netCDF4') because TRJ was already imported
            # during setup() and already sits in the global module list so we just
            # set it to None because that is what TRJ does if it cannot find netCDF4
            sys.modules['MDAnalysis.coordinates.TRJ'].netCDF4 = None
            assert MDAnalysis.coordinates.TRJ.netCDF4 is None  # should happen if netCDF4 not found
            return self._test_write_trajectory(universe, outfile)
        finally:
            sys.modules['MDAnalysis.coordinates.TRJ'].netCDF4 = loaded_netCDF4

    def test_OtherWriter(self, universe, outfile):
        t = universe.trajectory
        with t.OtherWriter(outfile) as W:
            self._copy_traj(W, universe)
        self._check_new_traj(universe, outfile)

    def _copy_traj(self, writer, universe):
        for ts in universe.trajectory:
            writer.write_next_timestep(ts)

    def _check_new_traj(self, universe, outfile):
        uw = mda.Universe(self.topology, outfile)

        # check that the trajectories are identical for each time step
        for orig_ts, written_ts in zip(universe.trajectory,
                                       uw.trajectory):
            assert_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between "
                                              "original and written trajectory at "
                                              "frame %d (orig) vs %d (written)" % (
                                                  orig_ts.frame,
                                                  written_ts.frame))
            # not a good test because in the example trajectory all times are 0
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the "
                                        "same.".format(orig_ts.frame))
            assert_almost_equal(written_ts.dimensions,
                                      orig_ts.dimensions,
                                      self.prec,
                                      err_msg="unitcells are not identical")
        # check that the NCDF data structures are the same
        nc_orig = universe.trajectory.trjfile
        nc_copy = uw.trajectory.trjfile

        # note that here 'dimensions' is a specific netcdf data structure and
        # not the unit cell dimensions in MDAnalysis
        for k, dim in nc_orig.dimensions.items():
            try:
                dim_new = nc_copy.dimensions[k]
            except KeyError:
                raise AssertionError("NCDFWriter did not write "
                                     "dimension '{0}'".format(k))
            else:
                assert_equal(dim, dim_new,
                             err_msg="Dimension '{0}' size mismatch".format(k))

        for k, v in nc_orig.variables.items():
            try:
                v_new = nc_copy.variables[k]
            except KeyError:
                raise AssertionError("NCDFWriter did not write "
                                     "variable '{0}'".format(k))
            else:
                try:
                    assert_almost_equal(v[:], v_new[:], self.prec,
                                              err_msg="Variable '{0}' not "
                                                      "written correctly".format(
                                                  k))
                except TypeError:
                    assert_equal(v[:], v_new[:],
                                       err_msg="Variable {0} not written "
                                               "correctly".format(k))

    def test_TRR2NCDF(self, outfile):
        trr = mda.Universe(GRO, TRR)
        with mda.Writer(outfile, trr.trajectory.n_atoms,
                        velocities=True, format="ncdf") as W:
            for ts in trr.trajectory:
                W.write_next_timestep(ts)

        uw = mda.Universe(GRO, outfile)

        for orig_ts, written_ts in zip(trr.trajectory,
                                       uw.trajectory):
            assert_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between "
                                              "original and written trajectory at "
                                              "frame {0} (orig) vs {1} (written)".format(
                                          orig_ts.frame, written_ts.frame))
            assert_almost_equal(written_ts._velocities,
                                      orig_ts._velocities, self.prec,
                                      err_msg="velocity mismatch between "
                                              "original and written trajectory at "
                                              "frame {0} (orig) vs {1} (written)".format(
                                          orig_ts.frame, written_ts.frame))
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the "
                                        "same.".format(orig_ts.frame))
            assert_almost_equal(written_ts.dimensions,
                                      orig_ts.dimensions,
                                      self.prec,
                                      err_msg="unitcells are not identical")
        del trr

    def test_write_AtomGroup(self, universe, outfile, outtop):
        """test to write NCDF from AtomGroup (Issue 116)"""
        p = universe.select_atoms("not resname WAT")
        p.write(outtop)
        with mda.Writer(outfile, n_atoms=p.n_atoms, format="ncdf") as W:
            for ts in universe.trajectory:
                W.write(p)

        uw = mda.Universe(outtop, outfile)
        pw = uw.atoms

        for orig_ts, written_ts in zip(universe.trajectory,
                                       uw.trajectory):
            assert_almost_equal(p.positions, pw.positions, self.prec,
                                      err_msg="coordinate mismatch between "
                                              "original and written trajectory at "
                                              "frame %d (orig) vs %d (written)" % (
                                                  orig_ts.frame,
                                                  written_ts.frame))
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the "
                                        "same.".format(orig_ts.frame))
            assert_almost_equal(written_ts.dimensions,
                                      orig_ts.dimensions,
                                      self.prec,
                                      err_msg="unitcells are not identical")


class TestNCDFWriter(_NCDFWriterTest, RefVGV):
    pass


class TestNCDFWriterTZ2(_NCDFWriterTest, RefTZ2):
    pass


class TestNCDFWriterVelsForces(object):
    """Test writing NCDF trajectories with a mixture of options"""
    prec = 3
    top = XYZ_mini
    n_atoms = 3

    @pytest.fixture()
    def ts1(self):
        ts = mda.coordinates.TRJ.Timestep(self.n_atoms, velocities=True,
                                          forces=True)
        ts._pos[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms, 3)
        ts._velocities[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms,
                                                                3) + 100
        ts._forces[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms,
                                                            3) + 200
        return ts

    @pytest.fixture()
    def ts2(self):
        ts = mda.coordinates.TRJ.Timestep(self.n_atoms, velocities=True,
                                          forces=True)
        ts._pos[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms, 3) + 300
        ts._velocities[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms,
                                                                3) + 400
        ts._forces[:] = np.arange(self.n_atoms * 3).reshape(self.n_atoms,
                                                            3) + 500
        return ts

    @pytest.mark.parametrize('pos, vel, force', (
            (True, False, False),
            (True, True, False),
            (True, False, True),
            (True, True, True),
    ))
    def test_write_ts(self, pos, vel, force, tmpdir, ts1, ts2):
        """Write the two reference timesteps, then open them up and check values

        pos vel and force are bools which define whether these properties
        should be in TS

        """
        outfile = str(tmpdir) + 'ncdf-write-vels-force.ncdf'
        with mda.Writer(outfile,
                        n_atoms=self.n_atoms,
                        velocities=vel,
                        forces=force) as w:
            w.write(ts1)
            w.write(ts2)

        u = mda.Universe(self.top, outfile)
        for ts, ref_ts in zip(u.trajectory, [ts1, ts2]):
            if pos:
                assert_almost_equal(ts._pos, ref_ts._pos, self.prec)
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(ts, 'positions')
            if vel:
                assert_almost_equal(ts._velocities, ref_ts._velocities,
                                    self.prec)
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(ts, 'velocities')
            if force:
                assert_almost_equal(ts._forces, ref_ts._forces, self.prec)
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(ts, 'forces')

        u.trajectory.close()


class TestNCDFWriterErrors(object):
    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + 'out.ncdf'

    def test_zero_atoms_VE(self, outfile):
        from MDAnalysis.coordinates.TRJ import NCDFWriter

        with pytest.raises(ValueError):
            NCDFWriter(outfile, 0)

    def test_wrong_n_atoms(self, outfile):
        from MDAnalysis.coordinates.TRJ import NCDFWriter

        with NCDFWriter(outfile, 100) as w:
            u = make_Universe(trajectory=True)
            with pytest.raises(IOError):
                w.write(u.trajectory.ts)

    def test_no_ts(self, outfile):
        # no ts supplied at any point
        from MDAnalysis.coordinates.TRJ import NCDFWriter

        with NCDFWriter(outfile, 100) as w:
            with pytest.raises(IOError):
                w.write_next_timestep()
