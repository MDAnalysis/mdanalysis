import pytest
from numpy.testing import assert_almost_equal, assert_array_equal, assert_equal
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.H5MD import HAS_H5PY
if HAS_H5PY:
    import h5py
from MDAnalysis.exceptions import NoDataError
from MDAnalysisTests import make_Universe
from MDAnalysisTests.datafiles import (H5MD_xvf, TPR_xvf,
                                       COORDINATES_TOPOLOGY,
                                       COORDINATES_H5MD)
from MDAnalysisTests.coordinates.base import (MultiframeReaderTest,
                                              BaseReference, BaseWriterTest,
                                              assert_timestep_almost_equal)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
class H5MDReference(BaseReference):
    """Reference synthetic trajectory that was
    copied from test_xdr.TRRReference"""

    def __init__(self):
        super(H5MDReference, self).__init__()
        self.trajectory = COORDINATES_H5MD
        self.topology = COORDINATES_TOPOLOGY
        self.reader = mda.coordinates.H5MD.H5MDReader
        self.writer = mda.coordinates.H5MD.H5MDWriter
        self.ext = 'h5md'
        self.prec = 3
        self.changing_dimensions = True

        self.first_frame.velocities = self.first_frame.positions / 10
        self.first_frame.forces = self.first_frame.positions / 100

        self.second_frame.velocities = self.second_frame.positions / 10
        self.second_frame.forces = self.second_frame.positions / 100

        self.last_frame.velocities = self.last_frame.positions / 10
        self.last_frame.forces = self.last_frame.positions / 100

        self.jump_to_frame.velocities = self.jump_to_frame.positions / 10
        self.jump_to_frame.forces = self.jump_to_frame.positions / 100

    def iter_ts(self, i):
        ts = self.first_frame.copy()
        ts.positions = 2**i * self.first_frame.positions
        ts.velocities = ts.positions / 10
        ts.forces = ts.positions / 100
        ts.time = i
        ts.frame = i
        return ts


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
class TestH5MDReader(MultiframeReaderTest):
    """Tests H5MDReader with MultiframeReaderTest."""
    @staticmethod
    @pytest.fixture()
    def ref():
        return H5MDReference()

    #@pytest.mark.xfail(reason='H5MD writer not implemented yet')
    def test_get_writer_1(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = 'test-writer.' + ref.ext
            with reader.Writer(outfile) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, reader.n_atoms)

    #@pytest.mark.xfail(reason='H5MD writer not implemented yet')
    def test_get_writer_2(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = 'test-writer.' + ref.ext
            with reader.Writer(outfile, n_atoms=100) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, 100)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
class TestH5MDWriter(BaseWriterTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return H5MDReference()

    def test_write_trajectory_atomgroup(self, ref,reader, universe, tmpdir):
        outfile = 'write-atoms-test.' + ref.ext
        with tmpdir.as_cwd():
            with ref.writer(outfile, universe.atoms.n_atoms,
                            velocities=True, forces=True) as w:
                for ts in universe.trajectory:
                    w.write(universe.atoms)
            self._check_copy(outfile, ref, reader)

    def test_write_trajectory_universe(self, ref, reader, universe, tmpdir):
        outfile = 'write-uni-test.' + ref.ext
        with tmpdir.as_cwd():
            with ref.writer(outfile, universe.atoms.n_atoms,
                            velocities=True, forces=True) as w:
                for ts in universe.trajectory:
                    w.write(universe)
            self._check_copy(outfile, ref, reader)

@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
class TestH5MDReader_2(object):

    prec = 3
    ext = 'h5md'

    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(TPR_xvf, H5MD_xvf)

    @pytest.fixture()
    def h5md_file(self):
        return h5py.File(H5MD_xvf, 'r')

    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir.join('h5md-reader-test.' + self.ext))

    def test_n_frames(self, universe):
        assert len(universe.trajectory) == 3

    def test_positions(self, universe):
        universe.trajectory[0]
        assert_almost_equal(universe.atoms.positions[0],
                            [32.309906, 13.77798, 14.372463],
                            decimal=self.prec)
        assert_almost_equal(universe.atoms.positions[42],
                            [28.116928, 19.405945, 19.647358],
                            decimal=self.prec)
        assert_almost_equal(universe.atoms.positions[10000],
                            [44.117805, 50.442093, 23.299038],
                            decimal=self.prec)

        universe.trajectory[1]
        assert_almost_equal(universe.atoms.positions[0],
                            [30.891968, 13.678971, 13.6000595],
                            decimal=self.prec)
        assert_almost_equal(universe.atoms.positions[42],
                            [27.163246, 19.846561, 19.3582],
                            decimal=self.prec)
        assert_almost_equal(universe.atoms.positions[10000],
                            [45.869278, 5.0342298, 25.460655],
                            decimal=self.prec)
        universe.trajectory[2]
        assert_almost_equal(universe.atoms.positions[0],
                            [31.276512, 13.89617, 15.015897],
                            decimal=self.prec)
        assert_almost_equal(universe.atoms.positions[42],
                            [28.567991, 20.56532, 19.40814],
                            decimal=self.prec)
        assert_almost_equal(universe.atoms.positions[10000],
                            [39.713223,  6.127234, 18.284992],
                            decimal=self.prec)

    def test_h5md_velocities(self, universe):
        universe.trajectory[0]
        assert_almost_equal(universe.atoms.velocities[0],
                            [-2.697732, 0.613568, 0.14334752],
                            decimal=self.prec)
        universe.trajectory[1]
        assert_almost_equal(universe.atoms.velocities[42],
                            [-6.8698354, 7.834235, -8.114698],
                            decimal=self.prec)
        universe.trajectory[2]
        assert_almost_equal(universe.atoms.velocities[10000],
                            [9.799492, 5.631466, 6.852126],
                            decimal=self.prec)

    def test_h5md_forces(self, universe):
        universe.trajectory[0]
        assert_almost_equal(universe.atoms.forces[0],
                            [20.071287, -155.2285, -96.72112],
                            decimal=self.prec)
        universe.trajectory[1]
        assert_almost_equal(universe.atoms.forces[42],
                            [-4.1959066, -31.31548, 22.663044],
                            decimal=self.prec)
        universe.trajectory[2]
        assert_almost_equal(universe.atoms.forces[10000],
                            [-41.43743, 83.35207, 62.94751],
                            decimal=self.prec)

    def test_h5md_dimensions(self, universe):
        universe.trajectory[0]
        assert_almost_equal(universe.trajectory.ts.dimensions,
                            [52.763, 52.763, 52.763, 90., 90., 90.],
                            decimal=self.prec)
        universe.trajectory[1]
        assert_almost_equal(universe.trajectory.ts.dimensions,
                            [52.807877, 52.807877, 52.807877, 90., 90., 90.],
                            decimal=self.prec)
        universe.trajectory[2]
        assert_almost_equal(universe.trajectory.ts.dimensions,
                            [52.839806, 52.839806, 52.839806, 90., 90., 90.],
                            decimal=self.prec)

    def test_h5md_data_step(self, universe):
        for ts, step in zip(universe.trajectory,
                           (0, 25000, 50000)):
            assert ts.data['step'] == step

    def test_rewind(self, universe):
        universe.trajectory[1]
        universe.trajectory.rewind()
        assert universe.trajectory.ts.frame == 0

    def test_next(self, universe):
        universe.trajectory.rewind()
        universe.trajectory.next()
        assert universe.trajectory.ts.frame == 1

    def test_jump_last_frame(self, universe):
        universe.trajectory[-1]
        assert universe.trajectory.ts.frame == 2

    @pytest.mark.parametrize("start, stop, step", ((0, 2, 1),
                                                   (1, 2, 1)))
    def test_slice(universe, start, stop, step):
        frames = [universe.trajectory.ts.frame
                  for ts in universe.trajectory[start:stop:step]]
        assert_array_equal(frames, np.arange(start, stop, step))

    @pytest.mark.parametrize("start, stop, step", ((0, 2, 1),
                                                   (1, 2, 1)))
    def test_slice(self, universe, start, stop, step):
        frames = [universe.trajectory.ts.frame
                  for ts in universe.trajectory[start:stop:step]]
        assert_array_equal(frames, np.arange(start, stop, step))

    @pytest.mark.parametrize("array_like", [list, np.array])
    def test_array_like(self, universe, array_like):
        array = array_like([0, 2])
        frames = [universe.trajectory.ts.frame
                  for ts in universe.trajectory[array]]
        assert_array_equal(frames, array)

    @pytest.mark.parametrize("indices", ([0, 1, 2, 1, 2, 2, 0]))
    def test_list_indices(self, universe, indices):
        frames = [universe.trajectory.ts.frame
                  for ts in universe.trajectory[indices]]
        assert_array_equal(frames, indices)

    @pytest.mark.parametrize('group, attr', (('position', 'positions'),
                                             ('velocity', 'velocities'),
                                             ('force', 'forces')))
    def test_no_group(self, h5md_file, outfile, attr, group):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g[f'particles/trajectory/{group}']
        u = mda.Universe(TPR_xvf, outfile)
        with pytest.raises(NoDataError):
            getattr(u.trajectory.ts, attr)

    @pytest.mark.parametrize('group', ('position', 'velocity', 'force'))
    def test_unknown_unit(self, h5md_file, outfile, group):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles'
                  '/trajectory'
                  f'/{group}'
                  '/value'].attrs['unit'] = 'random string'
        with pytest.raises(RuntimeError):
            u = mda.Universe(TPR_xvf, outfile)

    def test_unknown_time_unit(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles'
                  '/trajectory'
                  '/position'
                  '/time'].attrs['unit'] = 'random string'
        with pytest.raises(RuntimeError):
            u = mda.Universe(TPR_xvf, outfile)

    def test_length_unit_from_box(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/position']
        u = mda.Universe(TPR_xvf, outfile)
        assert u.trajectory.units['length'] == 'nm'

    @pytest.mark.parametrize('group', ('position', 'velocity', 'force'))
    def test_changing_n_atoms(self, h5md_file, outfile, group):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g[f'particles/trajectory/{group}/value'].resize((3, 10000, 3))
        with pytest.raises(ValueError):
            u = mda.Universe(TPR_xvf, outfile)

    def test_2D_box(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                new_box = np.ones(shape=(3, 2, 2))
                g['particles/trajectory/box'].attrs['dimension'] = 2
                del g['particles/trajectory/box/edges/value']
                g['particles/trajectory'
                  '/box/edges'].create_dataset('value', data=new_box)
        with pytest.raises(ValueError):
            u = mda.Universe(TPR_xvf, outfile)

    def test_no_box(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/box/edges']
        u = mda.Universe(TPR_xvf, outfile)
        assert u.trajectory.ts.dimensions is None

    def test_no_groups(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/position']
                del g['particles/trajectory/velocity']
                del g['particles/trajectory/force']
        with pytest.raises(NoDataError):
            u = mda.Universe(TPR_xvf, outfile)

    def test_no_convert_units(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                groups = ['position', 'velocity', 'force']
                for name in groups:
                    del g['particles/trajectory'][name]['value'].attrs['unit']
                del g['particles/trajectory/position/time'].attrs['unit']
        u = mda.Universe(TPR_xvf, outfile, convert_units=False)
        for unit in u.trajectory.units:
            assert u.trajectory.units[unit] is None

    def test_no_units(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                groups = ['position', 'velocity', 'force']
                for name in groups:
                    del g['particles/trajectory'][name]['value'].attrs['unit']
                del g['particles/trajectory/position/time'].attrs['unit']
                del g['particles/trajectory/box/edges/value'].attrs['unit']
        with pytest.raises(ValueError):
            u = mda.Universe(TPR_xvf, outfile)

    @pytest.mark.xfail(reason='Issue #2884')
    def test_open_filestream(self, h5md_file):
        with h5md_file as f:
            u = mda.Universe(TPR_xvf, h5md_file)

    def test_wrong_driver(self):
        with pytest.raises(ValueError):
            u = mda.Universe(TPR_xvf, H5MD_xvf,
                             driver='wrong_driver', comm="MPI.COMM_WORLD")

    def test_open_with_driver(self):
        u = mda.Universe(TPR_xvf, H5MD_xvf, driver="core")


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
class TestH5MDWriter_2(object):

    prec = 3

    @pytest.fixture()
    def universe(self):
        return mda.Universe(TPR_xvf, H5MD_xvf)

    @pytest.fixture()
    def Writer(self):
        return mda.coordinates.H5MD.H5MDWriter

    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + 'h5md-writer-test.h5md'

    @pytest.fixture()
    def outtop(self, tmpdir):
        return str(tmpdir) + 'h5md-writer-top.pdb'

    def test_no_atoms(self, universe, Writer, outfile):
        with pytest.raises(ValueError):
            with Writer(outfile, 0) as W:
                W.write(universe)

    @pytest.mark.parametrize('pos, vel, force', (
            (True, False, False),
            (True, True, False),
            (True, False, True),
            (True, True, True),
            (False, True, True),
            (False, False, True)
    ))
    def test_write_trajectory(self, universe, Writer, outfile,
                              pos, vel, force):
        with Writer(outfile,
                    universe.atoms.n_atoms,
                    positions=pos,
                    velocities=vel,
                    forces=force,
                    author='My Name',
                    author_email='my_email@asu.edu') as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile)

        # check the trajectory contents match reference universes
        for ts, ref_ts in zip(uw.trajectory, universe.trajectory):
            if pos:
                assert_almost_equal(ts._pos, ref_ts._pos, self.prec)
            else:
                with pytest.raises(NoDataError):
                    getattr(ts, 'positions')
            if vel:
                assert_almost_equal(ts._velocities, ref_ts._velocities,
                                    self.prec)
            else:
                with pytest.raises(NoDataError):
                    getattr(ts, 'velocities')
            if force:
                assert_almost_equal(ts._forces, ref_ts._forces, self.prec)
            else:
                with pytest.raises(NoDataError):
                    getattr(ts, 'forces')

    def test_write_AtomGroup(self, universe, outfile, outtop, Writer):
            """test to write H5MD from AtomGroup"""
            ca = universe.select_atoms("protein and name CA")
            ca.write(outtop)
            with Writer(outfile, n_atoms=ca.n_atoms) as W:
                for ts in universe.trajectory:
                    W.write(ca)

            uw = mda.Universe(outtop, outfile)
            caw = uw.atoms

            for orig_ts, written_ts in zip(universe.trajectory,
                                           uw.trajectory):
                assert_almost_equal(ca.positions, caw.positions, self.prec)
                assert_almost_equal(orig_ts.time, written_ts.time, self.prec)
                assert_almost_equal(written_ts.dimensions,
                                    orig_ts.dimensions,
                                    self.prec)

    @pytest.mark.parametrize('timeunit, lengthunit, velocityunit, forceunit', (
        ('fs', 'Angstrom', 'Angstrom/ps', 'kJ/(mol*Angstrom)'),
        ('s', 'pm', 'm/s', 'Newton'),
        ('ps', 'fm', 'Angstrom/fs', 'kcal/(mol*Angstrom)',),
        ('AKMA', 'nm', 'Angstrom/AKMA', 'kcal/(mol*Angstrom)')))
    def test_write_custom_units(self, universe, outfile, Writer,
                                timeunit, lengthunit,
                                velocityunit, forceunit):
        with Writer(outfile,
                    universe.atoms.n_atoms,
                    velocities=True,
                    forces=True,
                    lengthunit=lengthunit,
                    velocityunit=velocityunit,
                    forceunit=forceunit,
                    timeunit=timeunit) as W:
            for ts in universe.trajectory:
                W.write(universe)

        u = mda.Universe(TPR_xvf, outfile)
        for u_unit, custom_unit in zip (u.trajectory.units.values(),
                                       (timeunit, lengthunit,
                                       velocityunit, forceunit)):
            assert_equal(u_unit, custom_unit)

    @pytest.mark.parametrize('timeunit, lengthunit, velocityunit, forceunit', (
        ('imaginary time', None, None, None),
        (None, None, 'c', None),
        (None, None, None, 'HUGE FORCE',),
        (None, 'lightyear', None, None)))
    def test_write_bad_units(self, universe, outfile, Writer,
                             timeunit, lengthunit,
                             velocityunit, forceunit):
        with pytest.raises(RuntimeError):
            with Writer(outfile,
                        universe.atoms.n_atoms,
                        velocities=True,
                        forces=True,
                        lengthunit=lengthunit,
                        velocityunit=velocityunit,
                        forceunit=forceunit,
                        timeunit=timeunit) as W:
                for ts in universe.trajectory:
                    W.write(universe)

    #@pytest.mark.parametrize('chunks', (None, (1, 1000, 3)))
    def test_write_chunks(self, universe, outfile, Writer):
        chunks = (1, 1000, 3)
        with Writer(outfile,
                    universe.atoms.n_atoms,
                    chunks=chunks) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile)
        dset = uw.trajectory._particle_group['position/value']
        assert_equal(dset.chunks, chunks)

    @pytest.mark.parametrize('filter, opts', (('gzip', 1),
                                              ('gzip', 9),
                                              ('lzf', None)))
    def test_write_with_compression(self, universe,
                                    outfile, Writer,
                                    filter, opts):
        with Writer(outfile,
                    universe.atoms.n_atoms,
                    compression=filter,
                    compression_opts=opts) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile)
        dset = uw.trajectory._particle_group['position/value']
        assert_equal(dset.compression, filter)
        assert_equal(dset.compression_opts, opts)

    """def test_gaps(self, universe, Writer, outfile):
        ""Tests the writing and reading back of H5MDs with gaps in any of
        the coordinates/velocities properties.""
        with Writer(outfile, universe.atoms.n_atoms) as W:
            for ts in universe.trajectory:
                # Inset some gaps in the properties: coords every 4 steps, vels
                # every 2.
                if ts.frame == 1:
                    ts.has_positions = False
                if ts.frame == 2:
                    ts.has_velocities = False
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile)
        # check that the velocities are identical for each time step, except
        # for the gaps (that we must make sure to raise exceptions on).
        for orig_ts, written_ts in zip(universe.trajectory, uw.trajectory):
            if ts.frame == 1:
                assert_almost_equal(
                    written_ts.positions,
                    orig_ts.positions,
                    self.prec,
                    err_msg="coordinates mismatch "
                    "between original and written "
                    "trajectory at frame {} (orig) "
                    "vs {} (written)".format(orig_ts.frame, written_ts.frame))
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(written_ts, 'positions')

            if ts.frame == 2:
                assert_almost_equal(
                    written_ts.velocities,
                    orig_ts.velocities,
                    3,
                    err_msg="velocities mismatch "
                    "between original and written "
                    "trajectory at frame {} (orig) "
                    "vs {} (written)".format(orig_ts.frame, written_ts.frame))
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(written_ts, 'velocities')"""



    """def test_timestep_not_modified_by_writer(self, universe, Writer, outfile):
        trj = universe.trajectory
        ts = trj.ts

        trj[-1]  # last timestep (so that time != 0)
        x = ts._pos.copy()
        time = ts.time

        with Writer(outfile, trj.n_atoms) as W:
            # last timestep (so that time != 0) (say it again, just in case...)
            trj[-1]
            W.write(universe)

        assert_equal(
            ts._pos,
            x,
            err_msg="Positions in Timestep were modified by writer.")
        assert_equal(
            ts.time, time, err_msg="Time in Timestep was modified by writer.")"""
