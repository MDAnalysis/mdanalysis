import pytest
from numpy.testing import assert_almost_equal, assert_equal, assert_allclose
import numpy as np
import sys
import os
import MDAnalysis as mda
from MDAnalysis.coordinates.H5MD import HAS_H5PY

if HAS_H5PY:
    import h5py
    from MDAnalysis.coordinates.H5MD import H5MDReader
from MDAnalysis.exceptions import NoDataError
from MDAnalysisTests.datafiles import (
    H5MD_xvf,
    TPR_xvf,
    TRR_xvf,
    COORDINATES_TOPOLOGY,
    COORDINATES_H5MD,
    H5MD_energy,
)
from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest,
    BaseReference,
    BaseWriterTest,
)


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
class TestH5MDReaderBaseAPI(MultiframeReaderTest):
    """Tests H5MDReader with with synthetic trajectory."""
    @staticmethod
    @pytest.fixture()
    def ref():
        return H5MDReference()

    def test_get_writer_1(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = 'test-writer.' + ref.ext
            with reader.Writer(outfile) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, reader.n_atoms)

    def test_get_writer_2(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = 'test-writer.' + ref.ext
            with reader.Writer(outfile, n_atoms=100) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, 100)

    def test_copying(self, ref, reader):
        # Issue #3664 - test not done in test_copying due to dependencies
        original = mda.coordinates.H5MD.H5MDReader(
                ref.trajectory, convert_units=False, dt=2,
                time_offset=10, foo="bar")
        copy = original.copy()

        assert original.format not in ('MEMORY', 'CHAIN')
        assert original.convert_units is False
        assert copy.convert_units is False
        assert original._ts_kwargs['time_offset'] == 10
        assert copy._ts_kwargs['time_offset'] == 10
        assert original._ts_kwargs['dt'] == 2
        assert copy._ts_kwargs['dt'] == 2

        assert original.ts.data['time_offset'] == 10
        assert copy.ts.data['time_offset'] == 10

        assert original.ts.data['dt'] == 2
        assert copy.ts.data['dt'] == 2

        assert copy._kwargs['foo'] == 'bar'

        # check coordinates
        assert original.ts.frame == copy.ts.frame
        assert_allclose(original.ts.positions, copy.ts.positions)

        original.next()
        copy.next()

        assert original.ts.frame == copy.ts.frame
        assert_allclose(original.ts.positions, copy.ts.positions)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
class TestH5MDWriterBaseAPI(BaseWriterTest):
    """Tests H5MDWriter base API with synthetic trajectory"""
    @staticmethod
    @pytest.fixture()
    def ref():
        return H5MDReference()

    def test_write_trajectory_atomgroup(self, ref, reader, universe, tmpdir):
        outfile = 'write-atoms-test.' + ref.ext
        with tmpdir.as_cwd():
            with ref.writer(outfile, universe.atoms.n_atoms,
                            velocities=True, forces=True) as w:
                for ts in universe.trajectory:
                    w.write(universe.atoms)
            self._check_copy(outfile, ref, reader)

    @pytest.mark.xfail((os.name == 'nt' and sys.maxsize <= 2**32),
                       reason="occasional fail on 32-bit windows")
    def test_write_trajectory_universe(self, ref, reader, universe, tmpdir):
        outfile = 'write-uni-test.' + ref.ext
        with tmpdir.as_cwd():
            with ref.writer(outfile, universe.atoms.n_atoms,
                            velocities=True, forces=True) as w:
                for ts in universe.trajectory:
                    w.write(universe)
            self._check_copy(outfile, ref, reader)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
class TestH5MDReaderWithRealTrajectory(object):

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
        for ts, step in zip(universe.trajectory, (0, 25000, 50000)):
            assert_equal(ts.data['step'], step)

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
    def test_slice(self, universe, start, stop, step):
        frames = [universe.trajectory.ts.frame
                  for ts in universe.trajectory[start:stop:step]]
        assert_equal(frames, np.arange(start, stop, step))

    @pytest.mark.parametrize("array_like", [list, np.array])
    def test_array_like(self, universe, array_like):
        array = array_like([0, 2])
        frames = [universe.trajectory.ts.frame
                  for ts in universe.trajectory[array]]
        assert_equal(frames, array)

    def test_list_indices(self, universe):
        indices = [0, 1, 2, 1, 2, 2, 0]
        frames = [universe.trajectory.ts.frame
                  for ts in universe.trajectory[indices]]
        assert_equal(frames, indices)

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
        with pytest.raises(NoDataError,
                           match="This Timestep has no"):
            getattr(u.trajectory.ts, attr)

    @pytest.mark.parametrize('dset',
        ('position/value', 'position/time', 'velocity/value', 'force/value'))
    def test_unknown_unit(self, h5md_file, outfile, dset):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles'
                  '/trajectory'
                  f'/{dset}'].attrs['unit'] = 'random string'
        with pytest.raises(RuntimeError,
                           match=" is not recognized by H5MDReader."):
            u = mda.Universe(TPR_xvf, outfile)

    def test_length_unit_from_box(self, h5md_file, universe, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/position']
        ref_u = universe
        uw = mda.Universe(TPR_xvf, outfile)
        assert_equal(ref_u.trajectory.units['length'],
                     uw.trajectory.units['length'])
        for ref_ts, new_ts in zip(ref_u.trajectory, uw.trajectory):
            assert_equal(ref_ts.dimensions, new_ts.dimensions)
            assert_equal(ref_ts.triclinic_dimensions,
                         new_ts.triclinic_dimensions)

    @pytest.mark.parametrize('group', ('position', 'velocity', 'force'))
    def test_changing_n_atoms(self, h5md_file, outfile, group):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g[f'particles/trajectory/{group}/value'].resize((3, 10000, 3))
        with pytest.raises(ValueError,
                           match=" of either the postion, velocity, or force"):
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
        with pytest.raises(ValueError,
                           match="MDAnalysis only supports 3-dimensional"):
            u = mda.Universe(TPR_xvf, outfile)

    def test_box_vector(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                vector = [1, 2, 3]
                del g['particles/trajectory/box/edges']
                g['particles/trajectory/box/edges/value'] = [vector, vector, vector]
        u = mda.Universe(TPR_xvf, outfile)
        # values in vector are conveted from nm -> Angstrom
        assert_equal(u.trajectory.ts.dimensions, [10, 20, 30, 90, 90, 90])

    def test_no_box(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/box/edges']
        u = mda.Universe(TPR_xvf, outfile)
        assert_equal(u.trajectory.ts.dimensions, None)

    def test_no_groups(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/position']
                del g['particles/trajectory/velocity']
                del g['particles/trajectory/force']
        with pytest.raises(NoDataError,
                           match="Provide at least a position, velocity"):
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
            assert_equal(u.trajectory.units[unit], None)

    def test_no_units_but_convert_units_true_error(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                groups = ['position', 'velocity', 'force']
                for name in groups:
                    del g['particles/trajectory'][name]['value'].attrs['unit']
                del g['particles/trajectory/position/time'].attrs['unit']
                del g['particles/trajectory/box/edges/value'].attrs['unit']
        with pytest.raises(ValueError,
                           match="H5MD file must have readable units if"):
            u = mda.Universe(TPR_xvf, outfile, convert_units=True)

    @pytest.mark.xfail(reason='Issue #2884')
    def test_open_filestream(self, h5md_file, universe):
        with h5md_file as f:
            u = mda.Universe(TPR_xvf, h5md_file)
            for ts1, ts2 in zip(universe.trajectory, u.trajectory):
                assert_equal(ts1.positions, ts2.positions)
                assert_equal(ts1.velocities, ts2.velocities)
                assert_equal(ts1.forces, ts2.forces)

    def test_wrong_driver(self):
        with pytest.raises(ValueError,
                           match="If MPI communicator object is used to open"):
            u = mda.Universe(TPR_xvf, H5MD_xvf,
                             driver='wrong_driver',
                             comm="mock MPI.COMM_WORLD")

    def test_open_with_driver(self):
        u = mda.Universe(TPR_xvf, H5MD_xvf, driver="core")
        assert_equal(u.trajectory._file.driver, "core")

    @pytest.mark.parametrize('group1, group2', (('velocity', 'force'),
                                                ('position', 'force'),
                                                ('position', 'velocity')))
    def test_parse_n_atoms(self, h5md_file, outfile, group1, group2):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                traj_group = g['particles/trajectory']
                del traj_group[group1]
                del traj_group[group2]
                for dset in ('position/value', 'velocity/value',
                             'force/value'):
                    try:
                        n_atoms_in_dset = traj_group[dset].shape[1]
                        break
                    except KeyError:
                        continue

        u = mda.Universe(outfile)
        assert_equal(u.atoms.n_atoms, n_atoms_in_dset)

    def test_parse_n_atoms_error(self, h5md_file, outfile):
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                traj_group = g['particles/trajectory']
                del traj_group['position']
                del traj_group['velocity']
                del traj_group['force']

        errmsg = "Could not construct minimal topology"
        with pytest.raises(ValueError, match=errmsg):
            u = mda.Universe(outfile)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
class TestH5MDWriterWithRealTrajectory(object):

    prec = 3

    @pytest.fixture()
    def universe(self):
        return mda.Universe(TPR_xvf, H5MD_xvf)

    @pytest.fixture()
    def universe_no_units(self):
        u = mda.Universe(TPR_xvf, H5MD_xvf, convert_units=False)
        u.trajectory.units['time'] = None
        u.trajectory.units['length'] = None
        u.trajectory.units['velocity'] = None
        u.trajectory.units['force'] = None
        return u

    @pytest.fixture()
    def Writer(self):
        return mda.coordinates.H5MD.H5MDWriter

    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + 'h5md-writer-test.h5md'

    @pytest.fixture()
    def outtop(self, tmpdir):
        return str(tmpdir) + 'h5md-writer-top.pdb'

    @pytest.mark.parametrize('scalar, error, match',
        ((0, ValueError, "H5MDWriter: no atoms in output trajectory"),
        (0.5, IOError, "H5MDWriter: Timestep does not have")))
    def test_n_atoms_errors(self, universe, Writer, outfile,
                            scalar, error, match):
        n_atoms = universe.atoms.n_atoms * scalar
        with pytest.raises(error, match=match):
            with Writer(outfile, n_atoms) as W:
                W.write(universe)

    def test_chunk_error(self, universe, Writer, outfile):
        n_atoms = universe.atoms.n_atoms
        err = "H5MDWriter must know how many frames will be "
        with pytest.raises(ValueError, match=err):
            with Writer(outfile, n_atoms, chunks=False) as W:
                for ts in universe.trajectory:
                    W.write(universe)

    @pytest.mark.parametrize('dimensions', (None, 0))
    def test_no_dimensions(self, universe, Writer, outfile, dimensions):
        with Writer(outfile, universe.atoms.n_atoms) as W:
            for ts in universe.trajectory:
                ts.dimensions = dimensions
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile)
        box = uw.trajectory._particle_group['box']
        assert 'edges' not in box
        assert_equal(3*['none'], box.attrs['boundary'])
        assert_equal(3, box.attrs['dimension'])

    def test_step_not_monotonic(self, universe, Writer, outfile):
        with pytest.raises(ValueError,
                           match="The H5MD standard dictates that the step "):
            with Writer(outfile, universe.atoms.n_atoms) as W:
                for ts in universe.trajectory[[0, 1, 2, 1]]:
                    W.write(universe)

        with pytest.raises(ValueError,
                           match="The H5MD standard dictates that the step "):
            with Writer(outfile, universe.atoms.n_atoms) as W:
                for ts in universe.trajectory:
                    if ts.frame == 2:
                        ts.data['step'] = 0
                    W.write(universe)

    def test_step_from_frame(self, universe, Writer, outfile):
        with Writer(outfile, universe.atoms.n_atoms) as W:
            for ts in universe.trajectory:
                del ts.data['step']
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile)
        steps = [ts.data['step'] for ts in uw.trajectory]
        frames = [ts.frame for ts in universe.trajectory]
        for step, frame in zip(steps, frames):
            assert_equal(step, frame)

    def test_has_property(self, universe, Writer, outfile):
        with Writer(outfile, universe.atoms.n_atoms) as W:
            W.write(universe)
            # make sure property is pulled from _has dict
            assert W.has_positions == W._has['position']
            assert W.has_velocities == W._has['velocity']
            assert W.has_forces == W._has['force']
            # make sure the values are correct
            assert W.has_positions is True
            assert W.has_velocities is True
            assert W.has_forces is True

    @pytest.mark.parametrize('pos, vel, force', (
            (True, False, False),
            (True, True, False),
            (True, False, True),
            (True, True, True),
            (False, True, True),
            (False, False, True),
            (False, False, False)))
    def test_write_trajectory(self, universe, Writer, outfile,
                              pos, vel, force):
        try:
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
                assert_almost_equal(ts.dimensions, ref_ts.dimensions, self.prec)
                if pos:
                    assert_almost_equal(ts._pos, ref_ts._pos, self.prec)
                else:
                    with pytest.raises(NoDataError,
                                       match="This Timestep has no"):
                        getattr(ts, 'positions')
                if vel:
                    assert_almost_equal(ts._velocities, ref_ts._velocities,
                                        self.prec)
                else:
                    with pytest.raises(NoDataError,
                                       match="This Timestep has no"):
                        getattr(ts, 'velocities')
                if force:
                    assert_almost_equal(ts._forces, ref_ts._forces, self.prec)
                else:
                    with pytest.raises(NoDataError,
                                       match="This Timestep has no"):
                        getattr(ts, 'forces')

        # when (False, False, False)
        except(ValueError):
            with pytest.raises(ValueError,
                               match="At least one of positions, velocities"):
                with Writer(outfile,
                            universe.atoms.n_atoms,
                            positions=pos,
                            velocities=vel,
                            forces=force,
                            author='My Name',
                            author_email='my_email@asu.edu') as W:
                    for ts in universe.trajectory:
                        W.write(universe)

    def test_write_AtomGroup_with(self, universe, outfile, outtop, Writer):
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

    @pytest.mark.parametrize('frames, n_frames', ((None, 1),
                                                  ('all', 3)))
    def test_ag_write(self, universe, outfile, outtop,
                      Writer, frames, n_frames):
        """test to write with ag.write()"""
        ca = universe.select_atoms("protein and name CA")
        ca.write(outtop)

        ca.write(outfile, frames=frames, format='h5md')

        uw = mda.Universe(outtop, outfile)
        caw = uw.atoms

        assert_equal(n_frames, len(uw.trajectory))
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
                    lengthunit=lengthunit,
                    velocityunit=velocityunit,
                    forceunit=forceunit,
                    timeunit=timeunit) as W:
            for ts in universe.trajectory:
                W.write(universe)

        u = mda.Universe(TPR_xvf, outfile)
        for u_unit, custom_unit in zip(u.trajectory.units.values(),
                                       (timeunit, lengthunit,
                                       velocityunit, forceunit)):
            assert_equal(u_unit, custom_unit)

    @pytest.mark.parametrize('timeunit, lengthunit, velocityunit, forceunit', (
        ('imaginary time', None, None, None),
        (None, None, 'c', None),
        (None, None, None, 'HUGE FORCE',),
        (None, 'lightyear', None, None),))
    def test_write_bad_units(self, universe, outfile, Writer,
                             timeunit, lengthunit,
                             velocityunit, forceunit):
        with pytest.raises(ValueError, match=" is not a unit recognized by"):
            with Writer(outfile,
                        universe.atoms.n_atoms,
                        lengthunit=lengthunit,
                        velocityunit=velocityunit,
                        forceunit=forceunit,
                        timeunit=timeunit) as W:
                for ts in universe.trajectory:
                    W.write(universe)

    def test_no_units_w_convert_true(self, universe_no_units, outfile, Writer):
        # no units + convert_units = ValueError
        with pytest.raises(ValueError, match="The trajectory has no units,"):
            with Writer(outfile,
                        universe_no_units.atoms.n_atoms) as W:
                for ts in universe_no_units.trajectory:
                    W.write(universe_no_units)

    def test_no_units_w_convert_false(self, universe_no_units,
                                      outfile, Writer):
        with Writer(outfile,
                    universe_no_units.atoms.n_atoms,
                    convert_units=False) as W:
            for ts in universe_no_units.trajectory:
                W.write(universe_no_units)

            uw = mda.Universe(TPR_xvf, outfile, convert_units=False)
            for unit in uw.trajectory.units.values():
                assert_equal(unit, None)

    @pytest.mark.parametrize('convert_units', (True, False))
    def test_convert_units(self, universe, outfile, Writer,
                           convert_units):
        with Writer(outfile,
                    universe.atoms.n_atoms,
                    convert_units=convert_units) as W:
            for ts in universe.trajectory:
                W.write(universe)

        ref_units = universe.trajectory.units.items()
        uw = mda.Universe(TPR_xvf, outfile)
        uw_units = uw.trajectory.units.items()
        for u1, u2 in zip(ref_units, uw_units):
            assert_equal(u1, u2)

    @pytest.mark.parametrize('chunks', ((3, 1000, 1),
                                        (1, 1000, 3),
                                        (100, 100, 3)))
    def test_write_chunks(self, universe, outfile, Writer, chunks):
        with Writer(outfile,
                    universe.atoms.n_atoms,
                    chunks=chunks) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile)
        for dset in (uw.trajectory._particle_group['position/value'],
                     uw.trajectory._particle_group['velocity/value'],
                     uw.trajectory._particle_group['force/value']):
            assert_equal(dset.chunks, chunks)

        for ts1, ts2 in zip(universe.trajectory, uw.trajectory):
            assert_equal(ts1.positions, ts2.positions)
            assert_equal(ts1.velocities, ts2.velocities)
            assert_equal(ts1.forces, ts2.forces)

    def test_write_chunks_with_nframes(self, universe, outfile, Writer):
        n_atoms = universe.atoms.n_atoms
        n_frames = universe.trajectory.n_frames
        with Writer(outfile,
                    n_atoms=n_atoms,
                    n_frames=n_frames) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile)
        for dset in (uw.trajectory._particle_group['position/value'],
                     uw.trajectory._particle_group['velocity/value'],
                     uw.trajectory._particle_group['force/value']):
            assert_equal(dset.chunks, (1, n_atoms, 3))

        for ts1, ts2 in zip(universe.trajectory, uw.trajectory):
            assert_equal(ts1.positions, ts2.positions)
            assert_equal(ts1.velocities, ts2.velocities)
            assert_equal(ts1.forces, ts2.forces)

    def test_write_contiguous1(self, universe, Writer, outfile):
        n_atoms = universe.atoms.n_atoms
        n_frames = len(universe.trajectory)
        with Writer(outfile,
                    n_atoms=n_atoms,
                    n_frames=n_frames,
                    chunks=False) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile)
        for dset in (uw.trajectory._particle_group['position/value'],
                     uw.trajectory._particle_group['velocity/value'],
                     uw.trajectory._particle_group['force/value']):
            assert_equal(dset.chunks, None)

    def test_write_contiguous2(self, universe, Writer, outfile):
        ag = universe.select_atoms('all')
        n_frames = len(ag.universe.trajectory)
        ag.write(outfile, frames='all', n_frames=n_frames, chunks=False)

        uw = mda.Universe(TPR_xvf, outfile)
        for dset in (uw.trajectory._particle_group['position/value'],
                     uw.trajectory._particle_group['velocity/value'],
                     uw.trajectory._particle_group['force/value']):
            assert_equal(dset.chunks, None)

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

    @pytest.mark.xfail(os.name == 'nt',
                       reason="occasional PermissionError on windows")
    @pytest.mark.parametrize('driver', ('core', 'stdio'))
    def test_write_with_drivers(self, universe, outfile, Writer, driver):
        with Writer(outfile,
                    universe.atoms.n_atoms,
                    driver=driver) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(TPR_xvf, outfile, driver=driver)
        file = uw.trajectory._file
        assert_equal(file.driver, driver)

    def test_parallel_disabled(self, universe, Writer, outfile,
                               driver='mpio'):
        with pytest.raises(ValueError,
                           match="H5MDWriter: parallel writing with MPI I/O "):
            with Writer(outfile,
                        universe.atoms.n_atoms,
                        driver=driver) as W:
                for ts in universe.trajectory:
                    W.write(universe)

    def test_timestep_not_modified_by_writer(self, universe, Writer, outfile):
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
            ts.time, time, err_msg="Time in Timestep was modified by writer.")


class TestH5PYNotInstalled(object):
    """Tests RuntimeErrors when h5py not installed"""

    @pytest.fixture(autouse=True)
    def block_h5py(self, monkeypatch):
        monkeypatch.setattr(
            sys.modules['MDAnalysis.coordinates.H5MD'], 'HAS_H5PY', False)

    @pytest.fixture()
    def Writer(self):
        return mda.coordinates.H5MD.H5MDWriter

    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + 'h5md-writer-test.h5md'

    def test_reader_no_h5py(self):
        with pytest.raises(RuntimeError, match="Please install h5py"):
            u = mda.Universe(TPR_xvf, H5MD_xvf)

    def test_writer_no_h5py(self, Writer, outfile):
        u = mda.Universe(TPR_xvf, TRR_xvf)
        with pytest.raises(RuntimeError,
                           match="H5MDWriter: Please install h5py"):
            with Writer(outfile,
                        u.atoms.n_atoms) as W:
                for ts in u.trajectory:
                    W.write(universe)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
class TestH5MDReaderWithObservables(object):
    """Read H5MD file with 'observables/atoms/energy'."""

    prec = 3
    ext = 'h5md'

    def test_read_h5md_issue4598(self):
        """Read a H5MD file with observables.

        The reader will ignore the 'observables/atoms/energy'.
        """

        u = mda.Universe.empty(n_atoms=108, trajectory=True)
        reader = H5MDReader(H5MD_energy)
        u.trajectory = reader
        for ts in u.trajectory:
            assert "atoms/energy" in ts.data
