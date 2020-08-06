import pytest
from numpy.testing import assert_almost_equal, assert_array_equal
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.H5MD import HAS_H5PY
if HAS_H5PY:
    import h5py
from MDAnalysis.exceptions import NoDataError
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
    """Tests H5MDReader with MultiframeReaderTest.
    get_writer tests are marked with expected to fail since
    H5MDWriter is not implemented yet."""
    @staticmethod
    @pytest.fixture()
    def ref():
        return H5MDReference()

    @pytest.mark.xfail(reason='H5MD writer not implemented yet')
    def test_get_writer_1(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = 'test-writer.' + ref.ext
            with reader.Writer(outfile) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, reader.n_atoms)

    @pytest.mark.xfail(reason='H5MD writer not implemented yet')
    def test_get_writer_2(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = 'test-writer.' + ref.ext
            with reader.Writer(outfile, n_atoms=100) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, 100)


# The tests below test an example trajectory H5MD_xvf
@pytest.fixture
@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def h5md_universe():
    return mda.Universe(TPR_xvf, H5MD_xvf)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_h5md_n_frames(h5md_universe):
    assert len(h5md_universe.trajectory) == 3


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_h5md_positions(h5md_universe):
    # first timestep tests
    h5md_universe.trajectory[0]
    assert_almost_equal(h5md_universe.atoms.positions[0],
                        [32.309906, 13.77798, 14.372463],
                        decimal=6)
    assert_almost_equal(h5md_universe.atoms.positions[42],
                        [28.116928, 19.405945, 19.647358],
                        decimal=6)
    assert_almost_equal(h5md_universe.atoms.positions[10000],
                        [44.117805, 50.442093, 23.299038],
                        decimal=6)
    # second timestep tests
    h5md_universe.trajectory[1]
    assert_almost_equal(h5md_universe.atoms.positions[0],
                        [30.891968, 13.678971, 13.6000595],
                        decimal=6)
    assert_almost_equal(h5md_universe.atoms.positions[42],
                        [27.163246, 19.846561, 19.3582],
                        decimal=6)
    assert_almost_equal(h5md_universe.atoms.positions[10000],
                        [45.869278, 5.0342298, 25.460655],
                        decimal=6)
    # third timestep tests
    h5md_universe.trajectory[2]
    assert_almost_equal(h5md_universe.atoms.positions[0],
                        [31.276512, 13.89617, 15.015897],
                        decimal=6)
    assert_almost_equal(h5md_universe.atoms.positions[42],
                        [28.567991, 20.56532, 19.40814],
                        decimal=6)
    assert_almost_equal(h5md_universe.atoms.positions[10000],
                        [39.713223,  6.127234, 18.284992],
                        decimal=6)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_h5md_velocities(h5md_universe):
    h5md_universe.trajectory[0]
    assert_almost_equal(h5md_universe.atoms.velocities[0],
                        [-2.697732, 0.613568, 0.14334752],
                        decimal=6)
    h5md_universe.trajectory[1]
    assert_almost_equal(h5md_universe.atoms.velocities[42],
                        [-6.8698354, 7.834235, -8.114698],
                        decimal=6)
    h5md_universe.trajectory[2]
    assert_almost_equal(h5md_universe.atoms.velocities[10000],
                        [9.799492, 5.631466, 6.852126],
                        decimal=6)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_h5md_forces(h5md_universe):
    h5md_universe.trajectory[0]
    assert_almost_equal(h5md_universe.atoms.forces[0],
                        [20.071287, -155.2285, -96.72112],
                        decimal=5)
    h5md_universe.trajectory[1]
    assert_almost_equal(h5md_universe.atoms.forces[42],
                        [-4.1959066, -31.31548, 22.663044],
                        decimal=6)
    h5md_universe.trajectory[2]
    assert_almost_equal(h5md_universe.atoms.forces[10000],
                        [-41.43743, 83.35207, 62.94751],
                        decimal=5)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_h5md_dimensions(h5md_universe):
    # first timestep
    h5md_universe.trajectory[0]
    assert_almost_equal(h5md_universe.trajectory.ts.dimensions,
                        [52.763, 52.763, 52.763, 90., 90., 90.],
                        decimal=6)
    # second timestep
    h5md_universe.trajectory[1]
    assert_almost_equal(h5md_universe.trajectory.ts.dimensions,
                        [52.807877, 52.807877, 52.807877, 90., 90., 90.],
                        decimal=6)
    # third timestep
    h5md_universe.trajectory[2]
    assert_almost_equal(h5md_universe.trajectory.ts.dimensions,
                        [52.839806, 52.839806, 52.839806, 90., 90., 90.],
                        decimal=6)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_h5md_data_step(h5md_universe):
    h5md_universe.trajectory[0]
    assert h5md_universe.trajectory.ts.data['step'] == 0
    h5md_universe.trajectory[1]
    assert h5md_universe.trajectory.ts.data['step'] == 25000
    h5md_universe.trajectory[2]
    assert h5md_universe.trajectory.ts.data['step'] == 50000


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_rewind(h5md_universe):
    h5md_universe.trajectory.rewind()
    assert h5md_universe.trajectory.ts.frame == 0


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_next(h5md_universe):
    h5md_universe.trajectory.rewind()
    h5md_universe.trajectory.next()
    assert h5md_universe.trajectory.ts.frame == 1


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_jump_last_frame(h5md_universe):
    h5md_universe.trajectory[-1]
    assert h5md_universe.trajectory.ts.frame == 2


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
@pytest.mark.parametrize("start, stop, step", ((0, 2, 1),
                                               (1, 2, 1)))
def test_slice(h5md_universe, start, stop, step):
    frames = [h5md_universe.trajectory.ts.frame
              for ts in h5md_universe.trajectory[start:stop:step]]
    assert_array_equal(frames, np.arange(start, stop, step))


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
@pytest.mark.parametrize("array_like", [list, np.array])
def test_array_like(h5md_universe, array_like):
    array = array_like([0, 2])
    frames = [h5md_universe.trajectory.ts.frame
              for ts in h5md_universe.trajectory[array]]
    assert_array_equal(frames, array)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
@pytest.mark.parametrize("indices", ([0, 1, 2, 1, 2, 2, 0]))
def test_list_indices(h5md_universe, indices):
    frames = [h5md_universe.trajectory.ts.frame
              for ts in h5md_universe.trajectory[indices]]
    assert_array_equal(frames, indices)


@pytest.fixture
@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def ref():
    return H5MDReference()


@pytest.fixture
@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def h5md_file():
    return h5py.File(H5MD_xvf, 'r')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_no_positions(h5md_file, ref, tmpdir):
    outfile = 'test_no_positions' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/position']

        u = mda.Universe(TPR_xvf, outfile, format='H5MD')
        with pytest.raises(NoDataError):
            u.trajectory.ts.positions


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_no_velocities(h5md_file, ref, tmpdir):
    outfile = 'test_no_velocities' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/velocity']
        u = mda.Universe(TPR_xvf, outfile, format='H5MD')
        with pytest.raises(NoDataError):
            u.trajectory.ts.velocities


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_no_forces(h5md_file, ref, tmpdir):
    outfile = 'test_no_forces' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/force']
        u = mda.Universe(TPR_xvf, outfile, format='H5MD')
        with pytest.raises(NoDataError):
            u.trajectory.ts.forces


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_unknown_position_unit(h5md_file, ref, tmpdir):
    outfile = 'test_unknown_position_unit' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles'
                  '/trajectory'
                  '/position'
                  '/value'].attrs['unit'] = 'random string'
        with pytest.raises(RuntimeError):
            u = mda.Universe(TPR_xvf, outfile, format='H5MD')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_length_unit_from_box(h5md_file, ref, tmpdir):
    outfile = 'test_length_unit_from_box' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/position']
        u = mda.Universe(TPR_xvf, outfile, format='H5MD')
        assert u.trajectory.units['length'] == 'Angstrom'


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_unknown_time_unit(h5md_file, ref, tmpdir):
    outfile = 'test_unknown_time_unit' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles'
                  '/trajectory'
                  '/position'
                  '/time'].attrs['unit'] = 'random string'
        with pytest.raises(RuntimeError):
            u = mda.Universe(TPR_xvf, outfile, format='H5MD')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_unknown_velocity_unit(h5md_file, ref, tmpdir):
    outfile = 'test_unknown_velocity_unit' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles'
                  '/trajectory'
                  '/velocity'
                  '/value'].attrs['unit'] = 'random string'
        with pytest.raises(RuntimeError):
            u = mda.Universe(TPR_xvf, outfile, format='H5MD')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_unknown_force_unit(h5md_file, ref, tmpdir):
    outfile = 'test_unknown_force_unit' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles'
                  '/trajectory'
                  '/force'
                  '/value'].attrs['unit'] = 'random string'
        with pytest.raises(RuntimeError):
            u = mda.Universe(TPR_xvf, outfile, format='H5MD')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_changing_n_atoms1(h5md_file, ref, tmpdir):
    outfile = 'test_changing_n_atoms1' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles/trajectory/position/value'].resize((3, 10000, 3))
        with pytest.raises(ValueError):
            u = mda.Universe(TPR_xvf, outfile, format='H5MD')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_changing_n_atoms2(h5md_file, ref, tmpdir):
    outfile = 'test_changing_n_atoms2' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles/trajectory/velocity/value'].resize((3, 10000, 3))
        with pytest.raises(ValueError):
            u = mda.Universe(TPR_xvf, outfile, format='H5MD')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_changing_n_atoms3(h5md_file, ref, tmpdir):
    outfile = 'test_changing_n_atoms3' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                g['particles/trajectory/force/value'].resize((3, 10000, 3))
        with pytest.raises(ValueError):
            u = mda.Universe(TPR_xvf, outfile, format='H5MD')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_2D_box(h5md_file, ref, tmpdir):
    outfile = 'test_2D_box' + ref.ext
    with tmpdir.as_cwd():
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
            u = mda.Universe(TPR_xvf, outfile, format='H5MD')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_no_box(h5md_file, ref, tmpdir):
    outfile = 'test_no_box' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/box/edges']
        u = mda.Universe(TPR_xvf, outfile, format='H5MD')
        assert u.trajectory.ts.dimensions is None


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_no_groups(h5md_file, ref, tmpdir):
    outfile = 'test_no_groups' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                del g['particles/trajectory/position']
                del g['particles/trajectory/velocity']
                del g['particles/trajectory/force']
        with pytest.raises(NoDataError):
            u = mda.Universe(TPR_xvf, outfile, format='H5MD')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
@pytest.mark.xfail(reason='Issue #2884')
def test_open_filestream(h5md_file):
    with h5md_file as f:
        u = mda.Universe(TPR_xvf, h5md_file)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_wrong_driver():
    with pytest.raises(ValueError):
        u = mda.Universe(TPR_xvf, H5MD_xvf,
                         driver='wrong_driver', comm="MPI.COMM_WORLD")


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_no_convert_units(h5md_file, ref, tmpdir):
    outfile = 'test_no_convert_units' + ref.ext
    with tmpdir.as_cwd():
        with h5md_file as f:
            with h5py.File(outfile, 'w') as g:
                f.copy(source='particles', dest=g)
                f.copy(source='h5md', dest=g)
                groups = ['position', 'velocity', 'force']
                for name in groups:
                    del g['particles/trajectory'][name]['value'].attrs['unit']
                del g['particles/trajectory/position/time'].attrs['unit']
        u = mda.Universe(TPR_xvf, outfile, convert_units=False, format="H5MD")
        for unit in u.trajectory.units:
            assert u.trajectory.units[unit] is None


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_no_units(h5md_file, ref, tmpdir):
    outfile = 'test_no_units' + ref.ext
    with tmpdir.as_cwd():
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
            u = mda.Universe(TPR_xvf, outfile, format="H5MD")


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_open_with_driver():
    u = mda.Universe(TPR_xvf, H5MD_xvf, driver="core")
