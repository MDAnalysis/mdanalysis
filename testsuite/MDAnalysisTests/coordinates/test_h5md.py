import pytest
from numpy.testing import assert_almost_equal, assert_array_equal
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.H5MD import HAS_H5PY
from MDAnalysisTests.datafiles import H5MD_xvf, TPR_xvf, COORDINATES_TOPOLOGY, COORDINATES_H5MD
from MDAnalysisTests.coordinates.base import (MultiframeReaderTest,
                                              BaseReference, BaseWriterTest,
                                              assert_timestep_almost_equal)


h5py = pytest.importorskip('h5py')
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


"""The tests below test an example trajectory H5MD_xvf"""
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
                        [32.309906, 13.77798 , 14.372463],
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
                        [-6.8698354, 7.834235 , -8.114698],
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
    frames = [h5md_universe.trajectory.ts.frame for ts in h5md_universe.trajectory[start:stop:step]]
    assert_array_equal(frames, np.arange(start, stop, step))


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
@pytest.mark.parametrize("array_like", [list, np.array])
def test_array_like(h5md_universe, array_like):
    array = array_like([0, 2])
    frames = [h5md_universe.trajectory.ts.frame for ts in h5md_universe.trajectory[array]]
    assert_array_equal(frames, array)


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
@pytest.mark.parametrize("indices", ([0, 1, 2, 1, 2, 2, 0]))
def test_list_indices(h5md_universe, indices):
    frames = [h5md_universe.trajectory.ts.frame for ts in h5md_universe.trajectory[indices]]
    assert_array_equal(frames, indices)
