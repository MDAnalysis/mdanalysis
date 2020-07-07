import pytest
from numpy.testing import assert_almost_equal, assert_array_equal
import numpy as np
import MDAnalysis as mda
from MDAnalysisTests.datafiles import H5MD_xvf, TPR_xvf
from MDAnalysis.coordinates import H5MD


@pytest.fixture
def h5md_universe():
    return mda.Universe(TPR_xvf, H5MD_xvf)

def test_h5md_n_frames(h5md_universe):
    assert len(h5md_universe.trajectory) == 3

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

def test_h5md_dimensions(h5md_universe):
    # first timestep
    ts = h5md_universe.trajectory[0]
    assert_almost_equal(ts.dimensions,
                        [52.763, 52.763, 52.763, 90., 90., 90.],
                        decimal=6)
    # second timestep
    ts = h5md_universe.trajectory[1]
    assert_almost_equal(ts.dimensions,
                        [52.807877, 52.807877, 52.807877, 90., 90., 90.],
                        decimal=6)
    # third timestep
    ts = h5md_universe.trajectory[2]
    assert_almost_equal(ts.dimensions,
                        [52.839806, 52.839806, 52.839806, 90., 90., 90.],
                        decimal=6)

def test_h5md_data_step(h5md_universe):
    assert h5md_universe.trajectory[0].data['step'] == 0
    assert h5md_universe.trajectory[1].data['step'] == 1
    assert h5md_universe.trajectory[2].data['step'] == 2

def test_rewind(h5md_universe):
    h5md_universe.trajectory.rewind()
    assert h5md_universe.trajectory.ts.frame == 0

def test_next(h5md_universe):
    h5md_universe.trajectory.rewind()
    h5md_universe.trajectory.next()
    assert h5md_universe.trajectory.ts.frame == 1

def test_jump_last_frame(h5md_universe):
    h5md_universe.trajectory[-1]
    assert h5md_universe.trajectory.ts.frame == 2

@pytest.mark.parametrize("start, stop, step", ((0, 2, 1),
                                               (1, 2, 1)))
def test_slice(h5md_universe, start, stop, step):
    frames = [h5md_universe.trajectory.ts.frame for ts in h5md_universe.trajectory[start:stop:step]]
    assert_array_equal(frames, np.arange(start, stop, step))


@pytest.mark.parametrize("array_like", [list, np.array])
def test_array_like(h5md_universe, array_like):
    array = array_like([0, 2])
    frames = [h5md_universe.trajectory.ts.frame for ts in h5md_universe.trajectory[array]]
    assert_array_equal(frames, array)

@pytest.mark.parametrize("indices", ([0, 1, 2, 1, 2, 2, 0]))
def test_list_indices(h5md_universe, indices):
    frames = [h5md_universe.trajectory.ts.frame for ts in h5md_universe.trajectory[indices]]
    assert_array_equal(frames, indices)
