import numpy as np
import pytest
from numpy.testing import assert_allclose

import MDAnalysis as mda
from MDAnalysis.transformations import NoJump, wrap
from MDAnalysisTests import datafiles as data


@pytest.fixture()
def nojump_universes_fromfile():
    '''
    Create the universe objects for the tests.
    '''
    u = mda.Universe(data.PSF_TRICLINIC, data.DCD_TRICLINIC)
    transformation = NoJump()
    u.trajectory.add_transformations(transformation)
    return u


@pytest.fixture()
def nojump_universe():
    """
    Create the universe objects for the tests.
    """
    u = mda.Universe.empty(1, trajectory=True)
    coordinates = np.empty((100, u.atoms.n_atoms, 3))  # number of frames
    coordinates[::3, 0] = 0 * np.ones(3) / 3
    coordinates[1::3, 0] = 1 * np.ones(3) / 3
    coordinates[2::3, 0] = 2 * np.ones(3) / 3
    u.load_new(coordinates, order="fac")
    return u


@pytest.fixture()
def nojump_constantvel_universe():
    """
    Create the universe objects for the tests.
    """
    Natom = 1
    Nframe = 100
    coordinates = np.empty((Nframe, Natom, 3))  # number of frames
    coordinates[:, 0, 0] = np.linspace(0, 45, Nframe)
    coordinates[:, 0, 1] = np.linspace(0, 15, Nframe)
    coordinates[:, 0, 2] = np.linspace(0, 10, Nframe)
    reference = mda.Universe.empty(Natom, trajectory=True)
    reference.load_new(coordinates, order="fac")
    return reference


@pytest.fixture
def nojump_universe_npt_2nd_frame():
    """
    Create a Universe in which a single atom jumps across the periodic boundary in the
    x-dimensions at the second frame.

    Unwrapped coordinates should all be 97.5.
    """
    n_atoms = 1
    n_frames = 4
    u = mda.Universe.empty(n_atoms, trajectory=True)
    coordinates = np.empty((n_frames, u.atoms.n_atoms, 3))
    coordinates[0] = [97.5, 50.0, 50.0]
    coordinates[1] = [2.5, 50.0, 50.0]
    coordinates[2] = [2.5, 50.0, 50.0]
    coordinates[3] = [2.5, 50.0, 50.0]
    u.load_new(coordinates, order="fac")
    return u


@pytest.fixture
def nojump_universe_npt_3rd_frame():
    """
    Create a Universe in which a single atom jumps across the periodic boundary in the
    x-dimensions at the third frame.

    Unwrapped coordinates should all be 97.5.
    """
    n_atoms = 1
    n_frames = 4
    u = mda.Universe.empty(n_atoms, trajectory=True)
    coordinates = np.empty((n_frames, u.atoms.n_atoms, 3))
    coordinates[0] = [97.5, 50.0, 50.0]
    coordinates[1] = [97.5, 50.0, 50.0]
    coordinates[2] = [2.5, 50.0, 50.0]
    coordinates[3] = [2.5, 50.0, 50.0]
    u.load_new(coordinates, order="fac")
    return u


@pytest.fixture(scope="module")
def nojump_universe_npt_2nd_frame_from_file(tmp_path_factory):
    """
    Write the `nojump_universe_npt_2nd_frame` fixture to file, read it in and
    return the Universe.

    Used for testing that coordinates can be unwrapped correctly when iterating
    over the trajectory multiple times.

    We can't use an in-memory trajectory to test this because the
    transformation would only be applied once.

    Note, we use `tmp_path_factory` because this fixture requies `module` scope
    so we can read the file after the fixture has been created, and
    `tmp_path` has function-level scope.
    """
    n_atoms = 1
    n_frames = 4
    u = mda.Universe.empty(n_atoms, trajectory=True)
    coordinates = np.empty((n_frames, u.atoms.n_atoms, 3))
    coordinates[0] = [97.5, 50.0, 50.0]
    coordinates[1] = [2.5, 50.0, 50.0]
    coordinates[2] = [2.5, 50.0, 50.0]
    coordinates[3] = [2.5, 50.0, 50.0]
    u.load_new(coordinates, order="fac")
    dim = np.asarray([
        [100, 100, 100, 90, 90, 90],
        [95, 100, 100, 90, 90, 90],  # Box shrinks by 5 in the x-dimension
        [95, 100, 100, 90, 90, 90],
        [95, 100, 100, 90, 90, 90],
    ])
    workflow = [
        mda.transformations.boxdimensions.set_dimensions(dim),
    ]
    u.trajectory.add_transformations(*workflow)
    tmp_pdb = (tmp_path_factory.getbasetemp() / "nojump_npt_2nd_frame.pdb").as_posix()
    tmp_xtc = (tmp_path_factory.getbasetemp() / "nojump_npt_2nd_frame.xtc").as_posix()
    u.atoms.write(tmp_pdb)
    with mda.Writer(tmp_xtc) as f:
        for ts in u.trajectory:
            f.write(u.atoms)
    return mda.Universe(tmp_pdb, tmp_xtc)


def test_nojump_orthogonal_fwd(nojump_universe):
    """
    Test if the nojump transform is returning the correct
    values when iterating forwards over the sample trajectory.
    """
    u = nojump_universe
    dim = np.asarray([1, 1, 1, 90, 90, 90], np.float32)
    workflow = [mda.transformations.boxdimensions.set_dimensions(dim), NoJump()]
    u.trajectory.add_transformations(*workflow)
    transformed_coordinates = u.trajectory.timeseries()[0]
    # Step is 1 unit every 3 steps. After 99 steps from the origin,
    # we'll end up at 33.
    assert_allclose(
        transformed_coordinates, np.outer(np.linspace(0, 33, 100), np.ones(3))
    )


def test_nojump_nonorthogonal_fwd(nojump_universe):
    """
    Test if the nojump transform is returning the correct
    values when iterating forwards over the sample trajectory.
    """
    u = nojump_universe
    # Set a non-orthogonal box dimension. The box below works out to be this cell:
    # [[1.        0.        0.       ]
    # [0.        1.        0.       ]
    # [0.5       0.        0.8660254]]
    dim = np.asarray([1, 1, 1, 90, 60, 90], np.float32)
    workflow = [mda.transformations.boxdimensions.set_dimensions(dim), NoJump()]
    u.trajectory.add_transformations(*workflow)
    transformed_coordinates = u.trajectory.timeseries()[0]
    # After the transformation, you should end up in a repeating pattern, since you are
    # working in a hexagonal unit cell system. Since you jump every third timestep across
    # a periodic boundary, the shift in each axis is saved. As a consequence, the correct
    # jump every third step is just the original position + the size of the periodic cells.
    assert_allclose(
        transformed_coordinates[::3],
        np.outer(np.arange(33.5), np.array([0.5, 1, np.sqrt(3) / 2])),
    )
    assert_allclose(
        transformed_coordinates[1::3],
        np.outer(np.arange(32.5), np.array([0.5, 1, np.sqrt(3) / 2])) + 1 * np.ones(3) / 3,
        rtol=1.2e-7
    )
    assert_allclose(
        transformed_coordinates[2::3],
        np.outer(np.arange(32.5), np.array([0.5, 1, np.sqrt(3) / 2])) + 2 * np.ones(3) / 3,
        rtol=1.2e-7
    )


def test_nojump_constantvel(nojump_constantvel_universe):
    """
    Test if the nojump transform is returning the correct
    values when iterating forwards over the sample trajectory.
    """
    ref = nojump_constantvel_universe
    towrap = ref.copy() # This copy of the universe will be wrapped, then unwrapped,
    # and should be equal to ref.
    dim = np.asarray([5, 5, 5, 54, 60, 90], np.float32)
    workflow = [
        mda.transformations.boxdimensions.set_dimensions(dim),
        wrap(towrap.atoms),
        NoJump(),
    ]
    towrap.trajectory.add_transformations(*workflow)
    assert_allclose(
        towrap.trajectory.timeseries(),
        ref.trajectory.timeseries(),
        rtol=5e-07,
        atol=5e-06,
    )


def test_nojump_2nd_frame(nojump_universe_npt_2nd_frame):
    """
    Test if the nojump transform returns the correct values
    at all frames when iterating over an npt trajectory
    and an atom crosses the x-boundary at the second frame

    Wrapped coordinates are:
    coordinates = [
        [97.5, 50.0, 50.0],
        [2.5, 50.0, 50.0],
        [2.5, 50.0, 50.0],
        [2.5, 50.0, 50.0],
    ]
    where each row is a different timestep.

    Unwrapped coordinates are the same at each frame:
    unwrapped = [97.5, 50.0, 50.0]
    """
    u = nojump_universe_npt_2nd_frame
    dim = np.asarray([
        [100, 100, 100, 90, 90, 90],
        [95, 100, 100, 90, 90, 90],  # Box shrinks by 5 in the x-dimension
        [95, 100, 100, 90, 90, 90],
        [95, 100, 100, 90, 90, 90],
    ])
    workflow = [
        mda.transformations.boxdimensions.set_dimensions(dim),
        NoJump(),
    ]
    u.trajectory.add_transformations(*workflow)
    x_position = 97.5
    np.testing.assert_allclose(u.trajectory.timeseries()[:, 0, 0], x_position)


def test_nojump_3rd_frame(nojump_universe_npt_3rd_frame):
    """
    Test if the nojump transform returns the correct values
    at all frames when iterating over an npt trajectory
    and an atom crosses the x-boundary at the third frame.

    Wrapped coordinates are:
    coordinates = [
        [97.5, 50.0, 50.0],
        [97.5, 50.0, 50.0],
        [2.5, 50.0, 50.0],
        [2.5, 50.0, 50.0],
    ]
    where each row is a different timestep.

    Unwarpped coordinates are the same at each frame:
    unwrapped = [97.5, 50.0, 50.0]
    """
    u = nojump_universe_npt_3rd_frame
    dim = np.asarray([
        [100, 100, 100, 90, 90, 90],
        [100, 100, 100, 90, 90, 90],
        [95, 100, 100, 90, 90, 90],  # Box shrinks by 5 in the x-dimension
        [95, 100, 100, 90, 90, 90],
    ])
    workflow = [
        mda.transformations.boxdimensions.set_dimensions(dim),
        NoJump(),
    ]
    u.trajectory.add_transformations(*workflow)
    x_position = 97.5
    np.testing.assert_allclose(u.trajectory.timeseries()[:, 0, 0], x_position)


def test_nojump_iterate_twice(nojump_universe_npt_2nd_frame_from_file):
    """
    Test if the nojump transform always returns the correct values
    at all frames when iterating over multiple times.
    """
    u = nojump_universe_npt_2nd_frame_from_file
    u.trajectory.add_transformations(NoJump())
    timeseries_first_iteration = u.trajectory.timeseries()
    timeseries_second_iteration = u.trajectory.timeseries()
    np.testing.assert_allclose(timeseries_first_iteration, timeseries_second_iteration)


def test_nojump_constantvel_skip(nojump_universes_fromfile):
    """
    Test if the nojump transform warning is emitted.
    """
    with pytest.warns(UserWarning):
        u = nojump_universes_fromfile
        u.trajectory[0]
        u.trajectory[9] #Exercises the warning.


def test_nojump_constantvel_stride_2(nojump_universes_fromfile):
    """
    Test if the nojump transform warning is emitted.
    """
    match = "Currently jumping between frames with a step of more than 1."
    with pytest.warns(UserWarning, match=match):
        u = nojump_universes_fromfile
        for ts in u.trajectory[::2]:  # Exercises the warning.
            pass


def test_nojump_constantvel_jumparound(nojump_universes_fromfile):
    """
    Test if the nojump transform is emitting a warning.
    """
    with pytest.warns(UserWarning):
        u = nojump_universes_fromfile
        u.trajectory[0]
        u.trajectory[1]
        u.trajectory[9]


def test_missing_dimensions_init(nojump_universe):
    """
    Test if the nojump transform raises a NoDataError if there is no
    initial dimension for the periodic unit cell.
    """
    with pytest.raises(mda.exceptions.NoDataError):
        u = nojump_universe
        workflow = [NoJump()]
        u.trajectory.add_transformations(*workflow)
        transformed_coordinates = u.trajectory.timeseries()[0]


def test_missing_dimensions(nojump_universe):
    """
    Test if the nojump transform raises a NoDataError if there is no
    dimension for the periodic unit cell in a subsequent timestep.
    """
    with pytest.raises(mda.exceptions.NoDataError):
        u = nojump_universe
        u.dimensions = [73, 73, 73, 90, 90, 90]
        workflow = [NoJump()]
        u.trajectory.add_transformations(*workflow)
        transformed_coordinates = u.trajectory.timeseries()[0]


def test_notinvertible(nojump_universe):
    """
    Test if the nojump transform raises a NoDataError if the dimensions
    are invalid for the periodic unit cell.
    """
    with pytest.raises(mda.exceptions.NoDataError):
        u = nojump_universe
        dim = [1, 0, 0, 90, 90, 90]
        workflow = [mda.transformations.boxdimensions.set_dimensions(dim),NoJump()]
        u.trajectory.add_transformations(*workflow)
        transformed_coordinates = u.trajectory.timeseries()[0]
