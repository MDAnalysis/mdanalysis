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
    towrap = mda.Universe.empty(Natom, trajectory=True)
    towrap.load_new(coordinates, order="fac")
    return reference, towrap


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
    dim = np.asarray([1, 1, 1, 90, 60, 90], np.float32)
    workflow = [mda.transformations.boxdimensions.set_dimensions(dim), NoJump()]
    u.trajectory.add_transformations(*workflow)
    transformed_coordinates = u.trajectory.timeseries()[0]
    # Step is 1 unit every 3 steps. After 99 steps from the origin,
    # we'll end up at 33. However, since the unit cell is non-orthogonal,
    # we'll end up at a distorted place.
    assert_allclose(
        transformed_coordinates[::3],
        np.outer(np.arange(33.5), np.array([0.5, 1, np.sqrt(3) / 2])),
    )


def test_nojump_constantvel(nojump_constantvel_universe):
    """
    Test if the nojump transform is returning the correct
    values when iterating forwards over the sample trajectory.
    """
    ref, towrap = nojump_constantvel_universe
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


def test_nojump_constantvel_skip(nojump_universes_fromfile):
    """
    Test if the nojump transform warning is emitted.
    """
    with pytest.warns(UserWarning):
        u = nojump_universes_fromfile
        u.trajectory[0]
        u.trajectory[9] #Exercises the warning.


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
    with pytest.raises(mda.exceptions.NoDataError):
        u = nojump_universe
        workflow = [NoJump()]
        u.trajectory.add_transformations(*workflow)
        transformed_coordinates = u.trajectory.timeseries()[0]


def test_missing_dimensions(nojump_universe):
    with pytest.raises(mda.exceptions.NoDataError):
        u = nojump_universe
        u.dimensions = [73, 73, 73, 90, 90, 90]
        workflow = [NoJump()]
        u.trajectory.add_transformations(*workflow)
        transformed_coordinates = u.trajectory.timeseries()[0]
