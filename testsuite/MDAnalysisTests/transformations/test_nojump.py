import numpy as np
import pytest
from numpy.testing import assert_allclose

import MDAnalysis as mda
from MDAnalysis.transformations import NoJump
from MDAnalysisTests import datafiles as data


@pytest.fixture()
def nojump_universe():
    """
    Create the universe objects for the tests.
    """
    u = mda.Universe.empty(N, trajectory=True)
    coordinates = np.empty((100,  # number of frames
                            u.atoms.n_atoms,
                            3))
    coordinates[::3,0] = 0 * np.ones(3) / 3
    coordinates[1::3,0] = 1 * np.ones(3) / 3
    coordinates[2::3,0] = 2 * np.ones(3) / 3
    u.load_new(coordinates, order='fac')
    return u


def test_nojump_orthogonal_fwd(nojump_universe):
    """
    Test if the nojump transform is returning the correct
    values when iterating forwards over the sample trajectory.
    """
    u = nojump_universe
    dim = np.asarray([1, 1, 1, 90, 90, 90], np.float32)
    workflow = [mda.transformations.boxdimensions.set_dimensions(dim), NoJump()]
    u.add_transformations(*workflow)
    transformed_coordinates = u.trajectory.timeseries()[0]
    # Step is 1 unit every 3 steps. After 99 steps from the origin,
    # we'll end up at 33.
    assert_allclose(transformed_coordinates, np.outer(np.linspace(0,33,100),np.ones(3)))


def test_nojump_nonorthogonal_fwd(nojump_universe):
    """
    Test if the nojump transform is returning the correct
    values when iterating forwards over the sample trajectory.
    """
    u = nojump_universe
    dim = np.asarray([1, 1, 1, 90, 60, 90], np.float32)
    workflow = [mda.transformations.boxdimensions.set_dimensions(dim), NoJump()]
    u.add_transformations(*workflow)
    transformed_coordinates = u.trajectory.timeseries()[0]
    # Step is 1 unit every 3 steps. After 99 steps from the origin,
    # we'll end up at 33. However, since the unit cell is non-orthogonal,
    # we'll end up at a distorted place.
    assert_allclose(transformed_coordinates[::3], np.outer(np.arange(33.5),np.array([0.5, 1, np.sqrt(3)/2])))
