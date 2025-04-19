from MDAnalysisTests.datafiles import (
    TPR2024_4_bonded,
    TPR_EXTRA_2024_4,
    TPR2024_4,
    TPR2024,
    TPR2023,
    TPR_xvf_2024_4,
    TPR_NNPOT_2025_0,
)
import MDAnalysis as mda


import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal


@pytest.mark.parametrize(
    "tpr_file, exp_first_atom, exp_last_atom, exp_shape, exp_vel_first_atom, exp_vel_last_atom",
    [
        # this case is an alanine dipeptide
        # with neural network potential active
        # and nonzero velocities
        (
            TPR_NNPOT_2025_0,  # tpx 137
            [2.36700e00, 2.30000e-02, 9.20000e-02],
            [2.95100e00, 2.00000e-01, 2.41000e-01],
            (23, 3),
            [-4.72100e-01, -2.20900e-01, -2.42800e-01],
            [-1.11900e-01, -3.69300e-01, -6.10000e-03],
        ),
        (
            TPR2024_4_bonded,  # tpx 134
            [4.446, 4.659, 2.384],
            [4.446, 4.659, 2.384],
            (14, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        # same coordinates, different shape
        (
            TPR_EXTRA_2024_4,  # tpx 134
            [4.446, 4.659, 2.384],
            [4.446, 4.659, 2.384],
            (18, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        # different coordinates and different shape
        (
            TPR2024_4,  # tpx 134
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        # nonzero velocities
        (
            TPR_xvf_2024_4,  # tpx 134
            [3.19900e00, 1.62970e00, 1.54480e00],
            [3.39350e00, 3.49420e00, 3.02400e00],
            (19385, 3),
            [-0.20668714, 0.26678202, -0.10564042],
            [-3.38010e-02, -3.22064e-01, -1.9863836e-01],
        ),
        (
            TPR2024,  # tpx 133
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2023,  # tpx 129
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
    ],
)
def test_basic_read_tpr(
    tpr_file,
    exp_first_atom,
    exp_last_atom,
    exp_shape,
    exp_vel_first_atom,
    exp_vel_last_atom,
):
    # verify basic ability to read positions and
    # velocities from GMX .tpr files
    # expected values are from gmx dump
    u = mda.Universe(tpr_file)
    assert_allclose(u.atoms.positions[0, ...], exp_first_atom)
    assert_allclose(u.atoms.positions[-1, ...], exp_last_atom)
    assert_equal(u.atoms.positions.shape, exp_shape)
    assert_allclose(u.atoms.velocities[0, ...], exp_vel_first_atom)
    assert_allclose(u.atoms.velocities[-1, ...], exp_vel_last_atom)
    assert_equal(u.atoms.velocities.shape, exp_shape)
