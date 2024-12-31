from MDAnalysisTests.datafiles import (TPR2024_4_bonded,
                                       TPR_EXTRA_2024_4,
                                       TPR2024_4,
                                       TPR2024,
                                       TPR2023,
                                       TPR_xvf_2024_4)
import MDAnalysis as mda


import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal


@pytest.mark.parametrize("tpr_file, exp_first_atom, exp_last_atom, exp_shape, exp_vel_first_atom, exp_vel_last_atom", [
    (TPR2024_4_bonded, # tpx 134
     [4.446, 4.659, 2.384],
     [4.446, 4.659, 2.384],
     (14, 3),
     np.zeros(3),
     np.zeros(3),
    ),
    # same coordinates, different shape
    (TPR_EXTRA_2024_4, # tpx 134
     [4.446, 4.659, 2.384],
     [4.446, 4.659, 2.384],
     (18, 3),
     np.zeros(3),
     np.zeros(3),
    ),
    # different coordinates and different shape
    (TPR2024_4, # tpx 134
     [3.25000e-01,  1.00400e+00,  1.03800e+00],
     [-2.56000e-01,  1.37300e+00,  3.59800e+00],
     (2263, 3),
     np.zeros(3),
     np.zeros(3),
    ),
    # nonzero velocities
    (TPR_xvf_2024_4, # tpx 134
     [3.19900e+00,  1.62970e+00,  1.54480e+00],
     [3.39350e+00,  3.49420e+00,  3.02400e+00],
     (19385, 3),
     [-2.06687e-01,  2.66782e-01, -1.05640e-01],
     [-3.38010e-02, -3.22064e-01, -1.98638e-01],
    ),
    (TPR2024, # tpx 133
     [3.25000e-01,  1.00400e+00,  1.03800e+00],
     [-2.56000e-01,  1.37300e+00,  3.59800e+00],
     (2263, 3),
     np.zeros(3),
     np.zeros(3),
    ),
    (TPR2023, # tpx 129
     [3.25000e-01,  1.00400e+00,  1.03800e+00],
     [-2.56000e-01,  1.37300e+00,  3.59800e+00],
     (2263, 3),
     np.zeros(3),
     np.zeros(3),
    ),
])
def test_basic_read_tpr(tpr_file,
                        exp_first_atom,
                        exp_last_atom,
                        exp_shape,
                        exp_vel_first_atom,
                        exp_vel_last_atom):
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
