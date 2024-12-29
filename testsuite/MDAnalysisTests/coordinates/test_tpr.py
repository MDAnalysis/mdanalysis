from MDAnalysisTests.datafiles import (TPR2024_4_bonded,
                                       TPR_EXTRA_2024_4,
                                       TPR2024_4)
import MDAnalysis as mda


import pytest
from numpy.testing import assert_allclose, assert_equal


@pytest.mark.parametrize("tpr_file, exp_first_atom, exp_last_atom, exp_shape", [
    (TPR2024_4_bonded,
     [4.446, 4.659, 2.384],
     [4.446, 4.659, 2.384],
     (14, 3),
    ),
    # same coordinates, different shape
    (TPR_EXTRA_2024_4,
     [4.446, 4.659, 2.384],
     [4.446, 4.659, 2.384],
     (18, 3),
    ),
    # different coordinates and different shape
    (TPR2024_4,
     [3.25000e-01,  1.00400e+00,  1.03800e+00],
     [-2.56000e-01,  1.37300e+00,  3.59800e+00],
     (2263, 3),
    ),
])
def test_basic_read_tpr(tpr_file, exp_first_atom, exp_last_atom, exp_shape):
    u = mda.Universe(tpr_file)
    assert_allclose(u.atoms.positions[0, ...], exp_first_atom)
    assert_allclose(u.atoms.positions[-1, ...], exp_last_atom)
    assert_equal(u.atoms.positions.shape, exp_shape)
