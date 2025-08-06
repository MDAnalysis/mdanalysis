# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from MDAnalysisTests.datafiles import (
    TPR2024_4_bonded,
    TPR_EXTRA_2024_4,
    TPR2024_4,
    TPR2024,
    TPR2023,
    TPR2022RC1,
    TPR2022RC1_bonded,
    TPR2021,
    TPR2021Double,
    TPR2021_bonded,
    TPR2020,
    TPR2020Double,
    TPR2020_bonded,
    TPR2019B3,
    TPR2019B3_bonded,
    TPR2018,
    TPR2018_bonded,
    TPR2016,
    TPR2016_bonded,
    TPR510,
    TPR510_bonded,
    TPR502,
    TPR504,
    TPR505,
    TPR460,
    TPR461,
    TPR455,
    TPR454,
    TPR453,
    TPR452,
    TPR451,
    TPR450,
    TPR400,
    TPR402,
    TPR403,
    TPR404,
    TPR405,
    TPR406,
    TPR407,
    TPR455Double,
    TPR_xvf_2024_4,
    TPR_NNPOT_2025_0,
    TPR2020B2,
    INPCRD,
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
        (
            TPR2022RC1,  # tpx 127
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2022RC1_bonded,  # tpx 127
            [4.446, 4.659, 2.384],
            [4.446, 4.659, 2.384],
            (14, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2021,  # tpx 122
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2021Double,  # tpx 122
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2021_bonded,  # tpx 122
            [4.446, 4.659, 2.384],
            [4.446, 4.659, 2.384],
            (14, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2020,  # tpx 119
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2020Double,  # tpx 119
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2020_bonded,  # tpx 119
            [4.446, 4.659, 2.384],
            [4.446, 4.659, 2.384],
            (14, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2019B3,  # tpx 116
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2019B3_bonded,  # tpx 116
            [4.446, 4.659, 2.384],
            [4.446, 4.659, 2.384],
            (14, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2018,  # tpx 112
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2018_bonded,  # tpx 112
            [4.446, 4.659, 2.384],
            [4.446, 4.659, 2.384],
            (14, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2016,  # tpx 110
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR2016_bonded,  # tpx 110
            [4.446, 4.659, 2.384],
            [4.446, 4.659, 2.384],
            (14, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR510,  # tpx 103
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR510_bonded,  # tpx 103
            [4.446, 4.659, 2.384],
            [4.446, 4.659, 2.384],
            (14, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR502,  # tpx 100
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR504,  # tpx 100
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR505,  # tpx 100
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR460,  # tpx 83
            [1.68200e00, 1.59800e00, 3.32200e00],
            [6.11800e00, 3.62100e00, 5.94800e00],
            (44_052, 3),
            [-0.07961109, -0.21509364, -0.42340994],
            [0.30528614, 0.42575037, 0.3387295],
        ),
        # same tpx version, same positions,
        # but different velocities confirmed by
        # gmx dump
        (
            TPR461,  # tpx 83
            [1.68200e00, 1.59800e00, 3.32200e00],
            [6.11800e00, 3.62100e00, 5.94800e00],
            (44_052, 3),
            [0.90468854, -0.33816168, -0.21799089],
            [-0.22830929, -0.15923975, 0.10457583],
        ),
        (
            TPR455,  # tpx 73
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR454,  # tpx 73
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR453,  # tpx 73
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR452,  # tpx 73
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR451,  # tpx 73
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR450,  # tpx 73
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR455Double,  # tpx 73
            [1.48700e00, 4.02600e00, 7.95400e00],
            [2.90100e00, 2.67900e00, 2.23600e00],
            (21692, 3),
            [0.10503311, 0.22275597, 0.002303541],
            [0.089259855, 1.0131412, 2.0914786],
        ),
        (
            TPR400,  # tpx 58
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR402,  # tpx 58
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR403,  # tpx 58
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR404,  # tpx 58
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR405,  # tpx 58
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR406,  # tpx 58
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
        (
            TPR407,  # tpx 58
            [3.25000e-01, 1.00400e00, 1.03800e00],
            [-2.56000e-01, 1.37300e00, 3.59800e00],
            (2263, 3),
            np.zeros(3),
            np.zeros(3),
        ),
    ],
)
@pytest.mark.parametrize("double_incantation", [True, False])
def test_basic_read_tpr(
    tpr_file,
    exp_first_atom,
    exp_last_atom,
    exp_shape,
    exp_vel_first_atom,
    exp_vel_last_atom,
    double_incantation,
):
    # verify basic ability to read positions and
    # velocities from GMX .tpr files
    # expected values are from gmx dump
    if double_incantation:
        u = mda.Universe(tpr_file, tpr_file)
    else:
        u = mda.Universe(tpr_file)
    assert_allclose(u.atoms.positions[0, ...], exp_first_atom)
    assert_allclose(u.atoms.positions[-1, ...], exp_last_atom)
    assert_equal(u.atoms.positions.shape, exp_shape)
    assert_allclose(u.atoms.velocities[0, ...], exp_vel_first_atom)
    assert_allclose(u.atoms.velocities[-1, ...], exp_vel_last_atom)
    assert_equal(u.atoms.velocities.shape, exp_shape)


def test_error_handling_incorrect_format():
    with pytest.raises(IOError, match="Invalid tpr coordinate file"):
        u = mda.Universe(TPR2024_4_bonded, INPCRD, format="TPR")


def test_2020B2_unsupported():
    with pytest.raises(IOError, match="beta versions of gromacs 2020"):
        u = mda.Universe(TPR2020, TPR2020B2)


def test_different_versions():
    # for now, we only have crude testing for
    # reading topology and positions/velocities
    # for the same system with different TPR
    # file versions
    exp_first_atom = [3.25000e-01, 1.00400e00, 1.03800e00]
    exp_last_atom = [-2.56000e-01, 1.37300e00, 3.59800e00]
    exp_shape = (2263, 3)
    exp_vel = np.zeros(3)
    u = mda.Universe(TPR2020, TPR2024_4)
    assert_allclose(u.atoms.positions[0, ...], exp_first_atom)
    assert_allclose(u.atoms.positions[-1, ...], exp_last_atom)
    assert_allclose(u.atoms.velocities[0, ...], exp_vel)
    assert_allclose(u.atoms.velocities[-1, ...], exp_vel)
    assert_equal(u.atoms.positions.shape, exp_shape)
    assert_equal(u.atoms.velocities.shape, exp_shape)
