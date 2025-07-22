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


from MDAnalysis.analysis.msd import EinsteinMSD as MSD
import MDAnalysis as mda

from numpy.testing import assert_almost_equal, assert_equal
import numpy as np

from MDAnalysisTests.datafiles import (
    PSF,
    DCD,
    RANDOM_WALK,
    RANDOM_WALK_TOPO,
    LAMMPSDUMP_non_linear,
)
from MDAnalysisTests.util import block_import, import_not_available

import pytest


@pytest.fixture(scope="module")
def SELECTION():
    selection = "backbone and name CA and resid 1-10"
    return selection


@pytest.fixture(scope="module")
def u():
    return mda.Universe(PSF, DCD)


@pytest.fixture(scope="module")
def NSTEP():
    nstep = 5000
    return nstep


@pytest.fixture(scope="module")
def random_walk_u():
    # 100x100
    return mda.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)


@pytest.fixture(scope="module")
def msd(u, SELECTION):
    # non fft msd
    m = MSD(u, SELECTION, msd_type="xyz", fft=False)
    m.run()
    return m


@pytest.fixture(scope="module")
def step_traj(NSTEP):  # constant velocity
    x = np.arange(NSTEP)
    traj = np.vstack([x, x, x]).T
    traj_reshape = traj.reshape([NSTEP, 1, 3])
    u = mda.Universe.empty(1)
    u.load_new(traj_reshape)
    return u


@block_import("tidynamics")
def test_notidynamics(u, SELECTION):
    with pytest.raises(ImportError, match="tidynamics was not found"):
        u = mda.Universe(PSF, DCD)
        msd = MSD(u, SELECTION)
        msd.run()


def characteristic_poly(n, d):
    # polynomial that describes unit step traj MSD
    x = np.arange(0, n)
    y = d * x * x
    return y


class TestMSDSimple(object):

    def test_selection_works(self, msd):
        # test some basic size and shape things
        assert_equal(msd.n_particles, 10)

    def test_ag_accepted(self, u):
        ag = u.select_atoms("resid 1")
        m = MSD(ag, msd_type="xyz", fft=False)

    def test_updating_ag_rejected(self, u):
        updating_ag = u.select_atoms("around 3.5 resid 1", updating=True)
        errmsg = "UpdatingAtomGroups are not valid"
        with pytest.raises(TypeError, match=errmsg):
            m = MSD(updating_ag, msd_type="xyz", fft=False)

    @pytest.mark.parametrize("msdtype", ["foo", "bar", "yx", "zyx"])
    def test_msdtype_error(self, u, SELECTION, msdtype):
        errmsg = f"invalid msd_type: {msdtype}"
        with pytest.raises(ValueError, match=errmsg):
            m = MSD(u, SELECTION, msd_type=msdtype)

    @pytest.mark.parametrize(
        "dim, dim_factor",
        [
            ("xyz", 3),
            ("xy", 2),
            ("xz", 2),
            ("yz", 2),
            ("x", 1),
            ("y", 1),
            ("z", 1),
        ],
    )
    def test_simple_step_traj_all_dims(self, step_traj, NSTEP, dim, dim_factor):
        # testing the "simple" algorithm on constant velocity trajectory
        # should fit the polynomial y=dim_factor*x**2
        m_simple = MSD(step_traj, "all", msd_type=dim, fft=False)
        m_simple.run()
        poly = characteristic_poly(NSTEP, dim_factor)
        assert_almost_equal(m_simple.results.timeseries, poly, decimal=4)

    @pytest.mark.parametrize(
        "dim, dim_factor",
        [
            ("xyz", 3),
            ("xy", 2),
            ("xz", 2),
            ("yz", 2),
            ("x", 1),
            ("y", 1),
            ("z", 1),
        ],
    )
    def test_simple_start_stop_step_all_dims(self, step_traj, NSTEP, dim, dim_factor):
        # testing the "simple" algorithm on constant velocity trajectory
        # test start stop step is working correctly
        m_simple = MSD(step_traj, "all", msd_type=dim, fft=False)
        m_simple.run(start=10, stop=1000, step=10)
        poly = characteristic_poly(NSTEP, dim_factor)
        # polynomial must take offset start into account
        assert_almost_equal(m_simple.results.timeseries, poly[0:990:10], decimal=4)

    def test_random_walk_u_simple(self, random_walk_u):
        # regress against random_walk test data
        msd_rw = MSD(random_walk_u, "all", msd_type="xyz", fft=False)
        msd_rw.run()
        norm = np.linalg.norm(msd_rw.results.timeseries)
        val = 3932.39927487146
        assert_almost_equal(norm, val, decimal=5)


@pytest.mark.skipif(
    import_not_available("tidynamics"),
    reason="Test skipped because tidynamics not found",
)
class TestMSDFFT(object):

    @pytest.fixture(scope="class")
    def msd_fft(self, u, SELECTION):
        # fft msd
        m = MSD(u, SELECTION, msd_type="xyz", fft=True)
        m.run()
        return m

    def test_fft_vs_simple_default(self, msd, msd_fft):
        # testing on the  PSF, DCD trajectory
        timeseries_simple = msd.results.timeseries
        timeseries_fft = msd_fft.results.timeseries
        assert_almost_equal(timeseries_simple, timeseries_fft, decimal=4)

    def test_fft_vs_simple_default_per_particle(self, msd, msd_fft):
        # check fft and simple give same result per particle
        per_particle_simple = msd.results.msds_by_particle
        per_particle_fft = msd_fft.results.msds_by_particle
        assert_almost_equal(per_particle_simple, per_particle_fft, decimal=4)

    @pytest.mark.parametrize("dim", ["xyz", "xy", "xz", "yz", "x", "y", "z"])
    def test_fft_vs_simple_all_dims(self, u, SELECTION, dim):
        # check fft and simple give same result for each dimensionality
        m_simple = MSD(u, SELECTION, msd_type=dim, fft=False)
        m_simple.run()
        timeseries_simple = m_simple.results.timeseries
        m_fft = MSD(u, SELECTION, msd_type=dim, fft=True)
        m_fft.run()
        timeseries_fft = m_fft.results.timeseries
        assert_almost_equal(timeseries_simple, timeseries_fft, decimal=4)

    @pytest.mark.parametrize("dim", ["xyz", "xy", "xz", "yz", "x", "y", "z"])
    def test_fft_vs_simple_all_dims_per_particle(self, u, SELECTION, dim):
        # check fft and simple give same result for each particle in each
        # dimension
        m_simple = MSD(u, SELECTION, msd_type=dim, fft=False)
        m_simple.run()
        per_particle_simple = m_simple.results.msds_by_particle
        m_fft = MSD(u, SELECTION, msd_type=dim, fft=True)
        m_fft.run()
        per_particle_fft = m_fft.results.msds_by_particle
        assert_almost_equal(per_particle_simple, per_particle_fft, decimal=4)

    @pytest.mark.parametrize(
        "dim, dim_factor",
        [
            ("xyz", 3),
            ("xy", 2),
            ("xz", 2),
            ("yz", 2),
            ("x", 1),
            ("y", 1),
            ("z", 1),
        ],
    )
    def test_fft_step_traj_all_dims(self, step_traj, NSTEP, dim, dim_factor):
        # testing the fft algorithm on constant velocity trajectory
        # this should fit the polynomial y=dim_factor*x**2
        # fft based tests require a slight decrease in expected prescision
        # primarily due to roundoff in fft(ifft()) calls.
        # relative accuracy expected to be around ~1e-12
        m_simple = MSD(step_traj, "all", msd_type=dim, fft=True)
        m_simple.run()
        poly = characteristic_poly(NSTEP, dim_factor)
        # this was relaxed from decimal=4 for numpy=1.13 test
        assert_almost_equal(m_simple.results.timeseries, poly, decimal=3)

    @pytest.mark.parametrize(
        "dim, dim_factor",
        [
            ("xyz", 3),
            ("xy", 2),
            ("xz", 2),
            ("yz", 2),
            ("x", 1),
            ("y", 1),
            ("z", 1),
        ],
    )
    def test_fft_start_stop_step_all_dims(self, step_traj, NSTEP, dim, dim_factor):
        # testing the fft algorithm on constant velocity trajectory
        # test start stop step is working correctly
        m_simple = MSD(step_traj, "all", msd_type=dim, fft=True)
        m_simple.run(start=10, stop=1000, step=10)
        poly = characteristic_poly(NSTEP, dim_factor)
        # polynomial must take offset start into account
        assert_almost_equal(m_simple.results.timeseries, poly[0:990:10], decimal=3)

    def test_random_walk_u_fft(self, random_walk_u):
        # regress against random_walk test data
        msd_rw = MSD(random_walk_u, "all", msd_type="xyz", fft=True)
        msd_rw.run()
        norm = np.linalg.norm(msd_rw.results.timeseries)
        val = 3932.39927487146
        assert_almost_equal(norm, val, decimal=5)


def test_msd_non_linear():
    u = mda.Universe(LAMMPSDUMP_non_linear, format="LAMMPSDUMP")

    msd = MSD(u, select="all", msd_type="xyz", non_linear=True)
    msd.run()

    result_msd = msd.results.timeseries
    result_delta_t = msd.results.delta_t_values

    assert result_msd.ndim == 1
    assert result_msd.shape[0] > 0

    assert result_delta_t.ndim == 1
    assert result_delta_t.shape[0] > 0

    expected_msd = np.array(
        [
            0.00000000e00,
            7.70976963e-05,
            2.90842662e-04,
            6.55040347e-04,
            1.20610926e-03,
            2.52547250e-03,
            3.31645965e-03,
            5.38852795e-03,
            1.01941562e-02,
            1.24745603e-02,
            1.35380300e-02,
            1.57475527e-02,
            2.85165801e-02,
            3.50591021e-02,
            3.81292797e-02,
            3.96176470e-02,
            3.83551274e-02,
            5.51041371e-02,
            5.95049433e-02,
            6.07026502e-02,
            6.14434181e-02,
            6.19512436e-02,
            6.61293773e-02,
            9.46607497e-02,
            1.01300585e-01,
            9.96583811e-02,
            9.81112279e-02,
            9.72780657e-02,
            9.69221886e-02,
            1.29442431e-01,
            1.80752226e-01,
            1.86358673e-01,
            1.98140564e-01,
            2.00603000e-01,
            1.99094789e-01,
            1.97272787e-01,
            1.96156023e-01,
            2.67664446e-01,
            4.50987076e-01,
            4.02344442e-01,
            3.91458056e-01,
            4.10370922e-01,
            4.22997445e-01,
            4.26217251e-01,
            4.26484034e-01,
            4.26360794e-01,
            6.91315347e-01,
            9.94317423e-01,
            1.19622365e00,
            1.04919180e00,
            1.06437594e00,
            1.09426432e00,
            1.10194082e00,
            1.10275424e00,
            1.10383947e00,
            1.10493159e00,
        ]
    )

    expected_delta_t = np.array(
        [
            0.000e00,
            1.000e00,
            2.000e00,
            3.000e00,
            4.000e00,
            6.000e00,
            7.000e00,
            8.000e00,
            1.200e01,
            1.400e01,
            1.500e01,
            1.600e01,
            2.400e01,
            2.800e01,
            3.000e01,
            3.100e01,
            3.200e01,
            4.800e01,
            5.600e01,
            6.000e01,
            6.200e01,
            6.300e01,
            6.400e01,
            9.600e01,
            1.120e02,
            1.200e02,
            1.240e02,
            1.260e02,
            1.270e02,
            1.280e02,
            1.920e02,
            2.240e02,
            2.400e02,
            2.480e02,
            2.520e02,
            2.540e02,
            2.550e02,
            2.560e02,
            3.840e02,
            4.480e02,
            4.800e02,
            4.960e02,
            5.040e02,
            5.080e02,
            5.100e02,
            5.110e02,
            5.120e02,
            7.680e02,
            8.960e02,
            9.600e02,
            9.920e02,
            1.008e03,
            1.016e03,
            1.020e03,
            1.022e03,
            1.023e03,
        ]
    )

    np.testing.assert_allclose(result_msd, expected_msd, rtol=1e-5)
    np.testing.assert_allclose(result_delta_t, expected_delta_t, rtol=1e-5)
