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

from numpy.testing import assert_almost_equal, assert_equal, assert_allclose
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
def u_nonlinear():
    return mda.Universe(LAMMPSDUMP_non_linear, format="LAMMPSDUMP")


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
    def test_simple_step_traj_all_dims(
        self, step_traj, NSTEP, dim, dim_factor
    ):
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
    def test_simple_start_stop_step_all_dims(
        self, step_traj, NSTEP, dim, dim_factor
    ):
        # testing the "simple" algorithm on constant velocity trajectory
        # test start stop step is working correctly
        m_simple = MSD(step_traj, "all", msd_type=dim, fft=False)
        m_simple.run(start=10, stop=1000, step=10)
        poly = characteristic_poly(NSTEP, dim_factor)
        # polynomial must take offset start into account
        assert_almost_equal(
            m_simple.results.timeseries, poly[0:990:10], decimal=4
        )

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
    def test_fft_start_stop_step_all_dims(
        self, step_traj, NSTEP, dim, dim_factor
    ):
        # testing the fft algorithm on constant velocity trajectory
        # test start stop step is working correctly
        m_simple = MSD(step_traj, "all", msd_type=dim, fft=True)
        m_simple.run(start=10, stop=1000, step=10)
        poly = characteristic_poly(NSTEP, dim_factor)
        # polynomial must take offset start into account
        assert_almost_equal(
            m_simple.results.timeseries, poly[0:990:10], decimal=3
        )

    def test_random_walk_u_fft(self, random_walk_u):
        # regress against random_walk test data
        msd_rw = MSD(random_walk_u, "all", msd_type="xyz", fft=True)
        msd_rw.run()
        norm = np.linalg.norm(msd_rw.results.timeseries)
        val = 3932.39927487146
        assert_almost_equal(norm, val, decimal=5)


class TestMSDNonLinear:

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
    def test_all_msd_types(self, u_nonlinear, dim, dim_factor):
        msd = MSD(u_nonlinear, select="all", msd_type=dim, non_linear=True)
        msd.run()
        result_msd = msd.results.timeseries
        result_delta_t = msd.results.delta_t_values
        result_msd_per_particle = msd.results.msds_by_particle
        expected_results_msd = {
            "xyz": np.array(
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
            ),
            "xy": np.array(
                [
                    0.00000000e00,
                    4.71534353e-05,
                    2.00284753e-04,
                    4.34113725e-04,
                    9.36004413e-04,
                    1.92444776e-03,
                    2.48128885e-03,
                    3.74869519e-03,
                    7.45690928e-03,
                    9.32883150e-03,
                    1.02001429e-02,
                    1.18451152e-02,
                    2.16722754e-02,
                    2.69028318e-02,
                    2.92942088e-02,
                    3.04040224e-02,
                    3.09111456e-02,
                    4.32505110e-02,
                    4.71624323e-02,
                    4.84247750e-02,
                    4.88898846e-02,
                    4.91026665e-02,
                    5.18382213e-02,
                    7.64373669e-02,
                    8.16001408e-02,
                    7.79840283e-02,
                    7.61576850e-02,
                    7.51756311e-02,
                    7.46696051e-02,
                    6.91138291e-02,
                    1.32156093e-01,
                    1.32978164e-01,
                    1.52259070e-01,
                    1.47420312e-01,
                    1.46571506e-01,
                    1.46636604e-01,
                    1.46716952e-01,
                    1.49788156e-01,
                    2.34954789e-01,
                    2.15475722e-01,
                    1.99697304e-01,
                    2.35851828e-01,
                    2.39591438e-01,
                    2.41774318e-01,
                    2.43132227e-01,
                    2.43853085e-01,
                    4.73633489e-01,
                    8.32733529e-01,
                    1.01977895e00,
                    8.72773844e-01,
                    8.63705171e-01,
                    9.10979795e-01,
                    9.22629515e-01,
                    9.24681213e-01,
                    9.26321374e-01,
                    9.27520965e-01,
                ]
            ),
            "xz": np.array(
                [
                    0.00000000e00,
                    4.16965032e-05,
                    1.32580865e-04,
                    3.15741675e-04,
                    5.00108758e-04,
                    1.02934249e-03,
                    1.37025814e-03,
                    3.14244287e-03,
                    5.31340742e-03,
                    6.09128877e-03,
                    6.42306131e-03,
                    1.20010188e-02,
                    1.98178137e-02,
                    2.30804487e-02,
                    2.44381472e-02,
                    2.50944884e-02,
                    2.58745154e-02,
                    3.91995598e-02,
                    4.33821430e-02,
                    4.39535754e-02,
                    4.42724505e-02,
                    4.45635874e-02,
                    3.41294299e-02,
                    4.98672274e-02,
                    5.84413430e-02,
                    6.25792031e-02,
                    6.41470811e-02,
                    6.45891296e-02,
                    6.47562187e-02,
                    1.09986645e-01,
                    1.43454970e-01,
                    1.33856061e-01,
                    1.39832186e-01,
                    1.43856499e-01,
                    1.43524801e-01,
                    1.42333307e-01,
                    1.41539074e-01,
                    1.92375956e-01,
                    3.44024779e-01,
                    3.28460560e-01,
                    3.16161165e-01,
                    3.18352125e-01,
                    3.24388442e-01,
                    3.26201142e-01,
                    3.26221383e-01,
                    3.26069834e-01,
                    3.57054087e-01,
                    4.65715415e-01,
                    5.33167435e-01,
                    5.24342584e-01,
                    5.54527169e-01,
                    5.31925547e-01,
                    5.20297822e-01,
                    5.15474719e-01,
                    5.14206075e-01,
                    5.14110889e-01,
                ]
            ),
            "yz": np.array(
                [
                    0.00000000e00,
                    6.53454542e-05,
                    2.48819706e-04,
                    5.60225294e-04,
                    9.76105353e-04,
                    2.09715474e-03,
                    2.78137231e-03,
                    3.88591783e-03,
                    7.61799579e-03,
                    9.52900038e-03,
                    1.04528558e-02,
                    7.64897138e-03,
                    1.55430710e-02,
                    2.01349237e-02,
                    2.25262033e-02,
                    2.37367832e-02,
                    1.99245940e-02,
                    2.77582034e-02,
                    2.84653112e-02,
                    2.90269501e-02,
                    2.97245011e-02,
                    3.02362334e-02,
                    4.62911034e-02,
                    6.30169052e-02,
                    6.25596862e-02,
                    5.87535308e-02,
                    5.59176898e-02,
                    5.47913707e-02,
                    5.44185534e-02,
                    7.97843873e-02,
                    8.58933884e-02,
                    1.05883121e-01,
                    1.04189871e-01,
                    1.09929188e-01,
                    1.08093271e-01,
                    1.05575663e-01,
                    1.04056021e-01,
                    1.93164780e-01,
                    3.22994584e-01,
                    2.60752603e-01,
                    2.67057643e-01,
                    2.66537890e-01,
                    2.82015010e-01,
                    2.84459043e-01,
                    2.83614457e-01,
                    2.82798669e-01,
                    5.51943117e-01,
                    6.90185902e-01,
                    8.39500901e-01,
                    7.01267168e-01,
                    7.10519535e-01,
                    7.45623290e-01,
                    7.60954298e-01,
                    7.65352545e-01,
                    7.67151495e-01,
                    7.68231328e-01,
                ]
            ),
            "x": np.array(
                [
                    0.00000000e00,
                    1.17522422e-05,
                    4.20229558e-05,
                    9.48150526e-05,
                    2.30003909e-04,
                    4.28317752e-04,
                    5.35087342e-04,
                    1.50261012e-03,
                    2.57616045e-03,
                    2.94555994e-03,
                    3.08517422e-03,
                    8.09858129e-03,
                    1.29735090e-02,
                    1.49241784e-02,
                    1.56030763e-02,
                    1.58808638e-02,
                    1.84305335e-02,
                    2.73459337e-02,
                    3.10396321e-02,
                    3.16757002e-02,
                    3.17189170e-02,
                    3.17150102e-02,
                    1.98382739e-02,
                    3.16438446e-02,
                    3.87408988e-02,
                    4.09048503e-02,
                    4.21935381e-02,
                    4.24866950e-02,
                    4.25036352e-02,
                    4.96580433e-02,
                    9.48588375e-02,
                    8.04755522e-02,
                    9.39506926e-02,
                    9.06738115e-02,
                    9.10015182e-02,
                    9.16971239e-02,
                    9.21000022e-02,
                    7.44996658e-02,
                    1.27992492e-01,
                    1.41591839e-01,
                    1.24400413e-01,
                    1.43833032e-01,
                    1.40982435e-01,
                    1.41758209e-01,
                    1.42869576e-01,
                    1.43562125e-01,
                    1.39372230e-01,
                    3.04131521e-01,
                    3.56722744e-01,
                    3.47924630e-01,
                    3.53856403e-01,
                    3.48641025e-01,
                    3.40986519e-01,
                    3.37401693e-01,
                    3.36687977e-01,
                    3.36700263e-01,
                ]
            ),
            "y": np.array(
                [
                    0.00000000e00,
                    3.54011931e-05,
                    1.58261797e-04,
                    3.39298672e-04,
                    7.06000504e-04,
                    1.49613000e-03,
                    1.94620151e-03,
                    2.24608508e-03,
                    4.88074882e-03,
                    6.38327155e-03,
                    7.11496870e-03,
                    3.74653386e-03,
                    8.69876636e-03,
                    1.19786534e-02,
                    1.36911325e-02,
                    1.45231586e-02,
                    1.24806121e-02,
                    1.59045773e-02,
                    1.61228002e-02,
                    1.67490749e-02,
                    1.71709676e-02,
                    1.73876563e-02,
                    3.19999474e-02,
                    4.47935223e-02,
                    4.28592420e-02,
                    3.70791780e-02,
                    3.39641469e-02,
                    3.26889361e-02,
                    3.21659699e-02,
                    1.94557858e-02,
                    3.72972556e-02,
                    5.25026117e-02,
                    5.83083776e-02,
                    5.67465008e-02,
                    5.55699879e-02,
                    5.49394802e-02,
                    5.46169493e-02,
                    7.52884900e-02,
                    1.06962297e-01,
                    7.38838824e-02,
                    7.52968914e-02,
                    9.20187962e-02,
                    9.86090028e-02,
                    1.00016109e-01,
                    1.00262650e-01,
                    1.00290960e-01,
                    3.34261260e-01,
                    5.28602008e-01,
                    6.63056210e-01,
                    5.24849214e-01,
                    5.09848769e-01,
                    5.62338769e-01,
                    5.81642995e-01,
                    5.87279520e-01,
                    5.89633397e-01,
                    5.90820702e-01,
                ]
            ),
            "z": np.array(
                [
                    0.00000000e00,
                    2.99442610e-05,
                    9.05579089e-05,
                    2.20926622e-04,
                    2.70104849e-04,
                    6.01024741e-04,
                    8.35170802e-04,
                    1.63983276e-03,
                    2.73724696e-03,
                    3.14572883e-03,
                    3.33788709e-03,
                    3.90243752e-03,
                    6.84430468e-03,
                    8.15627028e-03,
                    8.83507087e-03,
                    9.21362461e-03,
                    7.44398187e-03,
                    1.18536261e-02,
                    1.23425110e-02,
                    1.22778752e-02,
                    1.25535335e-02,
                    1.28485771e-02,
                    1.42911560e-02,
                    1.82233829e-02,
                    1.97004442e-02,
                    2.16743528e-02,
                    2.19535430e-02,
                    2.21024346e-02,
                    2.22525835e-02,
                    6.03286015e-02,
                    4.85961328e-02,
                    5.33805088e-02,
                    4.58814937e-02,
                    5.31826873e-02,
                    5.25232828e-02,
                    5.06361832e-02,
                    4.94390715e-02,
                    1.17876290e-01,
                    2.16032287e-01,
                    1.86868721e-01,
                    1.91760752e-01,
                    1.74519093e-01,
                    1.83406007e-01,
                    1.84442934e-01,
                    1.83351807e-01,
                    1.82507709e-01,
                    2.17681858e-01,
                    1.61583894e-01,
                    1.76444690e-01,
                    1.76417954e-01,
                    2.00670767e-01,
                    1.83284521e-01,
                    1.79311302e-01,
                    1.78073026e-01,
                    1.77518098e-01,
                    1.77410626e-01,
                ]
            ),
        }
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
        expected_msd = expected_results_msd[dim]
        assert result_msd.shape == expected_msd.shape
        assert result_delta_t.shape == expected_delta_t.shape
        assert_allclose(result_msd, expected_msd, rtol=1e-5)
        assert_allclose(result_delta_t, expected_delta_t, rtol=1e-5)

    def test_msds_per_particle(self, u_nonlinear):
        msd = MSD(u_nonlinear, select="all", msd_type="xyz", non_linear=True)
        msd.run()
        result_msd_per_particle = msd.results.msds_by_particle
        expected_msd_per_particle = np.array(
            [
                [
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                ],
                [
                    3.29547038e-05,
                    5.48541720e-06,
                    1.82997159e-04,
                    5.86525130e-06,
                    2.00792165e-05,
                    2.19020443e-04,
                    1.07918182e-04,
                    1.40443055e-04,
                    1.31939923e-05,
                    6.32751799e-05,
                    9.77258060e-05,
                    9.01585336e-05,
                    1.60915810e-04,
                    5.26360600e-06,
                    1.11690897e-05,
                ],
                [
                    1.11020650e-04,
                    2.95432983e-05,
                    8.14488783e-04,
                    5.48629032e-05,
                    6.14942110e-05,
                    7.83963249e-04,
                    3.56422974e-04,
                    2.87535921e-04,
                    4.34042777e-05,
                    2.12105602e-04,
                    4.61613050e-04,
                    5.22839753e-04,
                    5.17809243e-04,
                    2.41769339e-05,
                    8.13590741e-05,
                ],
                [
                    2.63813577e-04,
                    5.56142040e-05,
                    1.76485543e-03,
                    8.62745947e-05,
                    1.34402031e-04,
                    1.83009465e-03,
                    8.41153775e-04,
                    8.22875827e-04,
                    1.03974497e-04,
                    4.98132900e-04,
                    9.78943150e-04,
                    1.04020679e-03,
                    1.25445600e-03,
                    4.90201928e-05,
                    1.01787577e-04,
                ],
                [
                    1.95492961e-04,
                    3.94420946e-04,
                    3.48568599e-03,
                    9.03695662e-04,
                    9.07771875e-04,
                    2.06086619e-03,
                    1.42940438e-03,
                    3.74201296e-04,
                    1.60095136e-04,
                    9.84441944e-04,
                    2.28692875e-03,
                    2.65435437e-03,
                    7.15362199e-04,
                    1.98882743e-04,
                    1.34003449e-03,
                ],
                [
                    5.98546074e-04,
                    6.27848276e-04,
                    7.62002755e-03,
                    1.39002870e-03,
                    1.26137879e-03,
                    5.28286425e-03,
                    3.03226555e-03,
                    9.07422454e-04,
                    3.49758022e-04,
                    1.89043468e-03,
                    4.75146946e-03,
                    5.44308133e-03,
                    2.35025800e-03,
                    3.44575110e-04,
                    2.03212918e-03,
                ],
                [
                    9.06110930e-04,
                    7.07433705e-04,
                    1.01047072e-02,
                    1.49807159e-03,
                    1.31984901e-03,
                    7.59302815e-03,
                    4.07066270e-03,
                    1.53285325e-03,
                    4.83971744e-04,
                    2.45190217e-03,
                    6.14914359e-03,
                    6.84060395e-03,
                    3.68796213e-03,
                    4.04686522e-04,
                    1.99590813e-03,
                ],
                [
                    9.90922374e-04,
                    4.06796762e-03,
                    9.20377615e-03,
                    7.26651738e-03,
                    8.90031873e-03,
                    7.36212283e-03,
                    4.05492755e-03,
                    6.15373092e-03,
                    2.51891323e-03,
                    1.05901041e-02,
                    9.11065105e-03,
                    3.45879569e-03,
                    3.63240972e-03,
                    1.60079393e-03,
                    1.91596798e-03,
                ],
                [
                    3.81979640e-04,
                    6.55204030e-03,
                    2.38731745e-02,
                    1.25663249e-02,
                    1.54929261e-02,
                    1.37795487e-02,
                    9.91460410e-03,
                    7.60385905e-03,
                    3.77585457e-03,
                    1.72474579e-02,
                    1.95383917e-02,
                    1.09266188e-02,
                    2.78197088e-03,
                    2.57883530e-03,
                    5.89875714e-03,
                ],
                [
                    2.16478569e-04,
                    7.45433617e-03,
                    3.32705393e-02,
                    1.42417599e-02,
                    1.67581842e-02,
                    1.84839473e-02,
                    1.37015482e-02,
                    6.23839337e-03,
                    4.31364206e-03,
                    1.92358004e-02,
                    2.54031426e-02,
                    1.53514380e-02,
                    2.24383608e-03,
                    2.89311686e-03,
                    7.31224184e-03,
                ],
                [
                    2.42450827e-04,
                    7.81072883e-03,
                    3.81783636e-02,
                    1.47116919e-02,
                    1.67693081e-02,
                    2.14513914e-02,
                    1.58174735e-02,
                    5.32206634e-03,
                    4.59639255e-03,
                    1.99454904e-02,
                    2.82963822e-02,
                    1.72861295e-02,
                    2.30458753e-03,
                    2.99275187e-03,
                    7.34524157e-03,
                ],
                [
                    9.57952908e-03,
                    2.40431962e-02,
                    4.26129746e-03,
                    5.37208999e-02,
                    1.25439480e-02,
                    1.26534763e-02,
                    8.24814401e-03,
                    1.13230807e-02,
                    2.28689772e-02,
                    1.16296926e-02,
                    9.84728085e-03,
                    3.53700031e-02,
                    1.83373150e-03,
                    4.66237709e-03,
                    1.36276561e-02,
                ],
                [
                    1.42789247e-02,
                    4.49616170e-02,
                    2.28301928e-02,
                    9.79386640e-02,
                    3.74167756e-02,
                    4.56202170e-03,
                    1.36675322e-02,
                    2.59318735e-02,
                    3.33500859e-02,
                    4.25564267e-02,
                    3.52603202e-02,
                    1.88316717e-02,
                    1.03127103e-02,
                    1.16787130e-02,
                    1.41711717e-02,
                ],
                [
                    1.13383621e-02,
                    5.28644852e-02,
                    4.34575366e-02,
                    1.12140195e-01,
                    4.91246783e-02,
                    1.13449105e-02,
                    1.81233278e-02,
                    3.14404287e-02,
                    3.75384545e-02,
                    5.34189770e-02,
                    5.37724765e-02,
                    1.67492262e-02,
                    8.32656867e-03,
                    1.37936572e-02,
                    1.24532467e-02,
                ],
                [
                    9.28934888e-03,
                    5.52802469e-02,
                    5.60822376e-02,
                    1.16430995e-01,
                    5.20844630e-02,
                    1.80384823e-02,
                    2.28397570e-02,
                    3.08590754e-02,
                    3.98280457e-02,
                    5.52886446e-02,
                    6.35076370e-02,
                    1.87768302e-02,
                    6.26772415e-03,
                    1.43155756e-02,
                    1.30501318e-02,
                ],
                [
                    8.24139445e-03,
                    5.60812412e-02,
                    6.26075889e-02,
                    1.17942565e-01,
                    5.26841868e-02,
                    2.22288270e-02,
                    2.59414039e-02,
                    3.00244323e-02,
                    4.10809178e-02,
                    5.54656041e-02,
                    6.82505538e-02,
                    2.01289939e-02,
                    5.43876016e-03,
                    1.44529277e-02,
                    1.36953083e-02,
                ],
                [
                    7.26572104e-03,
                    3.78830548e-04,
                    2.66685334e-03,
                    1.44832527e-03,
                    8.15123480e-04,
                    8.44108736e-02,
                    9.38697999e-02,
                    2.46033816e-02,
                    1.24227774e-01,
                    4.34564373e-02,
                    7.66600866e-03,
                    6.66236940e-02,
                    3.82861185e-02,
                    5.00133501e-02,
                    2.95946204e-02,
                ],
                [
                    3.04077466e-02,
                    2.90646969e-02,
                    1.34089199e-02,
                    4.13755194e-02,
                    7.07925708e-03,
                    1.18413483e-01,
                    1.06534420e-01,
                    4.38493195e-02,
                    8.91309851e-02,
                    6.21699068e-02,
                    1.35562900e-02,
                    1.39625085e-01,
                    4.84668860e-02,
                    4.57241120e-02,
                    3.77554297e-02,
                ],
                [
                    3.63440618e-02,
                    5.20567581e-02,
                    3.95215313e-02,
                    8.01595816e-02,
                    2.73821643e-02,
                    7.13045214e-02,
                    7.52595929e-02,
                    4.34848262e-02,
                    1.07567023e-01,
                    8.76068082e-02,
                    3.45693293e-02,
                    1.01547952e-01,
                    6.11471364e-02,
                    4.83617131e-02,
                    2.62611499e-02,
                ],
                [
                    3.23117802e-02,
                    6.07097117e-02,
                    6.51664345e-02,
                    9.29985226e-02,
                    3.76280686e-02,
                    5.77555913e-02,
                    5.97845004e-02,
                    5.12750264e-02,
                    1.12746984e-01,
                    9.16156305e-02,
                    4.77102418e-02,
                    7.31038407e-02,
                    5.89875917e-02,
                    5.02876992e-02,
                    1.84581297e-02,
                ],
                [
                    2.92695918e-02,
                    6.33243799e-02,
                    8.04641860e-02,
                    9.68530395e-02,
                    4.01675329e-02,
                    5.37861708e-02,
                    5.62332763e-02,
                    5.47512690e-02,
                    1.14313454e-01,
                    8.87141741e-02,
                    5.51365068e-02,
                    6.34530888e-02,
                    5.58837453e-02,
                    5.17867676e-02,
                    1.75140888e-02,
                ],
                [
                    2.76117514e-02,
                    6.41844460e-02,
                    8.83072363e-02,
                    9.81920130e-02,
                    4.06447266e-02,
                    5.24859771e-02,
                    5.56812634e-02,
                    5.66825173e-02,
                    1.14935578e-01,
                    8.62270814e-02,
                    5.89627860e-02,
                    6.01585961e-02,
                    5.43943312e-02,
                    5.27664145e-02,
                    1.80339359e-02,
                ],
                [
                    4.55190605e-02,
                    1.23504744e-01,
                    6.31111118e-02,
                    4.17455815e-02,
                    6.71836578e-02,
                    1.21421180e-01,
                    4.68796988e-02,
                    2.13966409e-01,
                    5.23277565e-02,
                    4.31326111e-02,
                    9.74844111e-03,
                    5.50657852e-02,
                    2.28663271e-02,
                    2.74040840e-02,
                    5.80642113e-02,
                ],
                [
                    1.77127374e-02,
                    1.15823873e-01,
                    7.24949145e-02,
                    5.57014738e-02,
                    5.77815115e-02,
                    4.00573632e-01,
                    1.26863097e-01,
                    1.84102374e-01,
                    1.01559988e-01,
                    1.21625711e-01,
                    2.57286919e-02,
                    5.51329473e-02,
                    1.83150224e-02,
                    1.99033312e-02,
                    4.65919405e-02,
                ],
                [
                    7.87452796e-03,
                    1.30542354e-01,
                    9.76530467e-02,
                    5.44756028e-02,
                    1.07494350e-01,
                    4.76297096e-01,
                    1.30635730e-01,
                    1.29260977e-01,
                    4.35271580e-02,
                    1.16218722e-01,
                    3.80001051e-02,
                    8.48773590e-02,
                    1.91909411e-02,
                    7.19048976e-03,
                    7.62703155e-02,
                ],
                [
                    1.30890397e-02,
                    1.26521264e-01,
                    1.56178831e-01,
                    6.92059728e-02,
                    1.59238430e-01,
                    3.74156481e-01,
                    9.34240399e-02,
                    1.31599576e-01,
                    3.97558461e-02,
                    1.06259429e-01,
                    7.07275694e-02,
                    6.70161003e-02,
                    2.13083230e-02,
                    4.26700552e-03,
                    6.21278092e-02,
                ],
                [
                    1.12121919e-02,
                    1.29750834e-01,
                    1.98096468e-01,
                    7.25404134e-02,
                    1.79567170e-01,
                    3.36764030e-01,
                    7.21105243e-02,
                    1.24657417e-01,
                    4.02441055e-02,
                    9.76995542e-02,
                    8.58065894e-02,
                    5.37411193e-02,
                    1.74258020e-02,
                    4.33662320e-03,
                    4.77155764e-02,
                ],
                [
                    1.03108396e-02,
                    1.30292806e-01,
                    2.20632668e-01,
                    7.38579068e-02,
                    1.84430439e-01,
                    3.21562339e-01,
                    6.36867632e-02,
                    1.16896847e-01,
                    4.05065880e-02,
                    9.00290369e-02,
                    9.35966376e-02,
                    4.87967467e-02,
                    1.52111116e-02,
                    4.74256226e-03,
                    4.46176945e-02,
                ],
                [
                    1.00592664e-02,
                    1.29907989e-01,
                    2.31809837e-01,
                    7.46272090e-02,
                    1.85303272e-01,
                    3.14663344e-01,
                    6.02779135e-02,
                    1.12402415e-01,
                    4.06291999e-02,
                    8.55696231e-02,
                    9.75423195e-02,
                    4.67323443e-02,
                    1.44511090e-02,
                    5.04185116e-03,
                    4.48151347e-02,
                ],
                [
                    1.31030153e-01,
                    3.34567995e-02,
                    4.94929071e-02,
                    1.30020064e-01,
                    1.00758336e-01,
                    2.23511078e-01,
                    1.53074901e-01,
                    7.92003807e-02,
                    4.75530940e-01,
                    9.11369062e-02,
                    1.93861038e-01,
                    2.68156812e-02,
                    1.86812581e-01,
                    2.85957388e-02,
                    3.83389544e-02,
                ],
                [
                    2.08690899e-01,
                    2.22975564e-01,
                    2.01470834e-02,
                    8.46321639e-02,
                    1.02934983e-01,
                    5.73348566e-01,
                    4.72739914e-02,
                    2.95948289e-01,
                    5.70657575e-01,
                    1.61901615e-02,
                    1.89428914e-01,
                    1.31983466e-01,
                    2.05563366e-01,
                    2.70351174e-02,
                    1.44732501e-02,
                ],
                [
                    1.83727647e-01,
                    2.15285427e-01,
                    2.95229664e-02,
                    7.35485119e-02,
                    9.11102616e-02,
                    1.08956321e00,
                    1.80128030e-01,
                    1.89540307e-01,
                    1.63483965e-01,
                    4.23285485e-02,
                    2.12547660e-01,
                    1.49284235e-01,
                    1.01318334e-01,
                    3.67877221e-02,
                    3.72032659e-02,
                ],
                [
                    1.62696106e-01,
                    2.22651521e-01,
                    4.68729804e-02,
                    1.75213126e-01,
                    1.37109099e-01,
                    1.14878177e00,
                    1.70241994e-01,
                    1.22589894e-01,
                    2.52989602e-01,
                    8.49382109e-02,
                    1.41906248e-01,
                    1.46673198e-01,
                    9.22477204e-02,
                    3.58142677e-02,
                    3.13827178e-02,
                ],
                [
                    1.76945602e-01,
                    2.06974203e-01,
                    9.67079167e-02,
                    2.25356711e-01,
                    2.15244723e-01,
                    1.01736314e00,
                    1.33952431e-01,
                    1.38972351e-01,
                    2.57512694e-01,
                    1.33546320e-01,
                    1.19004732e-01,
                    1.37112426e-01,
                    9.38492970e-02,
                    3.86824317e-02,
                    1.78200142e-02,
                ],
                [
                    1.74459726e-01,
                    2.10759539e-01,
                    1.36906376e-01,
                    2.33344971e-01,
                    2.43938160e-01,
                    9.50278264e-01,
                    1.21714452e-01,
                    1.28965883e-01,
                    2.60797540e-01,
                    1.43691799e-01,
                    1.05566832e-01,
                    1.22476748e-01,
                    1.05252358e-01,
                    3.75688869e-02,
                    1.07003001e-02,
                ],
                [
                    1.73637650e-01,
                    2.10793205e-01,
                    1.58467571e-01,
                    2.37048831e-01,
                    2.48592263e-01,
                    9.15971117e-01,
                    1.18626835e-01,
                    1.17965100e-01,
                    2.63760533e-01,
                    1.41710643e-01,
                    9.82716860e-02,
                    1.14969218e-01,
                    1.12527771e-01,
                    3.76243294e-02,
                    9.12505777e-03,
                ],
                [
                    1.73830044e-01,
                    2.09755401e-01,
                    1.68991063e-01,
                    2.39327940e-01,
                    2.48315000e-01,
                    8.98760969e-01,
                    1.17768543e-01,
                    1.11329969e-01,
                    2.65614224e-01,
                    1.39234806e-01,
                    9.45818872e-02,
                    1.11335212e-01,
                    1.16589571e-01,
                    3.79905707e-02,
                    8.91514683e-03,
                ],
                [
                    9.26454066e-01,
                    6.38855875e-01,
                    4.21319358e-01,
                    1.79943843e-01,
                    3.11304484e-02,
                    1.87329255e-01,
                    1.47447593e-01,
                    1.17698910e-01,
                    1.91043994e-01,
                    2.46418008e-01,
                    1.06673877e-01,
                    1.53158258e-01,
                    1.36744591e-01,
                    3.80844184e-01,
                    1.49904431e-01,
                ],
                [
                    4.74179240e-01,
                    5.50400254e-01,
                    4.65655469e-01,
                    5.04632657e-01,
                    1.39034165e-01,
                    7.37636338e-01,
                    5.62117094e-01,
                    1.63353578e-01,
                    9.85585836e-01,
                    6.14565368e-01,
                    2.24728007e-01,
                    1.45824731e-01,
                    5.47694577e-01,
                    5.01761979e-01,
                    1.47636846e-01,
                ],
                [
                    6.72023675e-01,
                    1.85966247e-01,
                    2.62832251e-01,
                    3.06822210e-01,
                    1.82682081e-01,
                    1.12373842e00,
                    3.48811862e-01,
                    3.14729488e-01,
                    9.38860768e-01,
                    3.48026214e-01,
                    1.61845874e-01,
                    1.79069251e-01,
                    4.95264141e-01,
                    4.19238174e-01,
                    9.52559806e-02,
                ],
                [
                    5.59203213e-01,
                    2.00428494e-01,
                    2.55943599e-01,
                    2.73965972e-01,
                    1.71434793e-01,
                    1.73716901e00,
                    4.59214613e-01,
                    2.67625422e-01,
                    4.56803123e-01,
                    3.57023355e-01,
                    1.52496550e-01,
                    1.65939820e-01,
                    4.13889405e-01,
                    3.54383108e-01,
                    4.63503606e-02,
                ],
                [
                    4.62597323e-01,
                    2.69023630e-01,
                    2.59106725e-01,
                    3.93478778e-01,
                    2.13155863e-01,
                    1.73348644e00,
                    4.08347091e-01,
                    1.70164570e-01,
                    6.63972600e-01,
                    4.81783044e-01,
                    8.84504943e-02,
                    7.99825428e-02,
                    4.05785779e-01,
                    4.37286861e-01,
                    8.89420795e-02,
                ],
                [
                    4.38002714e-01,
                    3.21722017e-01,
                    2.12020297e-01,
                    4.66967165e-01,
                    3.04185042e-01,
                    1.62204964e00,
                    3.82002047e-01,
                    1.34069738e-01,
                    7.05895905e-01,
                    6.26057896e-01,
                    4.75314559e-02,
                    8.59783925e-02,
                    4.09858051e-01,
                    4.89499864e-01,
                    9.91214509e-02,
                ],
                [
                    4.48419151e-01,
                    3.42292557e-01,
                    1.95040794e-01,
                    4.84747110e-01,
                    3.36990570e-01,
                    1.54289979e00,
                    3.85972715e-01,
                    1.20567067e-01,
                    7.19241938e-01,
                    6.63920059e-01,
                    3.83619000e-02,
                    7.79704502e-02,
                    4.37951102e-01,
                    4.99977747e-01,
                    9.89058194e-02,
                ],
                [
                    4.56856373e-01,
                    3.47945833e-01,
                    1.93607563e-01,
                    4.91481873e-01,
                    3.40634053e-01,
                    1.49659900e00,
                    3.88329643e-01,
                    1.18825632e-01,
                    7.26771488e-01,
                    6.69251679e-01,
                    3.52082740e-02,
                    7.33447198e-02,
                    4.54677321e-01,
                    5.02418614e-01,
                    1.01308440e-01,
                ],
                [
                    4.61455817e-01,
                    3.49767272e-01,
                    1.94792562e-01,
                    4.94530826e-01,
                    3.39094609e-01,
                    1.47242744e00,
                    3.89053698e-01,
                    1.19460086e-01,
                    7.30950347e-01,
                    6.68851649e-01,
                    3.37388871e-02,
                    7.12586137e-02,
                    4.63419021e-01,
                    5.03295792e-01,
                    1.03315279e-01,
                ],
                [
                    1.55385691e-01,
                    2.34697953e-01,
                    7.63966114e-01,
                    3.03825998e-01,
                    3.71506076e-02,
                    4.24745755e-01,
                    1.19922440e00,
                    1.14060964e00,
                    1.02442620e00,
                    1.06795685e00,
                    8.54842056e-01,
                    6.64783658e-01,
                    4.42815993e-01,
                    3.97895039e-01,
                    1.65740425e00,
                ],
                [
                    7.63425198e-01,
                    6.66793673e-01,
                    2.11022621e00,
                    4.30060123e-01,
                    1.36468186e-02,
                    2.91657502e-01,
                    1.45711098e00,
                    1.22988218e00,
                    1.77899207e00,
                    1.52833714e00,
                    7.05309099e-01,
                    1.17732544e00,
                    6.70122467e-01,
                    9.84776901e-01,
                    1.10709556e00,
                ],
                [
                    5.38195720e-01,
                    7.50513454e-01,
                    2.00575041e00,
                    4.88752468e-01,
                    9.39604234e-02,
                    4.33617754e-01,
                    1.99018819e00,
                    1.50717324e00,
                    2.47850173e00,
                    2.04821863e00,
                    1.26823956e00,
                    1.25395122e00,
                    7.92186959e-01,
                    9.99608541e-01,
                    1.29449637e00,
                ],
                [
                    7.95930368e-01,
                    5.04414627e-01,
                    1.77328308e00,
                    3.87578991e-01,
                    1.53210871e-01,
                    6.43157467e-01,
                    1.56569743e00,
                    7.65311059e-01,
                    2.19201079e00,
                    1.81063295e00,
                    1.23110073e00,
                    1.12163190e00,
                    6.56047003e-01,
                    1.06814145e00,
                    1.06972826e00,
                ],
                [
                    6.69417055e-01,
                    5.18959856e-01,
                    1.78786071e00,
                    3.88519787e-01,
                    1.35523899e-01,
                    1.06701697e00,
                    1.16882019e00,
                    6.99191934e-01,
                    1.85133116e00,
                    1.94968684e00,
                    1.36085853e00,
                    1.38243802e00,
                    7.97258569e-01,
                    7.39410939e-01,
                    1.44934461e00,
                ],
                [
                    5.24328863e-01,
                    5.47006165e-01,
                    1.83523282e00,
                    2.53374803e-01,
                    2.05757373e-01,
                    1.13835556e00,
                    1.20400992e00,
                    7.62286366e-01,
                    2.25817268e00,
                    2.24988226e00,
                    1.17011687e00,
                    1.18412898e00,
                    7.87787940e-01,
                    8.44074318e-01,
                    1.44944982e00,
                ],
                [
                    4.93969799e-01,
                    5.66602817e-01,
                    1.73937154e00,
                    2.40650627e-01,
                    2.99857097e-01,
                    1.03750989e00,
                    1.30688553e00,
                    7.67488098e-01,
                    2.35281508e00,
                    2.50052150e00,
                    1.01020887e00,
                    1.18375116e00,
                    7.63417659e-01,
                    9.08951744e-01,
                    1.35711084e00,
                ],
                [
                    5.10510237e-01,
                    5.87663492e-01,
                    1.69314196e00,
                    2.38866570e-01,
                    3.33667591e-01,
                    1.01024798e00,
                    1.37657881e00,
                    7.83644541e-01,
                    2.38828413e00,
                    2.54491213e00,
                    9.18408381e-01,
                    1.12638011e00,
                    8.01987300e-01,
                    9.24668474e-01,
                    1.30235187e00,
                ],
                [
                    5.23968337e-01,
                    5.91753645e-01,
                    1.68546560e00,
                    2.39636505e-01,
                    3.39810000e-01,
                    9.99568909e-01,
                    1.41895242e00,
                    7.96981286e-01,
                    2.40824620e00,
                    2.53419306e00,
                    8.78295300e-01,
                    1.09610536e00,
                    8.29086108e-01,
                    9.26152093e-01,
                    1.28937725e00,
                ],
                [
                    5.31575331e-01,
                    5.91889580e-01,
                    1.68629680e00,
                    2.40149073e-01,
                    3.39964249e-01,
                    9.94394263e-01,
                    1.44169989e00,
                    8.04796912e-01,
                    2.41905095e00,
                    2.52202592e00,
                    8.60755133e-01,
                    1.08341982e00,
                    8.43673187e-01,
                    9.25220791e-01,
                    1.28906196e00,
                ],
            ]
        )
        assert result_msd_per_particle.shape == expected_msd_per_particle.shape
        assert_allclose(
            result_msd_per_particle, expected_msd_per_particle, rtol=1e-5
        )

    def test_start_stop_step(self, u_nonlinear):
        msd = MSD(u_nonlinear, select="all", msd_type="xyz", non_linear=True)
        msd.run(start=3, stop=9, step=2)
        result_msd = msd.results.timeseries
        result_delta_t = msd.results.delta_t_values
        result_msd_per_particle = msd.results.msds_by_particle
        expected_msd = np.array([0.0, 0.02851658, 0.09466075, 0.09965838])
        expected_delta_t = np.array([0.0, 24.0, 96.0, 120.0])
        expected_msd_per_particle = np.array(
            [
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    0.01427892,
                    0.04496162,
                    0.02283019,
                    0.09793866,
                    0.03741678,
                    0.00456202,
                    0.01366753,
                    0.02593187,
                    0.03335009,
                    0.04255643,
                    0.03526032,
                    0.01883167,
                    0.01031271,
                    0.01167871,
                    0.01417117,
                ],
                [
                    0.01771274,
                    0.11582387,
                    0.07249491,
                    0.05570147,
                    0.05778151,
                    0.40057363,
                    0.1268631,
                    0.18410237,
                    0.10155999,
                    0.12162571,
                    0.02572869,
                    0.05513295,
                    0.01831502,
                    0.01990333,
                    0.04659194,
                ],
                [
                    0.01308904,
                    0.12652126,
                    0.15617883,
                    0.06920597,
                    0.15923843,
                    0.37415648,
                    0.09342404,
                    0.13159958,
                    0.03975585,
                    0.10625943,
                    0.07072757,
                    0.0670161,
                    0.02130832,
                    0.00426701,
                    0.06212781,
                ],
            ]
        )
        assert result_msd.shape == expected_msd.shape
        assert result_delta_t.shape == expected_delta_t.shape
        assert result_msd_per_particle.shape == expected_msd_per_particle.shape
        assert_allclose(result_msd, expected_msd, rtol=1e-5)
        assert_allclose(result_delta_t, expected_delta_t, rtol=1e-5)
        assert_allclose(
            result_msd_per_particle, expected_msd_per_particle, rtol=1e-5
        )
