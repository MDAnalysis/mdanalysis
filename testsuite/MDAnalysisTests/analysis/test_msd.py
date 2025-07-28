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

    def test_start_stop_step(self, u_nonlinear):
        msd = MSD(u_nonlinear, select="all", msd_type="xyz", non_linear=True)
        msd.run(start=3, stop=9, step=2)
        result_msd = msd.results.timeseries
        result_delta_t = msd.results.delta_t_values
        expected_msd = np.array([0.0, 0.02851658, 0.09466075, 0.09965838])
        expected_delta_t = np.array([0.0, 24.0, 96.0, 120.0])
        assert result_msd.shape == expected_msd.shape
        assert result_delta_t.shape == expected_delta_t.shape
        assert_allclose(result_msd, expected_msd, rtol=1e-5)
        assert_allclose(result_delta_t, expected_delta_t, rtol=1e-5)
