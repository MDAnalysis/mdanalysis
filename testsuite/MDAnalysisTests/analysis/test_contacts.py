# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
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
import warnings
import MDAnalysis as mda
import pytest
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.distances import distance_array

from numpy.testing import (
    assert_equal,
    assert_array_equal,
    assert_allclose,
)
import numpy as np

from MDAnalysisTests.datafiles import (
    PSF,
    DCD,
    TPR,
    XTC,
    contacts_villin_folded,
    contacts_villin_unfolded,
    contacts_file
)


def test_soft_cut_q():
    # just check some of the extremal points
    assert contacts.soft_cut_q([0], [0]) == .5
    assert_allclose(contacts.soft_cut_q([100], [0]), 0, rtol=0, atol=1.5e-7)
    assert_allclose(contacts.soft_cut_q([-100], [0]), 1, rtol=0, atol=1.5e-7)


def test_soft_cut_q_folded():
    u = mda.Universe(contacts_villin_folded)

    contacts_data = np.genfromtxt(contacts_file)
    # indices have been stored 1 indexed
    indices = contacts_data[:, :2].astype(int) - 1

    r = np.linalg.norm(u.atoms.positions[indices[:, 0]] -
                       u.atoms.positions[indices[:, 1]], axis=1)
    r0 = contacts_data[:, 2]

    beta = 5.0
    lambda_constant = 1.8
    Q = 1 / (1 + np.exp(beta * (r - lambda_constant * r0)))

    assert_allclose(Q.mean(), 1.0, rtol=0, atol=1.5e-3)


def test_soft_cut_q_unfolded():
    u = mda.Universe(contacts_villin_unfolded)

    contacts_data = np.genfromtxt(contacts_file)
    # indices have been stored 1 indexed
    indices = contacts_data[:, :2].astype(int) - 1

    r = np.linalg.norm(u.atoms.positions[indices[:, 0]] -
                       u.atoms.positions[indices[:, 1]], axis=1)
    r0 = contacts_data[:, 2]

    beta = 5.0
    lambda_constant = 1.8
    Q = 1 / (1 + np.exp(beta * (r - lambda_constant * r0)))

    assert_allclose(Q.mean(), 0.0, rtol=0, atol=1.5e-1)


@pytest.mark.parametrize('r, cutoff, expected_value', [
    ([1], 2, 1),
    ([2], 1, 0),
    ([2, 0.5], 1, 0.5),
    ([2, 3], [3, 4], 1),
    ([4, 5], [3, 4], 0)

])
def test_hard_cut_q(r, cutoff, expected_value):
    # just check some extremal points
    assert contacts.hard_cut_q(r, cutoff) == expected_value


@pytest.mark.parametrize('r, r0, radius, expected_value', [
    ([1], None, 2, 1),
    ([2], None, 1, 0),
    ([2, 0.5], None, 1, 0.5)

])
def test_radius_cut_q(r, r0, radius, expected_value):
    # check some extremal points
    assert contacts.radius_cut_q(r, r0, radius) == expected_value


def test_contact_matrix():
    d = np.arange(5)
    radius = np.ones(5) * 2.5

    out = contacts.contact_matrix(d, radius)
    assert_array_equal(out, [True, True, True, False, False])

    # check in-place update
    out = np.empty(out.shape)
    contacts.contact_matrix(d, radius, out=out)
    assert_array_equal(out, [True, True, True, False, False])


def test_new_selection():
    u = mda.Universe(PSF, DCD)
    selections = ('all', )
    sel = contacts._new_selections(u, selections, -1)[0]
    u.trajectory[-1]
    assert_array_equal(sel.positions, u.atoms.positions)


def soft_cut(ref, u, selA, selB, radius=4.5, beta=5.0, lambda_constant=1.8):
    """
        Reference implementation for testing
    """
    # reference groups A and B from selection strings
    refA, refB = ref.select_atoms(selA), ref.select_atoms(selB)

    # 2D float array, reference distances (r0)
    dref = distance_array(refA.positions, refB.positions)

    # 2D bool array, select reference distances that are less than the cutoff
    # radius
    mask = dref < radius

    # group A and B in a trajectory
    grA, grB = u.select_atoms(selA), u.select_atoms(selB)
    results = []

    for ts in u.trajectory:
        d = distance_array(grA.positions, grB.positions)
        r, r0 = d[mask], dref[mask]
        x = 1 / (1 + np.exp(beta * (r - lambda_constant * r0)))

        # average/normalize and append to results
        results.append((ts.time, x.sum() / mask.sum()))

    return np.asarray(results)


class TestContacts(object):
    sel_basic = "(resname ARG LYS) and (name NH* NZ)"
    sel_acidic = "(resname ASP GLU) and (name OE* OD*)"

    @staticmethod
    @pytest.fixture()
    def universe():
        return mda.Universe(PSF, DCD)

    def _run_Contacts(
        self, universe,
        start=None, stop=None, step=None, **kwargs
    ):
        acidic = universe.select_atoms(self.sel_acidic)
        basic = universe.select_atoms(self.sel_basic)
        return contacts.Contacts(
            universe,
            select=(self.sel_acidic, self.sel_basic),
            refgroup=(acidic, basic),
            radius=6.0,
            **kwargs).run(start=start, stop=stop, step=step)

    @pytest.mark.parametrize("seltxt", [sel_acidic, sel_basic])
    def test_select_valid_types(self, universe, seltxt):
        """Test if Contacts._get_atomgroup() can take both string and AtomGroup
         as selections.
        """
        ag = universe.select_atoms(seltxt)

        ag_from_string = contacts.Contacts._get_atomgroup(universe, seltxt)
        ag_from_ag = contacts.Contacts._get_atomgroup(universe, ag)

        assert ag_from_string == ag_from_ag

    def test_contacts_selections(self, universe):
        """Test if Contacts can take both string and AtomGroup as selections.
        """
        aga = universe.select_atoms(self.sel_acidic)
        agb = universe.select_atoms(self.sel_basic)

        cag = contacts.Contacts(
            universe, select=(aga, agb), refgroup=(aga, agb)
        )

        csel = contacts.Contacts(
            universe, select=(self.sel_acidic, self.sel_basic),
            refgroup=(aga, agb)
        )

        cag.run()
        csel.run()

        assert cag.grA == csel.grA
        assert cag.grB == csel.grB

    @pytest.mark.parametrize("ag", [1, [2], mda.Universe, "USE UPDATING AG"])
    def test_select_wrong_types(self, universe, ag):
        """Test that Contacts._get_atomgroup(u, sel) fails if sel is of the
        wrong type"""
        if ag == "USE UPDATING AG":
            ag = universe.select_atoms(self.sel_acidic, updating=True)

        with pytest.raises(
            TypeError, match="must be either string or a static AtomGroup"
        ) as te:
            contacts.Contacts._get_atomgroup(universe, ag)

    def test_startframe(self, universe):
        """test_startframe: TestContactAnalysis1: start frame set to 0 (resolution of
        Issue #624)

        """
        CA1 = self._run_Contacts(universe)
        assert len(CA1.results.timeseries) == universe.trajectory.n_frames

    def test_end_zero(self, universe):
        """test_end_zero: TestContactAnalysis1: stop frame 0 is not ignored"""
        CA1 = self._run_Contacts(universe, stop=0)
        assert len(CA1.results.timeseries) == 0

    def test_slicing(self, universe):
        start, stop, step = 10, 30, 5
        CA1 = self._run_Contacts(universe, start=start, stop=stop, step=step)
        frames = np.arange(universe.trajectory.n_frames)[start:stop:step]
        assert len(CA1.results.timeseries) == len(frames)

    def test_villin_folded(self):
        # one folded, one unfolded
        f = mda.Universe(contacts_villin_folded)
        u = mda.Universe(contacts_villin_unfolded)
        sel = "protein and not name H*"

        grF = f.select_atoms(sel)

        q = contacts.Contacts(u,
                              select=(sel, sel),
                              refgroup=(grF, grF),
                              method="soft_cut")
        q.run()

        results = soft_cut(f, u, sel, sel)
        assert_allclose(q.results.timeseries[:, 1], results[:, 1], rtol=0, atol=1.5e-7)

    def test_villin_unfolded(self):
        # both folded
        f = mda.Universe(contacts_villin_folded)
        u = mda.Universe(contacts_villin_folded)
        sel = "protein and not name H*"

        grF = f.select_atoms(sel)

        q = contacts.Contacts(u,
                              select=(sel, sel),
                              refgroup=(grF, grF),
                              method="soft_cut")
        q.run()

        results = soft_cut(f, u, sel, sel)
        assert_allclose(q.results.timeseries[:, 1], results[:, 1], rtol=0, atol=1.5e-7)

    def test_hard_cut_method(self, universe):
        ca = self._run_Contacts(universe)
        expected = [1., 0.58252427, 0.52427184, 0.55339806, 0.54368932,
                    0.54368932, 0.51456311, 0.46601942, 0.48543689, 0.52427184,
                    0.46601942, 0.58252427, 0.51456311, 0.48543689, 0.48543689,
                    0.48543689, 0.46601942, 0.51456311, 0.49514563, 0.49514563,
                    0.45631068, 0.47572816, 0.49514563, 0.50485437, 0.53398058,
                    0.50485437, 0.51456311, 0.51456311, 0.49514563, 0.49514563,
                    0.54368932, 0.50485437, 0.48543689, 0.55339806, 0.45631068,
                    0.46601942, 0.53398058, 0.53398058, 0.46601942, 0.52427184,
                    0.45631068, 0.46601942, 0.47572816, 0.46601942, 0.45631068,
                    0.47572816, 0.45631068, 0.48543689, 0.4368932, 0.4368932,
                    0.45631068, 0.50485437, 0.41747573, 0.4368932, 0.51456311,
                    0.47572816, 0.46601942, 0.46601942, 0.47572816, 0.47572816,
                    0.46601942, 0.45631068, 0.44660194, 0.47572816, 0.48543689,
                    0.47572816, 0.42718447, 0.40776699, 0.37864078, 0.42718447,
                    0.45631068, 0.4368932, 0.4368932, 0.45631068, 0.4368932,
                    0.46601942, 0.45631068, 0.48543689, 0.44660194, 0.44660194,
                    0.44660194, 0.42718447, 0.45631068, 0.44660194, 0.48543689,
                    0.48543689, 0.44660194, 0.4368932, 0.40776699, 0.41747573,
                    0.48543689, 0.45631068, 0.46601942, 0.47572816, 0.51456311,
                    0.45631068, 0.37864078, 0.42718447]
        assert len(ca.results.timeseries) == len(expected)
        assert_allclose(ca.results.timeseries[:, 1], expected, rtol=0, atol=1.5e-7)

    def test_radius_cut_method(self, universe):
        acidic = universe.select_atoms(self.sel_acidic)
        basic = universe.select_atoms(self.sel_basic)
        r = contacts.distance_array(acidic.positions, basic.positions)
        initial_contacts = contacts.contact_matrix(r, 6.0)
        expected = []
        for ts in universe.trajectory:
            r = contacts.distance_array(acidic.positions, basic.positions)
            expected.append(contacts.radius_cut_q(r[initial_contacts], None, radius=6.0))

        ca = self._run_Contacts(universe, method='radius_cut')
        assert_array_equal(ca.results.timeseries[:, 1], expected)

    @staticmethod
    def _is_any_closer(r, r0, dist=2.5):
        return np.any(r < dist)

    def test_own_method(self, universe):
        ca = self._run_Contacts(universe, method=self._is_any_closer)

        bound_expected = [1., 1., 0., 1., 1., 0., 0., 1., 0., 1., 1., 0., 0.,
                          1., 0., 0., 0., 0., 1., 1., 0., 0., 0., 1., 0., 1.,
                          0., 1., 1., 0., 1., 1., 1., 0., 0., 0., 0., 1., 0.,
                          0., 1., 0., 1., 1., 1., 0., 1., 0., 0., 1., 1., 1.,
                          0., 1., 0., 1., 1., 0., 0., 0., 1., 1., 1., 0., 0.,
                          1., 0., 1., 1., 1., 1., 1., 1., 0., 1., 1., 0., 1.,
                          0., 0., 1., 1., 0., 0., 1., 1., 1., 0., 1., 0., 0.,
                          1., 0., 1., 1., 1., 1., 1.]
        assert_array_equal(ca.results.timeseries[:, 1], bound_expected)

    @staticmethod
    def _weird_own_method(r, r0):
        return 'aaa'

    def test_own_method_no_array_cast(self, universe):
        with pytest.raises(ValueError):
            self._run_Contacts(universe, method=self._weird_own_method, stop=2)

    def test_non_callable_method(self, universe):
        with pytest.raises(ValueError):
            self._run_Contacts(universe, method=2, stop=2)

    @pytest.mark.parametrize("pbc,expected", [
    (True, [1., 0.43138152, 0.3989021, 0.43824337, 0.41948765,
            0.42223239, 0.41354071, 0.43641354, 0.41216834, 0.38334858]),
    (False, [1., 0.42327791, 0.39192399, 0.40950119, 0.40902613,
             0.42470309, 0.41140143, 0.42897862, 0.41472684, 0.38574822])
    ])
    def test_distance_box(self, pbc, expected):
        u = mda.Universe(TPR, XTC)
        sel_basic = "(resname ARG LYS)"
        sel_acidic = "(resname ASP GLU)"
        acidic = u.select_atoms(sel_acidic)
        basic = u.select_atoms(sel_basic)
        
        r = contacts.Contacts(u, select=(sel_acidic, sel_basic),
                        refgroup=(acidic, basic), radius=6.0, pbc=pbc)
        r.run()
        assert_allclose(r.results.timeseries[:, 1], expected,rtol=0, atol=1.5e-7)

    def test_warn_deprecated_attr(self, universe):
        """Test for warning message emitted on using deprecated `timeseries`
        attribute"""
        CA1 = self._run_Contacts(universe, stop=1)
        wmsg = "The `timeseries` attribute was deprecated in MDAnalysis"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(CA1.timeseries, CA1.results.timeseries)


def test_q1q2():
    u = mda.Universe(PSF, DCD)
    q1q2 = contacts.q1q2(u, 'name CA', radius=8)
    q1q2.run()

    q1_expected = [1., 0.98092643, 0.97366031, 0.97275204, 0.97002725,
                   0.97275204, 0.96276113, 0.96730245, 0.9582198, 0.96185286,
                   0.95367847, 0.96276113, 0.9582198, 0.95186194, 0.95367847,
                   0.95095368, 0.94187103, 0.95186194, 0.94277929, 0.94187103,
                   0.9373297, 0.93642144, 0.93097184, 0.93914623, 0.93278837,
                   0.93188011, 0.9373297, 0.93097184, 0.93188011, 0.92643052,
                   0.92824705, 0.92915531, 0.92643052, 0.92461399, 0.92279746,
                   0.92643052, 0.93278837, 0.93188011, 0.93369664, 0.9346049,
                   0.9373297, 0.94096276, 0.9400545, 0.93642144, 0.9373297,
                   0.9373297, 0.9400545, 0.93006358, 0.9400545, 0.93823797,
                   0.93914623, 0.93278837, 0.93097184, 0.93097184, 0.92733878,
                   0.92824705, 0.92279746, 0.92824705, 0.91825613, 0.92733878,
                   0.92643052, 0.92733878, 0.93278837, 0.92733878, 0.92824705,
                   0.93097184, 0.93278837, 0.93914623, 0.93097184, 0.9373297,
                   0.92915531, 0.93188011, 0.93551317, 0.94096276, 0.93642144,
                   0.93642144, 0.9346049, 0.93369664, 0.93369664, 0.93278837,
                   0.93006358, 0.93278837, 0.93006358, 0.9346049, 0.92824705,
                   0.93097184, 0.93006358, 0.93188011, 0.93278837, 0.93006358,
                   0.92915531, 0.92824705, 0.92733878, 0.92643052, 0.93188011,
                   0.93006358, 0.9346049, 0.93188011]
    assert_allclose(q1q2.results.timeseries[:, 1], q1_expected, rtol=0, atol=1.5e-7)

    q2_expected = [0.94649446, 0.94926199, 0.95295203, 0.95110701, 0.94833948,
                   0.95479705, 0.94926199, 0.9501845, 0.94926199, 0.95387454,
                   0.95202952, 0.95110701, 0.94649446, 0.94095941, 0.94649446,
                   0.9400369, 0.94464945, 0.95202952, 0.94741697, 0.94649446,
                   0.94188192, 0.94188192, 0.93911439, 0.94464945, 0.9400369,
                   0.94095941, 0.94372694, 0.93726937, 0.93819188, 0.93357934,
                   0.93726937, 0.93911439, 0.93911439, 0.93450185, 0.93357934,
                   0.93265683, 0.93911439, 0.94372694, 0.93911439, 0.94649446,
                   0.94833948, 0.95110701, 0.95110701, 0.95295203, 0.94926199,
                   0.95110701, 0.94926199, 0.94741697, 0.95202952, 0.95202952,
                   0.95202952, 0.94741697, 0.94741697, 0.94926199, 0.94280443,
                   0.94741697, 0.94833948, 0.94833948, 0.9400369, 0.94649446,
                   0.94741697, 0.94926199, 0.95295203, 0.94926199, 0.9501845,
                   0.95664207, 0.95756458, 0.96309963, 0.95756458, 0.96217712,
                   0.95756458, 0.96217712, 0.96586716, 0.96863469, 0.96494465,
                   0.97232472, 0.97140221, 0.9695572, 0.97416974, 0.9695572,
                   0.96217712, 0.96771218, 0.9704797, 0.96771218, 0.9695572,
                   0.97140221, 0.97601476, 0.97693727, 0.98154982, 0.98431734,
                   0.97601476, 0.9797048, 0.98154982, 0.98062731, 0.98431734,
                   0.98616236, 0.9898524, 1.]
    assert_allclose(q1q2.results.timeseries[:, 2], q2_expected, rtol=0, atol=1.5e-7)
