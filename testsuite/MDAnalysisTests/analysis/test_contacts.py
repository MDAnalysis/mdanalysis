# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.distances import distance_array

from numpy.testing import (dec, assert_almost_equal, assert_equal, raises,
                           assert_array_equal, assert_array_almost_equal)
import numpy as np

import os

from MDAnalysisTests.datafiles import (PSF,
                                       DCD,
                                       contacts_villin_folded,
                                       contacts_villin_unfolded,
                                       contacts_file, )

from MDAnalysisTests import parser_not_found, tempdir


def test_soft_cut_q():
    # just check some of the extremal points
    assert_equal(contacts.soft_cut_q([0], [0]), .5)
    assert_almost_equal(contacts.soft_cut_q([100], [0]), 0)
    assert_almost_equal(contacts.soft_cut_q([-100], [0]), 1)


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

    assert_almost_equal(Q.mean(), 1.0, decimal=3)


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

    assert_almost_equal(Q.mean(), 0.0, decimal=1)


def test_hard_cut_q():
    # just check some extremal points
    assert_equal(contacts.hard_cut_q([1], 2), 1)
    assert_equal(contacts.hard_cut_q([2], 1), 0)
    assert_equal(contacts.hard_cut_q([2, 0.5], 1), 0.5)
    assert_equal(contacts.hard_cut_q([2, 3], [3, 4]), 1)
    assert_equal(contacts.hard_cut_q([4, 5], [3, 4]), 0)


def test_radius_cut_q():
    # check some extremal points
    assert_equal(contacts.radius_cut_q([1], None, 2), 1)
    assert_equal(contacts.radius_cut_q([2], None, 1), 0)
    assert_equal(contacts.radius_cut_q([2, 0.5], None, 1), 0.5)


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
    @dec.skipif(
        parser_not_found('DCD'),
        'DCD parser not available. Are you using python 3?')
    def __init__(self):
        self.universe = mda.Universe(PSF, DCD)
        self.trajectory = self.universe.trajectory

        self.sel_basic = "(resname ARG LYS) and (name NH* NZ)"
        self.sel_acidic = "(resname ASP GLU) and (name OE* OD*)"

    def tearDown(self):
        # reset trajectory
        self.universe.trajectory[0]
        del self.universe

    def _run_Contacts(self, **kwargs):
        acidic = self.universe.select_atoms(self.sel_acidic)
        basic = self.universe.select_atoms(self.sel_basic)
        Contacts = contacts.Contacts(
            self.universe,
            selection=(self.sel_acidic, self.sel_basic),
            refgroup=(acidic, basic),
            radius=6.0,
            **kwargs)
        Contacts.run()
        return Contacts

    def test_startframe(self):
        """test_startframe: TestContactAnalysis1: start frame set to 0 (resolution of
        Issue #624)

        """
        CA1 = self._run_Contacts()
        assert_equal(len(CA1.timeseries), self.universe.trajectory.n_frames)

    def test_end_zero(self):
        """test_end_zero: TestContactAnalysis1: stop frame 0 is not ignored"""
        CA1 = self._run_Contacts(stop=0)
        assert_equal(len(CA1.timeseries), 0)

    def test_slicing(self):
        start, stop, step = 10, 30, 5
        CA1 = self._run_Contacts(start=start, stop=stop, step=step)
        frames = np.arange(self.universe.trajectory.n_frames)[start:stop:step]
        assert_equal(len(CA1.timeseries), len(frames))

    @staticmethod
    def test_villin_folded():
        # one folded, one unfolded
        f = mda.Universe(contacts_villin_folded)
        u = mda.Universe(contacts_villin_unfolded)
        sel = "protein and not name H*"

        grF = f.select_atoms(sel)

        q = contacts.Contacts(u,
                              selection=(sel, sel),
                              refgroup=(grF, grF),
                              method="soft_cut")
        q.run()

        results = soft_cut(f, u, sel, sel)
        assert_almost_equal(q.timeseries[:, 1], results[:, 1])

    @staticmethod
    def test_villin_unfolded():

        # both folded
        f = mda.Universe(contacts_villin_folded)
        u = mda.Universe(contacts_villin_folded)
        sel = "protein and not name H*"

        grF = f.select_atoms(sel)

        q = contacts.Contacts(u,
                              selection=(sel, sel),
                              refgroup=(grF, grF),
                              method="soft_cut")
        q.run()

        results = soft_cut(f, u, sel, sel)
        assert_almost_equal(q.timeseries[:, 1], results[:, 1])

    def test_hard_cut_method(self):
        ca = self._run_Contacts()
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
        assert_equal(len(ca.timeseries), len(expected))
        assert_array_almost_equal(ca.timeseries[:, 1], expected)

    @staticmethod
    def _is_any_closer(r, r0, dist=2.5):
        return np.any(r < dist)

    def test_own_method(self):
        ca = self._run_Contacts(method=self._is_any_closer)

        bound_expected = [1., 1., 0., 1., 1., 0., 0., 1., 0., 1., 1., 0., 0.,
                          1., 0., 0., 0., 0., 1., 1., 0., 0., 0., 1., 0., 1.,
                          0., 1., 1., 0., 1., 1., 1., 0., 0., 0., 0., 1., 0.,
                          0., 1., 0., 1., 1., 1., 0., 1., 0., 0., 1., 1., 1.,
                          0., 1., 0., 1., 1., 0., 0., 0., 1., 1., 1., 0., 0.,
                          1., 0., 1., 1., 1., 1., 1., 1., 0., 1., 1., 0., 1.,
                          0., 0., 1., 1., 0., 0., 1., 1., 1., 0., 1., 0., 0.,
                          1., 0., 1., 1., 1., 1., 1.]
        assert_array_equal(ca.timeseries[:, 1], bound_expected)

    @staticmethod
    def _weird_own_method(r, r0):
        return 'aaa'

    @raises(ValueError)
    def test_own_method_no_array_cast(self):
        self._run_Contacts(method=self._weird_own_method, stop=2)

    @raises(ValueError)
    def test_non_callable_method(self):
        self._run_Contacts(method=2, stop=2)

    def test_save(self):
        with tempdir.in_tempdir():
            ca = self._run_Contacts()
            ca.save('testfile.npy')
            saved = np.genfromtxt('testfile.npy')
            assert_array_almost_equal(ca.timeseries, saved)


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
    assert_array_almost_equal(q1q2.timeseries[:, 1], q1_expected)

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
    assert_array_almost_equal(q1q2.timeseries[:, 2], q2_expected)
