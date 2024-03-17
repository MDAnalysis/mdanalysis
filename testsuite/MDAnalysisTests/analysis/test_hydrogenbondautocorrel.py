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
import pytest

from MDAnalysisTests.datafiles import (
    TRZ, TRZ_psf,
    waterPSF, waterDCD,
    XYZ_mini,
)
from numpy.testing import assert_almost_equal
import numpy as np
from unittest import mock

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import (HydrogenBondAutoCorrel as HBAC,
                                               find_hydrogen_donors)


class TestHydrogenBondAutocorrel(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(TRZ_psf, TRZ)

    @pytest.fixture()
    def hydrogens(self, u):
        return u.atoms.select_atoms('name Hn')

    @pytest.fixture()
    def nitrogens(self, u):
        return u.atoms.select_atoms('name N')

    @pytest.fixture()
    def oxygens(self, u):
        return u.atoms.select_atoms('name O')

    # regression tests for different conditions
    def test_continuous(self, u, hydrogens, oxygens, nitrogens):
        hbond = HBAC(u,
                     hydrogens=hydrogens,
                     acceptors=oxygens,
                     donors=nitrogens,
                     bond_type='continuous',
                     sample_time=0.06,
        )
        hbond.run()

        assert_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.83137828,
                       0.74486804,  0.67741936,  0.60263932],
                     dtype=np.float32)
        )

    def test_continuous_excl(self, u, hydrogens, oxygens, nitrogens):
        hbond = HBAC(u,
                     hydrogens=hydrogens,
                     acceptors=oxygens,
                     donors=nitrogens,
                     bond_type='continuous',
                     exclusions=(np.arange(len(hydrogens)), np.array(
                         range(len(oxygens)))),
                     sample_time=0.06,
        )
        hbond.run()

        assert_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.83137828,
                       0.74486804,  0.67741936,  0.60263932],
                     dtype=np.float32)
        )


    def test_intermittent(self, u, hydrogens, oxygens, nitrogens):
        hbond = HBAC(u,
                     hydrogens=hydrogens,
                     acceptors=oxygens,
                     donors=nitrogens,
                     bond_type='intermittent',
                     sample_time=0.06,
        )
        hbond.run()

        assert_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.84310848,
                       0.79325515,  0.76392961,  0.72287393],
                     dtype=np.float32)
        )


    def test_intermittent_timecut(self, u, hydrogens, oxygens, nitrogens):
        hbond = HBAC(u,
                     hydrogens=hydrogens,
                     acceptors=oxygens,
                     donors=nitrogens,
                     bond_type='intermittent',
                     time_cut=0.01,  # time cut at traj.dt == continuous
                     sample_time=0.06,
        )
        hbond.run()

        assert_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.83137828,
                       0.74486804,  0.67741936,  0.60263932],
                     dtype=np.float32)
        )

    def test_intermittent_excl(self, u, hydrogens, oxygens, nitrogens):
        hbond = HBAC(u,
                     hydrogens=hydrogens,
                     acceptors=oxygens,
                     donors=nitrogens,
                     bond_type='intermittent',
                     exclusions=(np.arange(len(hydrogens)), np.array(
                         range(len(oxygens)))),
                     sample_time=0.06,
        )
        hbond.run()

        assert_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.84310848,
                       0.79325515,  0.76392961,  0.72287393],
                     dtype=np.float32)
        )

    # For `solve` the test trajectories aren't long enough
    # So spoof the results and check that solver finds solution
    def test_solve_continuous(self, u, hydrogens, oxygens, nitrogens):
        hbond = HBAC(u,
                     hydrogens=hydrogens,
                     acceptors=oxygens,
                     donors=nitrogens,
                     bond_type='continuous',
                     sample_time=0.06,
        )

        def actual_function_cont(t):
            A1 = 0.75
            A2 = 0.25
            tau1 = 0.5
            tau2 = 0.1
            return A1 * np.exp(-t/tau1) + A2 * np.exp(-t/tau2)
        hbond.solution['time'] = time = np.arange(0, 0.06, 0.001)
        hbond.solution['results'] = actual_function_cont(time)

        hbond.solve()

        assert_almost_equal(
            hbond.solution['fit'],
            np.array([0.75, 0.5, 0.1]),
        )

    def test_solve_intermittent(self, u, hydrogens, oxygens, nitrogens):
        hbond = HBAC(u,
                     hydrogens=hydrogens,
                     acceptors=oxygens,
                     donors=nitrogens,
                     bond_type='intermittent',
                     sample_time=0.06,
        )

        def actual_function_int(t):
            A1 = 0.33
            A2 = 0.33
            A3 = 0.34
            tau1 = 5
            tau2 = 1
            tau3 = 0.1
            return A1 * np.exp(-t/tau1) + A2 * np.exp(-t/tau2) + A3 * np.exp(-t/tau3)
        hbond.solution['time'] = time = np.arange(0, 6.0, 0.01)
        hbond.solution['results'] = actual_function_int(time)

        hbond.solve()

        assert_almost_equal(
            hbond.solution['fit'],
            np.array([0.33, 0.33, 5, 1, 0.1]),
        )

    # setup errors
    def test_wronglength_DA(self, u, hydrogens, oxygens, nitrogens):
        with pytest.raises(ValueError):
            HBAC(u,
                 hydrogens=hydrogens[:-1],
                 acceptors=oxygens,
                 donors=nitrogens,
                 bond_type='intermittent',
                 exclusions=(np.arange(len(hydrogens)), np.array(
                         range(len(oxygens)))),
                 sample_time=0.06,
        )

    def test_exclusions(self, u, hydrogens, oxygens, nitrogens):
        excl_list = (np.array(range(len(hydrogens))), np.array(
                         range(len(oxygens))))
        excl_list2 = excl_list[0], excl_list[1][:-1]
        with pytest.raises(ValueError):
            HBAC(u,
                 hydrogens=hydrogens,
                 acceptors=oxygens,
                 donors=nitrogens,
                 bond_type='intermittent',
                 exclusions=excl_list2,
                 sample_time=0.06,
        )

    def test_bond_type_VE(self, u, hydrogens, oxygens, nitrogens):
        with pytest.raises(ValueError):
            HBAC(u,
                 hydrogens=hydrogens,
                 acceptors=oxygens,
                 donors=nitrogens,
                 bond_type='marzipan',
                 exclusions=(np.arange(len(hydrogens)), np.array(range(
                     len(oxygens)))),
                 sample_time=0.06,
        )

    def test_solve_before_run_VE(self, u, hydrogens, oxygens, nitrogens):
        hbond = HBAC(u,
                     hydrogens=hydrogens,
                     acceptors=oxygens,
                     donors=nitrogens,
                     bond_type='continuous',
                     sample_time=0.06,
        )
        with pytest.raises(ValueError):
            hbond.solve()

    @mock.patch('MDAnalysis.coordinates.TRZ.TRZReader._read_frame')
    def test_unslicable_traj_VE(self, mock_read, u, hydrogens, oxygens, nitrogens):
        mock_read.side_effect = TypeError

        with pytest.raises(ValueError):
            HBAC(
                u,
                hydrogens=hydrogens,
                acceptors=oxygens,
                donors=nitrogens,
                bond_type='continuous',
                sample_time=0.06
        )

    def test_repr(self, u, hydrogens, oxygens, nitrogens):
        hbond = HBAC(u,
                     hydrogens=hydrogens,
                     acceptors=oxygens,
                     donors=nitrogens,
                     bond_type='continuous',
                     sample_time=0.06,
        )
        assert isinstance(repr(hbond), str)


def test_find_donors():
    u = mda.Universe(waterPSF, waterDCD)

    H = u.select_atoms('name H*')

    D = find_hydrogen_donors(H)

    assert len(H) == len(D)
    # check each O is bonded to the corresponding H
    for h_atom, o_atom in zip(H, D):
        assert o_atom in h_atom.bonded_atoms


def test_donors_nobonds():
    u = mda.Universe(XYZ_mini)

    with pytest.raises(mda.NoDataError):
        find_hydrogen_donors(u.atoms)
