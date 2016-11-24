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
import six
from six.moves import zip, range
from MDAnalysisTests.datafiles import TRZ, TRZ_psf, PRM, TRJ
from MDAnalysisTests import module_not_found, tempdir
from numpy.testing import assert_, assert_array_almost_equal, assert_raises, assert_, dec
import numpy as np
import mock

import MDAnalysis as mda
from MDAnalysis.analysis.hbonds import HydrogenBondAutoCorrel as HBAC


class TestHydrogenBondAutocorrel(object):
    def setUp(self):
        u = self.u = mda.Universe(TRZ_psf, TRZ)
        self.H = u.atoms.select_atoms('name Hn')
        self.N = u.atoms.select_atoms('name N')
        self.O = u.atoms.select_atoms('name O')
        self.excl_list = (np.array(range(len(self.H))),
                          np.array(range(len(self.O))))

    def tearDown(self):
        del self.H
        del self.N
        del self.O
        del self.u
        del self.excl_list

    # regression tests for different conditions
    def test_continuous(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='continuous',
                     sample_time=0.06,
        )
        hbond.run()

        assert_array_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.83137828,
                       0.74486804,  0.67741936,  0.60263932],
                     dtype=np.float32)
        )


    def test_continuous_excl(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='continuous',
                     exclusions=self.excl_list,
                     sample_time=0.06,
        )
        hbond.run()

        assert_array_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.83137828,
                       0.74486804,  0.67741936,  0.60263932],
                     dtype=np.float32)
        )


    def test_intermittent(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='intermittent',
                     sample_time=0.06,
        )
        hbond.run()

        assert_array_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.84310848,
                       0.79325515,  0.76392961,  0.72287393],
                     dtype=np.float32)
        )


    def test_intermittent_timecut(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='intermittent',
                     time_cut=0.01,  # time cut at traj.dt == continuous
                     sample_time=0.06,
        )
        hbond.run()

        assert_array_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.83137828,
                       0.74486804,  0.67741936,  0.60263932],
                     dtype=np.float32)
        )

    def test_intermittent_excl(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='intermittent',
                     exclusions=self.excl_list,
                     sample_time=0.06,
        )
        hbond.run()

        assert_array_almost_equal(
            hbond.solution['results'],
            np.array([ 1.        ,  0.92668623,  0.84310848,
                       0.79325515,  0.76392961,  0.72287393],
                     dtype=np.float32)
        )

    # For `solve` the test trajectories aren't long enough
    # So spoof the results and check that solver finds solution
    @dec.skipif(module_not_found('scipy'))
    def test_solve_continuous(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
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

        assert_array_almost_equal(
            hbond.solution['fit'],
            np.array([0.75, 0.5, 0.1]),
        )

    @dec.skipif(module_not_found('scipy'))
    def test_solve_intermittent(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
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

        assert_array_almost_equal(
            hbond.solution['fit'],
            np.array([0.33, 0.33, 5, 1, 0.1]),
        )

    def test_save(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='continuous',
                     sample_time=0.06,
        )
        hbond.run()

        with tempdir.in_tempdir():
            hbond.save_results('hbondout.npz')

            loaded = np.load('hbondout.npz')
            assert_('time' in loaded)
            assert_('results' in loaded)

    # setup errors
    def test_wronglength_DA(self):
        assert_raises(ValueError,
                      HBAC, self.u,
                      hydrogens=self.H[:-1],
                      acceptors=self.O,
                      donors=self.N,
                      bond_type='intermittent',
                      exclusions=self.excl_list,
                      sample_time=0.06,
        )

    def test_exclusions(self):
        excl_list2 = self.excl_list[0], self.excl_list[1][:-1]
        assert_raises(ValueError,
                      HBAC, self.u,
                      hydrogens=self.H,
                      acceptors=self.O,
                      donors=self.N,
                      bond_type='intermittent',
                      exclusions=excl_list2,
                      sample_time=0.06,
        )

    def test_bond_type_VE(self):
        assert_raises(ValueError,
                      HBAC, self.u,
                      hydrogens=self.H,
                      acceptors=self.O,
                      donors=self.N,
                      bond_type='marzipan',
                      exclusions=self.excl_list,
                      sample_time=0.06,
        )

    @dec.skipif(module_not_found('scipy'))
    def test_solve_before_run_VE(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='continuous',
                     sample_time=0.06,
        )
        assert_raises(ValueError, hbond.solve)

    @mock.patch('MDAnalysis.coordinates.TRZ.TRZReader._read_frame')
    def test_unslicable_traj_VE(self, mock_read):
        mock_read.side_effect = TypeError

        assert_raises(ValueError, HBAC,
                      self.u,
                      hydrogens=self.H,
                      acceptors=self.O,
                      donors=self.N,
                      bond_type='continuous',
                      sample_time=0.06
        )

    def test_save_without_run_VE(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='continuous',
                     sample_time=0.06,
        )
        assert_raises(ValueError, hbond.save_results)

    def test_repr(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='continuous',
                     sample_time=0.06,
        )
        assert_(isinstance(repr(hbond), six.string_types))
