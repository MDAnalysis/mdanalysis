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
from six.moves import zip, range
from MDAnalysisTests.datafiles import TRZ, TRZ_psf, PRM, TRJ
from MDAnalysisTests import module_not_found, tempdir
from numpy.testing import assert_, assert_array_almost_equal, assert_raises, dec
import numpy as np

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

    @dec.skipif(module_not_found('scipy'))
    def test_solve_continuous(self):
        hbond = HBAC(self.u,
                     hydrogens=self.H,
                     acceptors=self.O,
                     donors=self.N,
                     bond_type='continuous',
                     sample_time=0.06,
        )
        hbond.run()
        hbond.solve()

        assert_array_almost_equal(
            hbond.solution['fit'],
            np.array([ 0.52727374,  0.10231721,  0.10231605])
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
        hbond.run()
        hbond.solve()

        assert_array_almost_equal(
            hbond.solution['fit'],
            np.array([  6.13695140e-01,   3.30614269e-06,   5.90645176e+00,
                        2.00072251e-01,   4.15492813e-02])
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

    def test_unslicable_traj_VE(self):
        u = mda.Universe(PRM, TRJ)
        H = u.atoms[:10]
        O = u.atoms[10:20]
        N = u.atoms[20:30]
        assert_raises(ValueError, HBAC,
                      u,
                      hydrogens=H,
                      acceptors=O,
                      donors=N,
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
