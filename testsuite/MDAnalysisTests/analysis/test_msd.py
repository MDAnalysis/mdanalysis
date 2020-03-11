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
from __future__ import division, absolute_import, print_function



import MDAnalysis as mda
from MDAnalysis.analysis.msd import MeanSquaredDisplacement as MSD

from numpy.testing import (assert_array_less,
                           assert_almost_equal, assert_equal)
import numpy as np

from scipy import fft,ifft

from MDAnalysisTests.datafiles import PSF, DCD, DCD

import pytest

SELECTION = 'backbone and name CA and resid 1-10'
NSTEP = 1000

@pytest.fixture(scope='module')
def u():
    return mda.Universe(PSF, DCD)

@pytest.fixture(scope='module')
def msd(u):
    m = MSD(u, SELECTION, msd_type='xyz', position_treatment='atom', mass_weighted=False)
    m.run()
    return m

@pytest.fixture(scope='module')
def msd_fft(u):
    m = MSD(u, SELECTION, msd_type='xyz', position_treatment='atom', mass_weighted=False, fft=True)
    m.run()
    return m


@pytest.fixture(scope='module')
def dimension_list():
    dimensions = ['xyz', 'xy', 'xz', 'yz', 'x', 'y', 'z']
    return dimensions


@pytest.fixture(scope='module')
def step_traj():
    x = np.arange(NSTEP)
    traj = np.vstack([x,x,x]).T
    traj_reshape = traj.reshape([NSTEP,1,3])
    u = mda.Universe.empty(1)
    u.load_new(traj_reshape)
    return u

def characteristic_poly(n,d):
    x = np.arange(1,n+1)
    y = d*((x-1)*(x-1))
    return y
    


def test_fft_vs_simple_default(msd, msd_fft):
    timeseries_simple = msd.timeseries
    print(timeseries_simple)
    timeseries_fft = msd_fft.timeseries
    print(timeseries_fft)
    assert_almost_equal(timeseries_simple, timeseries_fft, decimal=5)


def test_fft_vs_simple_all_dims(dimension_list, u):
    for dim in dimension_list:
        print(dim)
        m_simple = MSD(u, SELECTION, msd_type=dim, position_treatment='atom', mass_weighted=False, fft=False)
        m_simple.run()
        timeseries_simple = m_simple.timeseries
        print(timeseries_simple)
        m_fft = MSD(u,SELECTION, msd_type=dim, position_treatment='atom', mass_weighted=False, fft=True)
        m_fft.run()
        timeseries_fft = m_fft.timeseries
        print(timeseries_fft)
        assert_almost_equal(timeseries_simple, timeseries_fft, decimal=5)

def test_simple_step_traj_3d(step_traj): # this should fit the polynomial 3x^2 - 6x +3
    m_simple = MSD(step_traj, 'all' , msd_type='xyz', position_treatment='atom', mass_weighted=False, fft=False)
    m_simple.run()
    poly = characteristic_poly(NSTEP,3)
    for i in range(len(poly)):
        if poly[i] != m_simple.timeseries[i]:
            print('MISMATCH')
            print(i)
            print(poly[i])
            print(m_simple.timeseries[i])
          # for some reason, this is much more prone to roundoff error?
    raise Exception

def test_fft_step_traj_3d(step_traj): # this should fit the polynomial 3(x-1)**2
    m_fft = MSD(step_traj, 'all' , msd_type='xyz', position_treatment='atom', mass_weighted=False, fft=True)
    m_fft.run()
    poly3 = characteristic_poly(NSTEP,3)
    assert_almost_equal(m_fft.timeseries, poly3)

def test_fft_step_traj_2d(step_traj): # this should fit the polynomial 2(x-1)**2
    m_fft = MSD(step_traj, 'all' , msd_type='xy', position_treatment='atom', mass_weighted=False, fft=True)
    m_fft.run()
    poly2 = characteristic_poly(NSTEP,2)
    assert_almost_equal(m_fft.timeseries, poly2)

def test_fft_step_traj_1d(step_traj): # this should fit the polynomial (x-1)**2
    m_fft = MSD(step_traj, 'all' , msd_type='x', position_treatment='atom', mass_weighted=False, fft=True)
    m_fft.run()
    poly1 = characteristic_poly(NSTEP,1)
    assert_almost_equal(m_fft.timeseries, poly1)

