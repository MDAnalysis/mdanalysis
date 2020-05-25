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
from MDAnalysis.analysis.msd import EinsteinMSD as MSD

from numpy.testing import (assert_array_less,
                           assert_almost_equal, assert_equal)
import numpy as np

from MDAnalysisTests.datafiles import PSF, DCD, RANDOM_WALK, RANDOM_WALK_TOPO

import pytest
import tidynamics

@pytest.fixture(scope='module')
def SELECTION():
    selection = 'backbone and name CA and resid 1-10'
    return selection

@pytest.fixture(scope='module')
def NSTEP():
    nstep = 5000
    return nstep

#universe
@pytest.fixture(scope='module')
def u():
    return mda.Universe(PSF, DCD)

@pytest.fixture(scope='module')
def random_walk_u():
    #100x100
    return mda.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)

#non fft msd
@pytest.fixture(scope='module')
def msd(u, SELECTION):
    m = MSD(u, SELECTION, msd_type='xyz', fft=False)
    m.run()
    return m

#fft msd
@pytest.fixture(scope='module')
def msd_fft(u, SELECTION):
    m = MSD(u, SELECTION, msd_type='xyz', fft=True)
    m.run()
    return m

#all possible dimensions
@pytest.fixture(scope='module')
def dimension_list():
    dimensions = ['xyz', 'xy', 'xz', 'yz', 'x', 'y', 'z']
    return dimensions

@pytest.fixture(scope='module')
def step_traj(NSTEP): # constant velocity
    x = np.arange(NSTEP)
    traj = np.vstack([x,x,x]).T
    traj_reshape = traj.reshape([NSTEP,1,3])
    u = mda.Universe.empty(1)
    u.load_new(traj_reshape)
    return u
@pytest.fixture(scope='module')
def step_traj_arr(NSTEP): # constant velocity
    x = np.arange(NSTEP)
    traj = np.vstack([x,x,x]).T
    return traj

def random_walk_3d(NSTEP):
    np.random.seed(1)
    steps = -1 + 2*np.random.randint(0, 2, size=(NSTEP, 3))
    traj = np.cumsum(steps, axis=0)
    traj_reshape = traj.reshape([NSTEP,1,3])
    u = mda.Universe.empty(1)
    u.load_new(traj_reshape)
    return u, traj

def characteristic_poly(n,d): #polynomial that describes unit step trajectory MSD
    x = np.arange(0,n)
    y = d*x*x
    return y

#test some basic size and shape things
def test_selection_works(msd):
    assert_equal(msd._n_particles, 10)

#testing on the  PSF, DCD trajectory
def test_fft_vs_simple_default(msd, msd_fft):
    timeseries_simple = msd.timeseries
    timeseries_fft = msd_fft.timeseries
    assert_almost_equal(timeseries_simple, timeseries_fft, decimal=4)

#check fft and simple give same result per particle
def test_fft_vs_simple_default_per_particle(msd, msd_fft):
    per_particle_simple = msd.msds_by_particle
    per_particle_fft = msd_fft.msds_by_particle
    assert_almost_equal(per_particle_simple, per_particle_fft, decimal=4)

#check fft and simple give same result for each dimensionality
@pytest.mark.parametrize("dim", ['xyz', 'xy', 'xz', 'yz', 'x', 'y', 'z'])
def test_fft_vs_simple_all_dims(u, SELECTION, dim):
    m_simple = MSD(u, SELECTION, msd_type=dim, fft=False)
    m_simple.run()
    timeseries_simple = m_simple.timeseries
    m_fft = MSD(u,SELECTION, msd_type=dim, fft=True)
    m_fft.run()
    timeseries_fft = m_fft.timeseries
    assert_almost_equal(timeseries_simple, timeseries_fft, decimal=4)

#check fft and simple give same result for each particle in each dimension
@pytest.mark.parametrize("dim", ['xyz', 'xy', 'xz', 'yz', 'x', 'y', 'z'])
def test_fft_vs_simple_all_dims_per_particle(u, SELECTION, dim):
    m_simple = MSD(u, SELECTION, msd_type=dim, fft=False)
    m_simple.run()
    per_particle_simple = m_simple.msds_by_particle
    m_fft = MSD(u,SELECTION, msd_type=dim, fft=True)
    m_fft.run()
    per_particle_fft = m_fft.msds_by_particle
    assert_almost_equal(per_particle_simple, per_particle_fft, decimal=4)

#testing the "simple" algorithm on constant velocity trajectory

def test_simple_step_traj_3d(step_traj, NSTEP): # this should fit the polynomial y=3x**2
    m_simple = MSD(step_traj, 'all' , msd_type='xyz', fft=False)
    m_simple.run()
    poly3 = characteristic_poly(NSTEP,3)
    assert_almost_equal(m_simple.timeseries, poly3, decimal=4)

def test_simple_step_traj_2d(step_traj, NSTEP): # this should fit the polynomial y=2x**2
    m_simple = MSD(step_traj, 'all' , msd_type='xy',  fft=False)
    m_simple.run()
    poly2 = characteristic_poly(NSTEP,2)
    assert_almost_equal(m_simple.timeseries, poly2, decimal=4)

def test_simple_step_traj_1d(step_traj, NSTEP): # this should fit the polynomial y=x**2
    m_simple = MSD(step_traj, 'all' , msd_type='x', fft=False)
    m_simple.run()
    poly1 = characteristic_poly(NSTEP,1)
    assert_almost_equal(m_simple.timeseries, poly1,decimal=4)  

#fft based tests require a slight decrease in expected prescision due to roundoff in fft(ifft()) calls
#relative accuracy expected to be around ~1e-12

def test_fft_step_traj_3d(step_traj, NSTEP): # this should fit the polynomial y=3x**2
    m_fft = MSD(step_traj, 'all' , msd_type='xyz', fft=True)
    m_fft.run()
    poly3 = characteristic_poly(NSTEP,3)
    assert_almost_equal(m_fft.timeseries, poly3, decimal=3) # this was relaxed from decimal=4 for numpy=1.13 test

def test_fft_step_traj_2d(step_traj, NSTEP): # this should fit the polynomial y=2x**2
    m_fft = MSD(step_traj, 'all' , msd_type='xy', fft=True)
    m_fft.run()
    poly2 = characteristic_poly(NSTEP,2)
    assert_almost_equal(m_fft.timeseries, poly2, decimal=3) # this was relaxed from decimal=4 for numpy=1.13 test

def test_fft_step_traj_1d(step_traj, NSTEP): # this should fit the polynomial y=x**2
    m_fft = MSD(step_traj, 'all' , msd_type='x', fft=True)
    m_fft.run()
    poly1 = characteristic_poly(NSTEP,1)
    assert_almost_equal(m_fft.timeseries, poly1, decimal=3) # this was relaxed from decimal=4 for numpy=1.13 test


#regress against random_walk test data
def test_random_walk_u_simple(random_walk_u):
    msd_rw = MSD(random_walk_u, 'all', msd_type='xyz', fft=False)
    msd_rw.run()
    norm = np.linalg.norm(msd_rw.timeseries)
    val = 3932.39927487146
    assert_almost_equal(norm, val, decimal=5)

#regress against random_walk test data
def test_random_walk_u_fft(random_walk_u):
    msd_rw = MSD(random_walk_u, 'all', msd_type='xyz', fft=True)
    msd_rw.run()
    norm = np.linalg.norm(msd_rw.timeseries)
    val = 3932.39927487146
    assert_almost_equal(norm, val, decimal=5)
