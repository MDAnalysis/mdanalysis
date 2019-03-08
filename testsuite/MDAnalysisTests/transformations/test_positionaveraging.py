###

from __future__ import absolute_import

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

import MDAnalysis as md
from MDAnalysis.transformations import PositionAverager
from MDAnalysisTests import datafiles

@pytest.fixture()
def posaveraging_universes():
    '''
    Create the universe objects for the tests.
    '''
    u = md.Universe(datafiles.XTC_multi_frame)
    transformation = PositionAverager(3)
    u.trajectory.add_transformations(transformation)
    return u
   

def test_posavging_fwd(posaveraging_universes):
    '''
    Test if the position averaging function is returning the correct
    values when iterating forwards over the trajectory.
    '''
    ref_matrix_fwd = np.asarray([80., 80., 80.])    
    size = (posaveraging_universes.trajectory.ts.positions.shape[0], 
            posaveraging_universes.trajectory.ts.positions.shape[1], 
            len(posaveraging_universes.trajectory))
    avgd = np.empty(size)    
    for ts in posaveraging_universes.trajectory:
        avgd[...,ts.frame] = ts.positions.copy()
        
    assert_array_almost_equal(ref_matrix_fwd, avgd[1,:,-1], decimal = 5)   

def test_posavging_bwd(posaveraging_universes):
    '''
    Test if the position averaging function is returning the correct
    values when iterating backwards over the trajectory.
    '''
    ref_matrix_bwd = np.asarray([10., 10., 10.])
    size = (posaveraging_universes.trajectory.ts.positions.shape[0], 
            posaveraging_universes.trajectory.ts.positions.shape[1], 
            len(posaveraging_universes.trajectory))
    back_avgd = np.empty(size)
    for ts in posaveraging_universes.trajectory[::-1]:
        back_avgd[...,9-ts.frame] = ts.positions.copy()
    assert_array_almost_equal(ref_matrix_bwd, back_avgd[1,:,-1], decimal = 5)

def test_posavging_reset(posaveraging_universes):
    '''
    Test if the automatic reset is working as intended.
    '''
    size = (posaveraging_universes.trajectory.ts.positions.shape[0], 
            posaveraging_universes.trajectory.ts.positions.shape[1], 
            len(posaveraging_universes.trajectory))
    avgd = np.empty(size)    
    for ts in posaveraging_universes.trajectory:
        avgd[...,ts.frame] = ts.positions.copy()
    after_reset = ts.positions.copy()
    assert_array_almost_equal(avgd[...,0], after_reset, decimal = 5)
    

