# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4


"""
Trajectory Coordinate Averaging --- 
=================================

Averages the coordinates of a given trajectory with the N previous frames.
For frames < N, the average of the frames iterated up to that point will be
returned. 

.. autoclass:: PositionAverager

"""
from __future__ import absolute_import, division

import numpy as np
import MDAnalysis as md
import warnings



class PositionAverager:  
    """
    
    Averages the coordinates of a given timestep with the N previous frames. 
    For frames < N, the average of the frames iterated up to that point will be
    returned.  
    
    Example
    -------
    
    u.trajectory.add_transformations(PositionAverager(N, check_reset=True))    
    
    
    Parameters
    ----------
    avg_frames: int
    Determines the number of frames to be used for the position averaging. 
    
    
    check_reset: bool, optional
    If 'True', position averaging will be reset, and a warning raised,
    if the trajectory iteration direction changes. If 'False', position
    averaging will not reset, regardless of the iteration.
    

    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    
    """
    
    def __init__(self, avg_frames, check_reset=True):
        self.avg_frames = avg_frames
        self.check_reset = check_reset
        self.current_avg = 0
        self.resetarrays()
    
    def resetarrays(self):
        self.idx_array = np.empty(self.avg_frames)
        self.idx_array[:] = np.nan
    
    def rollidx(self,ts):
        self.idx_array = np.roll(self.idx_array, 1)  
        self.idx_array[0] = ts.frame
        
    def rollposx(self,ts):        
        try:
            self.coord_array.size
        except AttributeError:
            size = (ts.positions.shape[0], ts.positions.shape[1],
                    self.avg_frames)
            self.coord_array = np.empty(size)
         
        self.coord_array = np.roll(self.coord_array, 1, axis=2)
        self.coord_array[...,0] = ts.positions.copy()
        
    
    def __call__(self, ts):
        self.rollidx(ts)
        test = ~np.isnan(self.idx_array)
        self.current_avg = sum(test)
        if self.current_avg == 1:
            return ts
          
        if self.check_reset:
            sign = np.sign(np.diff(self.idx_array[test]))
            
            if not (np.all(sign == 1) or np.all(sign==-1)):
                warnings.warn('Cannot average position for non sequential'
                              'iterations. Averager will be reset.',
                              Warning)
                self.resetarrays()
                return self(ts)
        
        self.rollposx(ts)
        ts.positions = np.mean(self.coord_array[...,test], axis=2)
            
        return ts
    
