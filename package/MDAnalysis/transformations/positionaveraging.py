# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

"""\
Trajectory Coordinate Averaging --- :mod:`MDAnalysis.transformations.positionaveraging` 
=======================================================================================

Averages the coordinates of a given trajectory with the N previous frames.
For frames < N, the average of the frames iterated up to that point will be
returned. 

.. autoclass:: PositionAverager

"""
from __future__ import absolute_import

import numpy as np
import warnings



class PositionAverager(object):  
    """
   
    Averages the coordinates of a given timestep so that the coordinates
    of the given atomgroup correspond to the average of the positions 
    of the N previous frames. 
    For frames < N, the average of the frames iterated up to that point will
    be returned.  
    
    Example
    -------

    e.g. Average the coordinates of a given AtomGroup over the course of
    the previous N frames. For N=3, the output will correspond to the
    average of the coordinates over the last 3 frames.
    When check_reset=True, the averager will be reset once the iteration is
    complete, or if the frames iterated are not sequential.

    .. code-block:: python
        
        N=3
        transformation = PositionAverager(N, check_reset=True)
        u.trajectory.add_transformations(transformation)   
        for ts in u.trajectory:
            print(ts.positions)
    
    e.g. When check_reset=False, the average of coordinates from non
    sequential timesteps can also be computed. However, the averager must be
    manually reset before restarting an iteration.
           
    .. code-block:: python
        
        N=3
        transformation = PositionAverager(N, check_reset=False)
        u.trajectory.add_transformations(transformation)
        frames = [0, 7, 1, 6]        
        transformation.resetarrays()
        for ts in u.trajectory[frames]:
            print(ts.positions)
    
    e.g. For frames < N, the average is calculated with the frames iterated up
    to that point and thus will not follow the same behaviour as for the
    frames > N. This behaviour can be followed by comparing current_avg to 
    avg_frames.For N=3, the first 2 frames should be skipped.

    .. code-block:: python
        
        N=3
        transformation = PositionAverager(N, check_reset=True)
        u.trajectory.add_transformations(transformation) 
        for ts in u.trajectory:
            if transformation.current_avg == transformation.avg_frames:
                print(ts.positions)
    
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
    MDAnalysis.coordinates.base.Timestep
    
    
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
    
