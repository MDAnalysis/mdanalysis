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
import numpy as np
import warnings

from .base import TransformationBase


class PositionAverager(TransformationBase):
    """
    Averages the coordinates of a given timestep so that the coordinates
    of the AtomGroup correspond to the average positions of the N previous
    frames. 
    For frames < N, the average of the frames iterated up to that point will
    be returned.  

    Example
    -------

    Average the coordinates of a given AtomGroup over the course of
    the previous N frames. For ``N=3``, the output will correspond to the
    average of the coordinates over the last 3 frames.
    When ``check_reset=True``, the averager will be reset once the iteration is
    complete, or if the frames iterated are not sequential.

    .. code-block:: python

        N=3
        transformation = PositionAverager(N, check_reset=True)
        u.trajectory.add_transformations(transformation)   
        for ts in u.trajectory:
            print(ts.positions)

    In this case, ``ts.positions`` will return the average coordinates of the
    last N iterated frames.


    When ``check_reset=False``, the average of coordinates from non
    sequential timesteps can also be computed. However, the averager must be
    manually reset before restarting an iteration. In this case,
    ``ts.positions`` will return the average coordinates of the last N
    iterated frames, despite them not being sequential
    (``frames = [0, 7, 1, 6]``). 

    .. code-block:: python
        
        N=3
        transformation = PositionAverager(N, check_reset=False)
        u.trajectory.add_transformations(transformation)
        frames = [0, 7, 1, 6]        
        transformation.resetarrays()
        for ts in u.trajectory[frames]:
            print(ts.positions)
    
    If ``check_reset=True``, the ``PositionAverager`` would have automatically
    reset after detecting a non sequential iteration (i.e. when iterating from
    frame 7 to frame 1 or when resetting the iterator from frame 6 back to
    frame 0).


    For frames < N, the average is calculated with the frames iterated up
    to that point and thus will not follow the same behaviour as for
    frames > N. This can be followed by comparing the number of frames being
    used to compute the current averaged frame (``current_avg``) to the one
    requested when calling ``PositionAverager`` (``avg_frames``) which in
    these examples corresponds to ``N=3``.

    .. code-block:: python
        
        N=3
        transformation = PositionAverager(N, check_reset=True)
        u.trajectory.add_transformations(transformation) 
        for ts in u.trajectory:
            if transformation.current_avg == transformation.avg_frames:
                print(ts.positions)
    
    In the case of ``N=3``, as the average is calculated with the frames
    iterated up to the current iteration, the first frame returned will
    not be averaged. During the first iteration no other frames are stored in
    memory thus far and, consequently, ``transformation.current_avg = 1``.
    The second frame iterated will return the average of frame 1 and frame 2,
    with ``transformation.current_avg = 2``. Only during the third and
    following iterations will ``ts.positions`` start returning the average of
    the last 3 frames and thus ``transformation.current_avg = 3``
    These initial frames are typically not desired during analysis, but one can
    easily avoid them, as seen in the previous example with 
    ``if transformation.current_avg == transformation.avg_frames:`` or by
    simply removing the first ``avg_frames-1`` frames from the analysis. 



    Parameters
    ----------
    avg_frames: int
        Determines the number of frames to be used for the position averaging. 
    check_reset: bool, optional
        If ``True``, position averaging will be reset and a warning raised
        when the trajectory iteration direction changes. If ``False``, position
        averaging will not reset, regardless of the iteration.
    

    Returns
    -------
    MDAnalysis.coordinates.base.Timestep


    .. versionchanged:: 2.0.0
       The transformation was changed to inherit from the base class for
       limiting threads and checking if it can be used in parallel analysis.
    """

    def __init__(self, avg_frames, check_reset=True,
                 max_threads=None,
                 parallelizable=False):
        super().__init__(max_threads=max_threads,
                         parallelizable=parallelizable)
        self.avg_frames = avg_frames
        self.check_reset = check_reset
        self.current_avg = 0
        self.resetarrays()
        self.current_frame = 0

    def resetarrays(self):
        self.idx_array = np.empty(self.avg_frames)
        self.idx_array[:] = np.nan

    def rollidx(self, ts):
        self.idx_array = np.roll(self.idx_array, 1)
        self.idx_array[0] = ts.frame

    def rollposx(self, ts):
        try:
            self.coord_array.size
        except AttributeError:
            size = (ts.positions.shape[0], ts.positions.shape[1],
                    self.avg_frames)
            self.coord_array = np.empty(size)

        self.coord_array = np.roll(self.coord_array, 1, axis=2)
        self.coord_array[..., 0] = ts.positions.copy()

    def _transform(self, ts):
        #  calling the same timestep will not add new data to coord_array
        #  This can prevent from getting different values when
        #  call `u.trajectory[i]` multiple times.
        if (ts.frame == self.current_frame
                and hasattr(self, 'coord_array')
                and not np.isnan(self.idx_array).all()):
            test = ~np.isnan(self.idx_array)
            ts.positions = np.mean(self.coord_array[..., test], axis=2)
            return ts
        else:
            self.current_frame = ts.frame

        self.rollidx(ts)
        test = ~np.isnan(self.idx_array)
        self.current_avg = sum(test)
        if self.check_reset:
            sign = np.sign(np.diff(self.idx_array[test]))
            if not (np.all(sign == 1) or np.all(sign == -1)):
                warnings.warn('Cannot average position for non sequential'
                              'iterations. Averager will be reset.',
                              Warning)
                self.resetarrays()
                return self(ts)

        self.rollposx(ts)
        ts.positions = np.mean(self.coord_array[..., test], axis=2)

        return ts
