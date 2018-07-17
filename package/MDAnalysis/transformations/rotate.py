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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""\
Rotate trajectory --- :mod:`MDAnalysis.transformations.translate`
=================================================================

Rotates the coordinates by a given angle arround an axis formed by a direction and a
point 
    
"""
from __future__ import absolute_import

import math
import numpy as np
from functools import partial

from ..lib.transformations import rotation_matrix

def rotateby(angle, direction, point=None, center="geometry", wrap=False, ag=None):
    '''
    Rotates the trajectory by a given angle on a given axis. The axis is defined by 
    the user, combining the direction vector and a point. This point can be the center
    of geometry or the center of mass of a user defined AtomGroup, or a list defining custom
    coordinates. 
    e.g. rotate the coordinates by 90 degrees on a x axis centered on a given atom group:
    
    .. code-block:: python
    
        ts = u.trajectory.ts
        angle = 90
        ag = u.atoms()
        rotated = MDAnalysis.transformations.rotate(angle, ag)(ts)
    
    e.g. rotate the coordinates by a custom axis:
    
    .. code-block:: python

        ts = u.trajectory.ts
        angle = 90
        p = [1,2,3]
        d = [0,0,1]
        rotated = MDAnalysis.transformations.rotate(angle, point=point, direction=d)(ts) 
    
    Parameters
    ----------
    angle: float
        rotation angle in degrees
    direction: array-like
        vector that will define the direction of a custom axis of rotation from the
        provided point.
    ag: AtomGroup, optional
        use this to define the center of mass or geometry as the point from where the
        rotation axis will be defined
    center: str, optional
        used to choose the method of centering on the given atom group. Can be 'geometry'
        or 'mass'
    wrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the center of mass or geometry. Default is `False`, no changes
        to the atom coordinates are done before calculating the center of the AtomGroup. 
    point: array-like, optional
        list of the coordinates of the point from where a custom axis of rotation will
        be defined. 

    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    Warning
    -------
    Wrapping/unwrapping the trajectory or performing PBC corrections may not be possible 
    after rotating the trajectory.

    '''
    pbc_arg = wrap
    angle = np.deg2rad(angle)
    if point:
        point = np.asarray(point, np.float32)
        if point.shape != (3, ) and point.shape != (1, 3):
            raise ValueError('{} is not a valid point'.format(point))
    elif ag:
        try:
            if center == 'geometry':
                center_method = partial(ag.center_of_geometry, pbc=pbc_arg)
            elif center == 'mass':
                center_method = partial(ag.center_of_mass, pbc=pbc_arg)
            else:
                raise ValueError('{} is not a valid argument for center'.format(center))
        except AttributeError:
            if center == 'mass':
                raise AttributeError('{} is not an AtomGroup object with masses'.format(ag))
            else:
                raise ValueError('{} is not an AtomGroup object'.format(ag))
    else:
        raise ValueError('A point or an AtomGroup must be specified')
    
    def wrapped(ts):
        if point is None:
            position = center_method()
        else:
            position = point
        rotation = rotation_matrix(angle, direction, position)[:3, :3]
        ts.positions= np.dot(ts.positions, rotation)
        
        return ts
    
    return wrapped
    