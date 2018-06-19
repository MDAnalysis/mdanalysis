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

from ..lib.transformations import rotation_matrix
from ..core.groups import AtomGroup

def rotateby(angle, direction, center="geometry", pbc=None, ag=None, position=None):
    '''
    Rotates the trajectory by a given angle on a given axis. The axis is defined by 
    the user, combining the direction vector and a position. This position can be the center
    of geometry or the center of mass of a user defined AtomGroup, or a list defining custom
    coordinates. 
    e.g. rotate the coordinates by pi/2 on a x axis centered on a given atom group:
    
    .. code-block:: python
    
        ts = u.trajectory.ts
        angle = math.pi/2
        ag = u.atoms()
        rotated = MDAnalysis.transformations.rotate(angle, ag)(ts)
    
    e.g. rotate the coordinates by a custom axis:
    
    .. code-block:: python

        ts = u.trajectory.ts
        angle = math.pi/2
        p = [1,2,3]
        d = [0,0,1]
        rotated = MDAnalysis.transformations.rotate(angle, point=point, direction=d)(ts) 
    
    Parameters
    ----------
    angle: float
        rotation angle in radians
    direction: list
        vector that will define the direction of a custom axis of rotation from the
        provided point.
    ag: AtomGroup, optional
        use this to define the center of mass or geometry as the point from where the
        rotation axis will be defined
    center: str, optional
        used to choose the method of centering on the given atom group. Can be 'geometry'
        or 'mass'
    pbc: bool or None, optional
        If True, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the center
    position: list, optional
        list of the coordinates of the point from where a custom axis of rotation will
        be defined. 

    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    '''
    pbc_arg = pbc
    if position and len(position)>2:
        position = position
    elif isinstance(ag, AtomGroup):
        if center == "geometry":
            position = ag.center_of_geometry(pbc=pbc_arg)
        elif center == "mass":
            position = ag.center_of_mass(pbc=pbc_arg)
        else:
            raise ValueError('{} is not a valid argument for center'.format(center))
    else:
        raise ValueError('A position or an AtomGroup must be specified') 
    
    def wrapped(ts):
        rotation = rotation_matrix(angle, direction, position)[:3, :3]
        ts.positions= np.dot(ts.positions, rotation)
        
        return ts
    
    return wrapped
    