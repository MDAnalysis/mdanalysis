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
Trajectory translation --- :mod:`MDAnalysis.transformations.translate`
======================================================================

Translate the coordinates of a given trajectory by a given vector.
The vector can either be user defined, using the function :func:`translate`
or defined by centering an AtomGroup in the unit cell using the function
:func:`center_in_box`
    
"""
from __future__ import absolute_import, division

import numpy as np
from numbers import Number
from six import string_types
from functools import partial

from ..lib.mdamath import triclinic_vectors

def translate(vector):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance by a given vector.
    
    Example
    -------
        ts = MDAnalysis.transformations.translate([1,2,3])(ts)    
    
    Parameters
    ----------
    vector: array-like
        coordinates of the vector to which the coordinates will be translated
        
    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    """
    if len(vector)>2:
        vector = np.float32(vector)
    else:
        raise ValueError("{} vector is too short".format(vector))
    
    def wrapped(ts):
        ts.positions += vector
        
        return ts
    
    return wrapped


def center_in_box(ag, center='geometry', point=None, wrap=False):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance so that the center of geometry/mass of the given :class:`~MDAnalysis.core.groups.AtomGroup`
    is centered on the unit cell. The unit cell dimensions are taken from the input Timestep object.
    If a point is given, the center of the atomgroup will be translated to this point instead. 
    
    Example
    -------
    
    .. code-block:: python
    
        ag = u.residues[1].atoms
        ts = MDAnalysis.transformations.center(ag,center='mass')(ts)
    
    Parameters
    ----------
    ag: AtomGroup
        atom group to be centered on the unit cell.
    center: str, optional
        used to choose the method of centering on the given atom group. Can be 'geometry'
        or 'mass'
    point: array-like, optional
        overrides the unit cell center - the coordinates of the Timestep are translated so
        that the center of mass/geometry of the given AtomGroup is aligned to this position
        instead. Defined as an array of size 3.
    wrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the center of mass or geometry. Default is `False`, no changes
        to the atom coordinates are done before calculating the center of the AtomGroup. 
    
    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    """
    
    pbc_arg = wrap
    if point:
        point = np.asarray(point, np.float32)
        if point.shape != (3, ) and point.shape != (1, 3):
            raise ValueError('{} is not a valid point'.format(point))
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
  
    def wrapped(ts):
        if point is None:
            boxcenter = np.sum(ts.triclinic_dimensions, axis=0) / 2
        else:
            boxcenter = point
    
        ag_center = center_method()

        vector = boxcenter - ag_center
        ts.positions += vector
        
        return ts
    
    return wrapped

    
def center_in_plane(ag, plane, coordinate=0, origin=None, center='geometry', wrap=False):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance so that the center of geometry/mass of the given :class:`~MDAnalysis.core.groups.AtomGroup`
    is centered on a plane. This plane can only be parallel to the x, y or z planes. 
    
    Example
    -------
    
    .. code-block:: python
    
        ag = u.residues[1].atoms
        ts = MDAnalysis.transformations.center_in_plane(ag,'x', 1, center='mass')(ts)
    
    Parameters
    ----------
    ag: AtomGroup
        atom group to be centered on the unit cell.
    plane: str
        used to define the plane on which the given AtomGroup will be centered. Defined as
        a string or an array-like of the plane and the value of its translation. Suported 
        planes are x, y and z.
    coordinate: scalar
        the coordinate on which the plane is defined ( 1 for the plane x=1, for example).
        Default is 0.
    origin: array-like or str
        coordinate on which the axes are centered. Can be an array of three coordinates or
        `center` to shift the origin to the center of the unit cell. Default is `None` 
        which sets [0, 0, 0] as the origin of the axes.
    center: str, optional
        used to choose the method of centering on the given atom group. Can be 'geometry'
        or 'mass'
    wrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the center of mass or geometry. Default is `False`, no changes
        to the atom coordinates are done before calculating the center of the AtomGroup. 
    
    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    """
    
    if isinstance(plane, string_types):
        if plane not in ('x', 'y', 'z'):
            raise ValueError('{} is not a valid plane'.format(plane))
    if not isinstance(coordinate, Number):
        raise ValueError('{} is not a valid coordinate number'.format(coordinate))
    if origin is not None:
        if isinstance(origin, string_types):
            if not origin == 'center':
                raise ValueError('{}is not a valid origin'.format(origin))
        else:
            origin = np.asarray(origin, np.float32)
            if origin.shape != (3, ) and origin.shape != (1, 3):
                raise ValueError('{} is not a valid origin'.format(origin))
    else:
        origin = np.asarray([0, 0, 0], np.float32)

    def wrapped(ts):
        boxcenter = ts.triclinic_dimensions.sum(axis=0) / 2.0
        if isinstance(origin, string_types) and origin == 'center':
            _origin = boxcenter
        else:
            _origin = origin
        if plane == 'x':
            shift = _origin[0]+coordinate
            position = [shift, boxcenter[1], boxcenter[2]]
        if plane == 'y':
            shift = _origin[1]+coordinate
            position = [boxcenter[0], shift, boxcenter[2]]
        if plane == 'z':
            shift = _origin[2]+coordinate
            position = [boxcenter[0], boxcenter[1], shift]
        ts = center_in_box(ag, center=center, point=position, wrap=wrap)(ts)
        
        return ts
    
    return wrapped


def center_in_axis(ag, axis, origin=None, center='geometry', wrap=False):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance so that the center of geometry/mass of the given :class:`~MDAnalysis.core.groups.AtomGroup`
    is centered on an axis. This axis can only be parallel to the x, y or z planes. 
    
    Example
    -------
    
    .. code-block:: python
    
        ag = u.residues[1].atoms
        ts = MDAnalysis.transformations.center_in_axis(ag,'x', [0,0,1], center='mass')(ts)
    
    Parameters
    ----------
    ag: AtomGroup
        atom group to be centered on the unit cell.
    axis: str
        used to define the plane on which the given AtomGroup will be centered. Defined as
        a string or an array-like of the plane and the value of its translation. Suported 
        planes are x, y and z.
    origin: array-like or str
        coordinate on which the axes are centered. Can be an array of three coordinates or
        `center` to shift the origin to the center of the unit cell. Default is `None` 
        which sets [0, 0, 0] as the origin of the axis. 
    center: str, optional
        used to choose the method of centering on the given atom group. Can be 'geometry'
        or 'mass'
    wrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the center of mass or geometry. Default is `False`, no changes
        to the atom coordinates are done before calculating the center of the AtomGroup. 
    
    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    """
    pbc_arg = wrap
    if isinstance(axis, string_types):
        if axis not in ('x', 'y', 'z'):
            raise ValueError('{} is not a valid axis'.format(axis))
    if origin is not None:
        if isinstance(origin, string_types):
            if not origin == 'center':
                raise ValueError('{}is not a valid origin'.format(origin))
        else:
            origin = np.asarray(origin, np.float32)
            if origin.shape != (3, ) and origin.shape != (1, 3):
                raise ValueError('{} is not a valid origin'.format(origin))
    else:
        origin = np.asarray([0, 0, 0], np.float32)
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
    
    def wrapped(ts):
        if isinstance(origin, string_types) and origin == 'center':
            _origin = ts.triclinic_dimensions.sum(axis=0) / 2.0
        else:
            _origin = origin
        ag_center = center_method()
        if axis == 'x':
            center = np.asarray([ag_center[0], _origin[1], _origin[2]], np.float32)
        if axis == 'y':
            center = np.asarray([_origin[0], ag_center[1], _origin[2]], np.float32)
        if axis == 'z':
            center = np.asarray([_origin[0], _origin[1], ag_center[2]], np.float32)
        vector = center - ag_center
        ts.positions += vector
        
        return ts
    
    return wrapped
