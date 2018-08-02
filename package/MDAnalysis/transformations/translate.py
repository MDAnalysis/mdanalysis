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

Translate the coordinates of a given trajectory. Coordinates may be 
translated by a vector, using the function :func:`translate` or defined 
by centering an AtomGroup in the unit cell, a plane or axis using the functions
:func:`center_in_box`, :func:`center_in_plane` or :func:`center_in_axis`,
respectively.

.. autofunction:: translate

.. autofunction:: center_in_box

.. autofunction:: center_in_plane

.. autofunction:: center_in_axis
    
"""
from __future__ import absolute_import, division

import numpy as np
from functools import partial

from ..lib.mdamath import triclinic_vectors
from ..lib.util import get_weights
from ..lib._cutil import make_whole

def translate(vector):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance by a given vector.
    
    Example
    -------
    
    Translate the coordinates of the system by the [1, 2, 3] vector:
    
    .. code-block:: python

        transform = mda.transformations.translate([1,2,3])
        u.trajectory.add_transformations(transform)
    
    Parameters
    ----------
    vector: array-like
        coordinates of the vector to which the coordinates will be translated
        
    Returns
    -------
    MDAnalysis.coordinates.base.Timestep
    
    """
    if len(vector)>2:
        vector = np.float32(vector)
    else:
        raise ValueError("{} vector is too short".format(vector))
    
    def wrapped(ts):
        ts.positions += vector
        
        return ts
    
    return wrapped


def center_in_box(ag, weights=None, center_to=None, wrap=False, unwrap=False):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance so that the center of geometry/mass of the given :class:`~MDAnalysis.core.groups.AtomGroup`
    is centered on the unit cell. The unit cell dimensions are taken from the input Timestep object.
    If a point is given, the center of the atomgroup will be translated to this point instead. 
    
    Example
    -------
    
    Translate the center of mass of of the second residue of the universe u to the center of the unit 
    cell:
    
    .. code-block:: python
    
        ag = u.residues[1].atoms
        transform = mda.transformations.center(ag, weights='mass')
        u.trajectory.add_transformations(transform)
    
    Parameters
    ----------
    ag: AtomGroup
        atom group to be centered on the unit cell.
    weights: {"mass", ``None``} or array_like, optional
        define the weights of the atoms when calculating the center of the AtomGroup.
        With ``"mass"`` uses masses as weights; with ``None`` weigh each atom equally.
        If a float array of the same length as `ag` is provided, use each element of
        the `array_like` as a weight for the corresponding atom in `ag`. Default is 
        None.
    center_to: array-like, optional
        overrides the unit cell center - the coordinates of the Timestep are translated so
        that the center of mass/geometry of the given AtomGroup is aligned to this position
        instead. Defined as an array of size 3.
    wrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the weighted center. Default is `False`, no changes to the atom
        coordinates are done before calculating the center of the AtomGroup. 
    unwrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved so as to not break 
        any bonds over periodic boundaries before calculating the weighted center. Default
        is `False`, no changes to the atom coordinates are done before calculating the center
        of the AtomGroup.
    
    Returns
    -------
    MDAnalysis.coordinates.base.Timestep
    
    """
    
    pbc_arg = wrap
    if center_to:
        center_to = np.asarray(center_to, np.float32)
        if center_to.shape != (3, ) and center_to.shape != (1, 3):
            raise ValueError('{} is not a valid point'.format(center_to))
        else: 
            center_to = center_to.reshape(3, )
    try:
        atoms = ag.atoms
    except AttributeError:
        raise ValueError('{} is not an AtomGroup object'.format(ag))
    else:
        try:
            weights = get_weights(atoms, weights=weights)
        except (ValueError, TypeError):
            raise TypeError("weights must be {'mass', None} or an iterable of the "
                        "same size as the atomgroup.")
    if unwrap and wrap:
        raise ValueError("wrap and unwrap can't be both True")
    center_method = partial(atoms.center, weights, pbc=wrap)
    
    def wrapped(ts):
        if center_to is None:
            boxcenter = np.sum(ts.triclinic_dimensions, axis=0) / 2
        else:
            boxcenter = center_to
        if unwrap:
            make_whole(ag)
        ag_center = center_method()

        vector = boxcenter - ag_center
        ts.positions += vector
        
        return ts
    
    return wrapped

    
def center_in_plane(ag, plane, center_to="center", weights=None, wrap=False, unwrap=False):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance so that the center of geometry/mass of the given :class:`~MDAnalysis.core.groups.AtomGroup`
    is centered on a plane. This plane can be the yz, xz or xy planes. 
    
    Example
    -------
    
    Translate the center of mass of the second residue of the universe u to the center of the
    xy plane:

    .. code-block:: python
    
        ag = u.residues[1].atoms
        transform = mda.transformations.center_in_plane(ag, xy, weights='mass')
        u.trajectory.add_transformations(transform)
    
    Parameters
    ----------
    ag: AtomGroup
        atom group to be centered on the unit cell.
    plane: str
        used to define the plane on which the given AtomGroup will be centered. Suported
        planes are yz, xz and xy planes.
    center_to: array-like or str
        coordinates to which the axes are centered. Can be an array of three coordinate 
        values or `center` which centers the AtomGroup coordinates on the center of the 
        unit cell. Default is `center`.
    weights: {"mass", ``None``} or array_like, optional
        define the weights of the atoms when calculating the center of the AtomGroup.
        With ``"mass"`` uses masses as weights; with ``None`` weigh each atom equally.
        If a float array of the same length as `ag` is provided, use each element of
        the `array_like` as a weight for the corresponding atom in `ag`. Default is 
        None.
    wrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the weighted center. Default is `False`, no changes to the atom
        coordinates are done before calculating the center of the AtomGroup. 
    unwrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved so as to not break 
        any bonds over periodic boundaries before calculating the weighted center. Default
        is `False`, no changes to the atom coordinates are done before calculating the center
        of the AtomGroup.
    
    Returns
    -------
    MDAnalysis.coordinates.base.Timestep
    
    """
    try:
        axes = {'yz' : 0, 'xz' : 1, 'xy' : 2}
        plane = axes[plane]
    except (KeyError, TypeError):
        raise ValueError('{} is not a valid plane'.format(plane))
    try:
        if center_to != "center":
            center_to = np.asarray(center_to, np.float32)
            if center_to.shape != (3, ) and center_to.shape != (1, 3):
                raise ValueError('{} is not a valid "center_to"'.format(center_to))
            center_to = center_to.reshape(3, )
    except ValueError:
        raise ValueError('{} is not a valid "center_to"'.format(center_to))
    try:
        atoms = ag.atoms
    except AttributeError:
        raise ValueError('{} is not an AtomGroup object'.format(ag))
    else:
        try:
            weights = get_weights(atoms, weights=weights)
        except (ValueError, TypeError):
            raise TypeError("weights must be {'mass', None} or an iterable of the "
                        "same size as the atomgroup.")
    if unwrap and wrap:
        raise ValueError("wrap and unwrap can't be both True")
    center_method = partial(atoms.center, weights, pbc=wrap)    

    def wrapped(ts):
        boxcenter = ts.triclinic_dimensions.sum(axis=0) / 2.0
        if center_to == 'center':
            _origin = boxcenter
        else:
            _origin = center_to
        if unwrap:
            make_whole(ag)
        position = center_method()
        position[plane] = _origin[plane]
        vector = np.asarray(position, np.float32) - center_method()
        ts.positions += vector
        
        return ts
    
    return wrapped


def center_in_axis(ag, axis, center_to="center", weights=None, wrap=False, unwrap=False):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance so that the center of geometry/mass of the given :class:`~MDAnalysis.core.groups.AtomGroup`
    is centered on an axis. This axis can only be parallel to the x, y or z planes. 
    
    Example
    -------
    
    Translate the center of mass of of the second residue of the universe u to the center of the line 
    parallel to the x axis that passes through the [0, 0, 1] point:
    
    .. code-block:: python
    
        ag = u.residues[1].atoms
        transform = mda.transformations.center_in_axis(ag, 'x', [0,0,1],
        weights='mass')
        u.trajectory.add_transformation(transform)
    
    Parameters
    ----------
    ag: AtomGroup
        atom group to be centered on the unit cell.
    axis: str
        used to define the plane on which the given AtomGroup will be centered. Defined as
        a string or an array-like of the plane and the value of its translation. Suported 
        planes are x, y and z.
    center_to: array-like or str
        coordinate on which the axes are centered. Can be an array of three coordinates or
        `center` to shift the origin to the center of the unit cell. Default is `None` 
        which sets [0, 0, 0] as the origin of the axis. 
    weights: {"mass", ``None``} or array_like, optional
        define the weights of the atoms when calculating the center of the AtomGroup.
        With ``"mass"`` uses masses as weights; with ``None`` weigh each atom equally.
        If a float array of the same length as `ag` is provided, use each element of
        the `array_like` as a weight for the corresponding atom in `ag`. Default is 
        None.
    wrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the weighted center. Default is `False`, no changes to the atom
        coordinates are done before calculating the center of the AtomGroup. 
    unwrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved so as to not break 
        any bonds over periodic boundaries before calculating the weighted center. Default
        is `False`, no changes to the atom coordinates are done before calculating the center
        of the AtomGroup.
        
    Returns
    -------
    MDAnalysis.coordinates.base.Timestep
    
    """
    try:
        axes = {'x' : 0, 'y': 1, 'z' : 2}
        axis = axes[axis]
    except (KeyError, TypeError):
        raise ValueError('{} is not a valid axis'.format(axis))
    try:
        if center_to != "center":
            center_to = np.asarray(center_to, np.float32)
            if center_to.shape != (3, ) and center_to.shape != (1, 3):
                raise ValueError('{} is not a valid "center_to"'.format(center_to))
            center_to = center_to.reshape(3, )
    except ValueError:
        raise ValueError('{} is not a valid "center_to"'.format(center_to))
    try:
        atoms = ag.atoms
    except AttributeError:
        raise ValueError('{} is not an AtomGroup object'.format(ag))
    else:
        try:
            weights = get_weights(atoms, weights=weights)
        except (ValueError, TypeError):
            raise TypeError("weights must be {'mass', None} or an iterable of the "
                            "same size as the atomgroup.")
    if unwrap and wrap:
        raise ValueError("wrap and unwrap can't be both True")
    center_method = partial(ag.atoms.center, weights, pbc=wrap)  
    
    def wrapped(ts):
        if center_to == 'center':
            _origin = ts.triclinic_dimensions.sum(axis=0) / 2.0
        else:
            _origin = center_to
        if unwrap:
            make_whole(ag)
        ag_center = center_method()
        center = _origin
        center[axis] = ag_center[axis]
        vector = center - ag_center
        ts.positions += vector
        
        return ts
    
    return wrapped
