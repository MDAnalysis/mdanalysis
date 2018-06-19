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
Translate trajectory --- :mod:`MDAnalysis.transformations.translate`
====================================================================

Translate the coordinates of a given trajectory by a given vector.
This is a wrapped function so usage is as the following example:

    ts = MDAnalysis.transformations.translate(vector)(timestep)
    
"""
from __future__ import absolute_import

import numpy as np

from ..lib.mdamath import triclinic_vectors
from ..core.groups import AtomGroup

def translate(vector):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance by a given vector.
    
    Example
    -------
        ts = MDAnalysis.transformations.translate([1,2,3])(ts)    
    
    Parameters
    ----------
    vector: list
        coordinates of the vector to which the coordinates will be translated
        
    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    """       
    def wrapped(ts):
        if len(vector)>2:
            v = np.float32(vector)
            ts.positions += v
        else:
            raise ValueError("{} vector is too short".format(vector))
        
        return ts
    
    return wrapped

def center_in_box(ag, center='geometry', box=None, pbc=None):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance so that the center of geometry/mass of the given :class:`~MDAnalysis.core.groups.AtomGroup`
    is centered on the unit cell.
    
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
    box: array
        Box dimensions, can be either orthogonal or triclinic information.
        Cell dimensions must be in an identical to format to those returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
        ``[lx, ly, lz, alpha, beta, gamma]``. If ``None``, uses these
        timestep dimensions.
    pbc: bool or None, optional
        If True, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the center of mass or geometry
    
    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    """
    
    pbc_arg = pbc
    if isinstance(ag, AtomGroup):
        if center == "geometry":
            ag_center = ag.center_of_geometry(pbc=pbc_arg)
        elif center == "mass":
            ag_center = ag.center_of_mass(pbc=pbc_arg)
        else:
            raise ValueError('{} is not a valid argument for center'.format(center))
    else:
        raise ValueError('{} is not an AtomGroup object'.format(ag))
    
    def wrapped(ts):
        if box:
            boxcenter = np.sum(triclinic_vectors(box), axis=0) / 2
        else:
            boxcenter = np.sum(ts.triclinic_dimensions, axis=0) / 2
        vector = boxcenter - ag_center
        translate(vector)(ts)
        
        return ts
    
    return wrapped
        
            