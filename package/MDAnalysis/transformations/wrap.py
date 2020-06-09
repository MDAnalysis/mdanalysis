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
Wrap/unwrap transformations --- :mod:`MDAnalysis.transformations.wrap`
======================================================================

Wrap/unwrap the atoms of a given AtomGroup in the unit cell. :func:`wrap` 
translates the atoms back in the unit cell. :func:`unwrap` translates the
atoms of each molecule so that bons don't split over images.

.. autofunction:: wrap

.. autofunction:: unwrap

"""

from __future__ import absolute_import

from ..lib._cutil import make_whole

def wrap(ag, compound='atoms'):
    """
    Shift the contents of a given AtomGroup back into the unit cell. ::
    
       +-----------+          +-----------+
       |           |          |           |
       |         3 | 6        | 6       3 |
       |         ! | !        | !       ! |
       |       1-2-|-5-8  ->  |-5-8   1-2-|
       |         ! | !        | !       ! |
       |         4 | 7        | 7       4 |
       |           |          |           |
       +-----------+          +-----------+
    
    Example
    -------
    
    .. code-block:: python
        
        ag = u.atoms 
        transform = mda.transformations.wrap(ag)
        u.trajectory.add_transformations(transform)
        
    Parameters
    ----------
    
    ag: Atomgroup
        Atomgroup to be wrapped in the unit cell
    compound : {'atoms', 'group', 'residues', 'segments', 'fragments'}, optional
        The group which will be kept together through the shifting process.
    
    Notes
    -----
    When specifying a `compound`, the translation is calculated based on
    each compound. The same translation is applied to all atoms
    within this compound, meaning it will not be broken by the shift.
    This might however mean that not all atoms from the compound are
    inside the unit cell, but rather the center of the compound is.
    
    Returns
    -------
    MDAnalysis.coordinates.base.Timestep
    
    """
    
    def wrapped(ts):
        ag.wrap(compound=compound)
        
        return ts
    
    return wrapped


def unwrap(ag):
    """
    Move all atoms in an AtomGroup so that bonds don't split over images

    Atom positions are modified in place.

    This function is most useful when atoms have been packed into the primary
    unit cell, causing breaks mid molecule, with the molecule then appearing
    on either side of the unit cell. This is problematic for operations
    such as calculating the center of mass of the molecule. ::
    
       +-----------+     +-----------+
       |           |     |           |
       | 6       3 |     |         3 | 6
       | !       ! |     |         ! | !
       |-5-8   1-2-| ->  |       1-2-|-5-8
       | !       ! |     |         ! | !
       | 7       4 |     |         4 | 7
       |           |     |           |
       +-----------+     +-----------+
    
    Example
    -------
    
    .. code-block:: python
        
        ag = u.atoms 
        transform = mda.transformations.unwrap(ag)
        u.trajectory.add_transformations(transform)
    
    Parameters
    ----------
    atomgroup : AtomGroup
        The :class:`MDAnalysis.core.groups.AtomGroup` to work with.
        The positions of this are modified in place.
    
    Returns
    -------
    MDAnalysis.coordinates.base.Timestep
    
    """
    
    try:
        ag.fragments
    except AttributeError:
        raise AttributeError("{} has no fragments".format(ag))
    
    def wrapped(ts):
        for frag in ag.fragments:
            make_whole(frag)
            
        return ts
    
    return wrapped

