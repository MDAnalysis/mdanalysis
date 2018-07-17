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
Fitting transformations --- :mod:`MDAnalysis.transformations.fit`
=================================================================

Translate and rotates the coordinates of a given trajectory to align
a given AtomGroup to a reference structure.
    
"""
from __future__ import absolute_import

import numpy as np
from functools import partial

from ..analysis import align  

def alignto(ag, reference, select='all', weights=None,
        subselection=None, tol_mass=0.1):
    
    """Perform a spatial superposition by minimizing the RMSD.

    Spatially align the group of atoms `ag` to `reference` by
    doing a RMSD fit on `select` atoms.

    A detailed description on the way the fitting is performed can be found in
    :func:`MDAnalysis.analysis.align.alignto`
    
    Example
    -------
    
    ..code-block::python
    
        ag = u.select_atoms("protein")
        ref = mda.Universe("reference.pdb")
        transform = MDAnalysis.transformations.alignto(ag, ref, wheights="mass")
        u.trajectory.add_transformations(transform)

    Parameters
    ----------
    ag : Universe or AtomGroup
       structure to be aligned, a
       :class:`~MDAnalysis.core.groups.AtomGroup` or a whole
       :class:`~MDAnalysis.core.universe.Universe`
    reference : Universe or AtomGroup
       reference structure, a :class:`~MDAnalysis.core.groups.AtomGroup`
       or a whole :class:`~MDAnalysis.core.universe.Universe`
    select : str or dict or tuple, optional
       The selection to operate on; can be one of:

       1. any valid selection string for
          :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` that
          produces identical selections in `ag` and `reference`; or

       2. a dictionary ``{'ag': sel1, 'reference': sel2}`` where *sel1*
          and *sel2* are valid selection strings that are applied to
          `ag` and `reference` respectively (the
          :func:`MDAnalysis.analysis.align.fasta2select` function returns such
          a dictionary based on a ClustalW_ or STAMP_ sequence alignment); or

       3. a tuple ``(sel1, sel2)``

       When using 2. or 3. with *sel1* and *sel2* then these selection strings
       are applied to `ag` and `reference` respectively and should
       generate *groups of equivalent atoms*.  *sel1* and *sel2* can each also
       be a *list of selection strings* to generate a
       :class:`~MDAnalysis.core.groups.AtomGroup` with defined atom order as
       described under :ref:`ordered-selections-label`).
    weights : {"mass", ``None``} or array_like, optional
       choose weights. With ``"mass"`` uses masses as weights; with ``None``
       weigh each atom equally. If a float array of the same length as
       `ag` is provided, use each element of the `array_like` as a
       weight for the corresponding atom in `mobile`.
    tol_mass: float, optional
       Reject match if the atomic masses for matched atoms differ by more than
       `tol_mass`, default [0.1]
    subselection : str or AtomGroup or None (optional)
       Apply the transformation only to this selection.

       ``None`` [default]
           Apply to ``mobile.universe.atoms`` (i.e., all atoms in the
           context of the selection from `ag` such as the rest of a
           protein, ligands and the surrounding water)
       *selection-string*
           Apply to ``mobile.select_atoms(selection-string)``
        :class:`~MDAnalysis.core.groups.AtomGroup`
           Apply to the arbitrary group of atoms

    Returns
    -------
    MDAnalysis.coordinates.base.Timestep
    
    """
    
    def wrapped(ts):
        align.alignto(ag, reference, select, weights, subselection, tol_mass)
        
        return ts

    return wrapped


def fit_translation(ag, reference, plane=None, center_of="geometry"):
    """
    Translates a given AtomGroup so that its center of geometry/mass matches 
    the respective center of the given reference. A plane can be given by the
    user using the option `plane`, and will result in the removal of
    the translation motions of the AtomGroup over that particular plane.
    
    Example
    -------
    
    ..code-block::python
    
        ag = u.select_atoms("protein")
        ref = mda.Universe("reference.pdb")
        transform = MDAnalysis.transformations.fit_translation(ag, ref, center_of="mass")
        u.trajectory.add_transformations(transform)
    
    Parameters
    ----------
    ag : Universe or AtomGroup
       structure to translate, a
       :class:`~MDAnalysis.core.groups.AtomGroup` or a whole 
       :class:`~MDAnalysis.core.universe.Universe`
    reference : Universe or AtomGroup
       reference structure, a :class:`~MDAnalysis.core.groups.AtomGroup` or a whole 
       :class:`~MDAnalysis.core.universe.Universe`
    plane: str, optional
        used to define the plane on which the translations will be removed. Defined as a 
        string of the plane. Suported planes are yz, xz and xy planes.
    center_of: str, optional
        used to choose the method of centering on the given atom group. Can be 'geometry'
        or 'mass'. Default is "geometry".
    
    Returns
    -------
    MDAnalysis.coordinates.base.Timestep
    """
    
    if plane is not None:
        if plane not in ('xy', 'yz', 'xz'):
            raise ValueError('{} is not a valid plane'.format(plane))
    try:
        if center_of == 'geometry':
            ref_center = partial(reference.atoms.center_of_geometry)
            ag_center = partial(ag.atoms.center_of_geometry)
        elif center_of == 'mass':
            ref_center = partial(reference.atoms.center_of_mass)
            ag_center = partial(ag.atoms.center_of_mass)
        else:
            raise ValueError('{} is not a valid argument for "center_of"'.format(center_of))
    except AttributeError:
        if center_of == 'mass':
            raise AttributeError('Either {} or {} is not an Universe/AtomGroup object with masses'.format(ag, reference))
        else:
            raise ValueError('Either {} or {} is not an Universe/AtomGroup object'.format(ag, reference))
    
    reference = np.asarray(ref_center(), np.float32)

    def wrapped(ts):
        center = np.asarray(ag_center(), np.float32)
        if plane == 'yz':
            vector = np.asarray([0, reference[1] - center[1], reference[2] - center[2]], np.float32)
        if plane == 'xz':
            vector = np.asarray([reference[0] - center[0], 0, reference[2] - center[2]], np.float32)
        if plane == 'xy':
            vector = np.asarray([reference[0] - center[0], reference[1] - center[1], 0], np.float32)
        else:
            vector = reference - center
        ts.positions += vector
        
        return ts
    
    return wrapped

   