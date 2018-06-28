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
Wrap/Unwrap trajectory --- :mod:`MDAnalysis.transformations.translate`
====================================================================

###
    
"""
from __future__ import absolute_import

import numpy as np

from ..lib.distances import apply_PBC
from ..lib.mdamath import triclinic_vectors, _is_contiguous

def make_whole(atomgroup, reference_atom=None):
    """Move all atoms in a single molecule so that bonds don't split over images

    Atoms are modified in place.

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


    Parameters
    ----------
    atomgroup : AtomGroup
        The :class:`MDAnalysis.core.groups.AtomGroup` to work with.
        The positions of this are modified in place.  All these atoms
        must belong in the same molecule or fragment.

    Raises
    ------
    NoDataError
        There are no bonds present.
        (See :func:`~MDAnalysis.topology.core.guess_bonds`)

    ValueError
        The algorithm fails to work.  This is usually
        caused by the atomgroup not being a single fragment.
        (ie the molecule can't be traversed by following bonds)


    Example
    -------
    Make fragments whole::

        from MDAnalysis.transformations.pbc import make_whole

        # This algorithm requires bonds, these can be guessed!
        u = mda.Universe(......, guess_bonds=True)

        # MDAnalysis can split molecules into their fragments
        # based on bonding information.
        # Note that this function will only handle a single fragment
        # at a time, necessitating a loop.
        for frag in u.fragments:
          make_whole(frag)


    """
    # check if box is valid
    box = atomgroup.dimensions[:3]
    if all(box == 0.0):
        raise ValueError("Supplied box had zero size")
    # this will fail if the atomgroup has no bonds
    try:
        b = atomgroup.bonds
    except (AttributeError, NoDataError):
        raise NoDataError("The atomgroup is required to have bonds")
    # check if molecule is contiguous
    if not _is_contiguous(atomgroup, atomgroup[0]):
        raise ValueError("atomgroup not contiguous from bonds")
    
    tric = triclinic_vectors(atomgroup.dimensions)
    halfbox = tric.sum(axis=0) / 2.0
    atgp = atomgroup
    diffs = None
    count = 0
    idx = None
    prev ={}
    while True:
        count +=1
        if diffs is not None:
            if idx is not None:
                prev[count] = idx
            idx = np.unique(np.where(diffs)[0])
            atgp = atgp.bonds.atom2[idx]
            # the following check is for cases when an atom in the atom2 group
            # is connected to two atoms in atom1. Since where changing the positions
            # of the atom2 group, suchs cases lead to the while loop being stuck:
            # this atom is moved back and forth between the two atoms it is conected
            # to. Since were tracking the atoms that need changing, if this group stays
            # the same over an iteration, we change atom1.positions instead.
            if prev.has_key(count) and np.array_equal(prev[count], idx):
                atgp.bonds.atom1.positions = atgp.bonds.atom2.positions - newvecs 
        vecs = atgp.bonds.atom2.positions - atgp.bonds.atom1.positions
        vecs += halfbox
        newvecs = apply_PBC(vecs, atgp.dimensions)
        diffs = np.any(~np.isclose(vecs, newvecs, 1e-6), axis=1)
        if not diffs.sum():
            break
        newvecs -= halfbox
        atgp.bonds.atom2.positions = atgp.bonds.atom1.positions + newvecs
        # for debugging purposes. Also makes a nice animation
        #atomgroup.write(filename='updated_%d.pdb'%count)


