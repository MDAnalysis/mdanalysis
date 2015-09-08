# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2014 Naveen Michaud-Agrawal,
# Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see AUTHORS for the full list)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Distance analysis --- :mod:`MDAnalysis.analysis.distances`
==========================================================

This module provides functions to rapidly compute distances between
atoms or groups of atoms.

:func:`dist` and :func:`between` can take atom groups that do not even
have to be from the same :class:`~MDAnalysis.core.AtomGroup.Universe`.

.. SeeAlso:: :mod:`MDAnalysis.lib.distances` and :mod:`MDAnalysis.lib.parallel.distances`
"""

__all__ = ['distance_array', 'self_distance_array', 'contact_matrix', 'dist']

import numpy as np
from scipy import sparse

from MDAnalysis.lib.distances import distance_array, self_distance_array
from MDAnalysis.lib._distances import contact_matrix_no_pbc, contact_matrix_pbc

import logging

logger = logging.getLogger("MDAnalysis.analysis.distances")


def contact_matrix(coord, cutoff=15.0, returntype="numpy", box=None):
    '''Calculates a matrix of contacts within a numpy array of type float32.

    There is a fast, high-memory-usage version for small systems
    (*returntype* = 'numpy'), and a slower, low-memory-usage version for
    larger systems (*returntype* = 'sparse').

    If *box* dimensions are passed (``box = [Lx, Ly, Lz]``), then
    periodic boundary conditions are applied.  Only orthorhombic boxes
    are currently supported.

    .. versionchanged:: 0.11.0
       Keyword *suppress_progmet* and *progress_meter_freq* were removed.
    '''
    if returntype == "numpy":
        adj = (distance_array(coord, coord, box=box) < cutoff)
        return adj
    elif returntype == "sparse":
        # Initialize square List of Lists matrix of dimensions equal to number of coordinates passed
        sparse_contacts = sparse.lil_matrix((len(coord), len(coord)), dtype='bool')
        if box is not None:
            # if PBC
            contact_matrix_pbc(coord, sparse_contacts, box, cutoff)
        else:
            # if no PBC
            contact_matrix_no_pbc(coord, sparse_contacts, cutoff)
        return sparse_contacts


def dist(A, B, offset=0):
    """Return distance between atoms in two atom groups.

    The distance is calculated atom-wise. The residue ids are also
    returned because a typical use case is to look at CA distances
    before and after an alignment. Using the *offset* keyword one can
    also add a constant offset to the resids which facilitates
    comparison with PDB numbering.

    :Arguments:
       *A*, *B*
          :class:`~MDAnalysis.core.AtomGroup.AtomGroup` with the
          same number of atoms

    :Keywords:
       *offset* : integer
          The *offset* is added to *resids_A* and *resids_B* (see
          below) in order to produce PDB numbers. The default is 0.

       *offset* : tuple
          *offset[0]* is added to *resids_A* and *offset[1]* to
          *resids_B*. Note that one can actually supply numpy arrays
          of the same length as the atom group so that an individual
          offset is added to each resid.

    :Returns: NumPy `array([resids_A, resids_B, distances])`

    """
    if A.atoms.n_atoms != B.atoms.n_atoms:
        raise ValueError("AtomGroups A and B do not have the same number of atoms")
    try:
        off_A, off_B = offset
    except (TypeError, ValueError):
        off_A = off_B = int(offset)
    residues_A = np.array(A.resids) + off_A
    residues_B = np.array(B.resids) + off_B
    r = A.coordinates() - B.coordinates()
    d = np.sqrt(np.sum(r * r, axis=1))
    return np.array([residues_A, residues_B, d])


def between(group, A, B, distance):
    """Return sub group of *group* that is within *distance* of both *A* and *B*.

    *group*, *A*, and *B* must be
    :class:`~MDAnalysis.core.AtomGroup.AtomGroup` instances.  Works best
    if *group* is bigger than either *A* or *B*. This function is not
    aware of periodic boundary conditions.

    Can be used to find bridging waters or molecules in an interface.

    Similar to "*group* and (AROUND *A* *distance* and AROUND *B* *distance*)".

    .. SeeAlso:: Makes use of :mod:`MDAnalysis.lib.NeighborSearch`.

    .. versionadded: 0.7.5
    """
    from MDAnalysis.core.AtomGroup import AtomGroup

    ns_group = AtomNeighborSearch(group)
    resA = set(ns_group.search_list(A, distance))
    resB = set(ns_group.search_list(B, distance))
    return AtomGroup(resB.intersection(resA))
