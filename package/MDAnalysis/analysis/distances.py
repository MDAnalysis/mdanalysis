# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
#


"""
Distance analysis --- :mod:`MDAnalysis.analysis.distances`
==========================================================

This module provides functions to rapidly compute distances between
atoms or groups of atoms.

:func:`dist` and :func:`between` can take atom groups that do not even
have to be from the same :class:`~MDAnalysis.core.universe.Universe`.

See Also
--------
:mod:`MDAnalysis.lib.distances`
"""

__all__ = [
    "distance_array",
    "self_distance_array",
    "contact_matrix",
    "dist",
    "between",
]

import numpy as np
import scipy.sparse

from MDAnalysis.lib.distances import (
    capped_distance,
    self_distance_array,
    distance_array,  # legacy reasons
)
from MDAnalysis.lib.c_distances import (
    contact_matrix_no_pbc,
    contact_matrix_pbc,
)
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from MDAnalysis.lib.distances import calc_bonds


import warnings
import logging

logger = logging.getLogger("MDAnalysis.analysis.distances")


def contact_matrix(coord, cutoff=15.0, returntype="numpy", box=None):
    """Calculates a matrix of contacts.

    There is a fast, high-memory-usage version for small systems
    (*returntype* = 'numpy'), and a slower, low-memory-usage version for
    larger systems (*returntype* = 'sparse').

    If *box* dimensions are passed then periodic boundary conditions
    are applied.

    Parameters
    ---------
    coord : array
       Array of coordinates of shape ``(N, 3)`` and dtype float32.
    cutoff : float, optional, default 15
       Particles within `cutoff` are considered to form a contact.
    returntype : string, optional, default "numpy"
       Select how the contact matrix is returned.
       * ``"numpy"``: return as an ``(N. N)`` :class:`numpy.ndarray`
       * ``"sparse"``: return as a :class:`scipy.sparse.lil_matrix`
    box : array-like or ``None``, optional, default ``None``
       Simulation cell dimensions in the form of
       :attr:`MDAnalysis.trajectory.timestep.Timestep.dimensions` when
       periodic boundary conditions should be taken into account for
       the calculation of contacts.

    Returns
    -------
    array or sparse matrix
       The contact matrix is returned in a format determined by the `returntype`
       keyword.

    See Also
    --------
    :mod:`MDAnalysis.analysis.contacts` for native contact analysis


    .. versionchanged:: 0.11.0
       Keyword *suppress_progmet* and *progress_meter_freq* were removed.
    """

    if returntype == "numpy":
        adj = np.full((len(coord), len(coord)), False, dtype=bool)
        pairs = capped_distance(
            coord, coord, max_cutoff=cutoff, box=box, return_distances=False
        )

        idx, idy = np.transpose(pairs)
        adj[idx, idy] = True

        return adj
    elif returntype == "sparse":
        # Initialize square List of Lists matrix of dimensions equal to number
        # of coordinates passed
        sparse_contacts = scipy.sparse.lil_matrix(
            (len(coord), len(coord)), dtype="bool"
        )
        if box is not None:
            # with PBC
            contact_matrix_pbc(coord, sparse_contacts, box, cutoff)
        else:
            # without PBC
            contact_matrix_no_pbc(coord, sparse_contacts, cutoff)
        return sparse_contacts


def dist(A, B, offset=0, box=None):
    """Return distance between atoms in two atom groups.

    The distance is calculated atom-wise. The residue ids are also
    returned because a typical use case is to look at CA distances
    before and after an alignment. Using the `offset` keyword one can
    also add a constant offset to the resids which facilitates
    comparison with PDB numbering.

    Arguments
    ---------
    A, B : AtomGroup
       :class:`~MDAnalysis.core.groups.AtomGroup` with the
       same number of atoms
    offset : integer or tuple, optional, default 0
       An integer `offset` is added to *resids_A* and *resids_B* (see
       below) in order to produce PDB numbers.

       If `offset` is :class:`tuple` then ``offset[0]`` is added to
       *resids_A* and ``offset[1]`` to *resids_B*. Note that one can
       actually supply numpy arrays of the same length as the atom
       group so that an individual offset is added to each resid.

    Returns
    -------
    resids_A : array
        residue ids of the `A` group (possibly changed with `offset`)
    resids_B : array
       residue ids of the `B` group (possibly changed with `offset`)
    distances : array
       distances between the atoms
    """

    if A.atoms.n_atoms != B.atoms.n_atoms:
        raise ValueError(
            "AtomGroups A and B do not have the same number of atoms"
        )
    try:
        off_A, off_B = offset
    except (TypeError, ValueError):
        off_A = off_B = int(offset)
    residues_A = np.array(A.resids) + off_A
    residues_B = np.array(B.resids) + off_B

    d = calc_bonds(A.positions, B.positions, box)
    return np.array([residues_A, residues_B, d])


def between(group, A, B, distance):
    """Return sub group of `group` that is within `distance` of both `A` and `B`

    This function is not aware of periodic boundary conditions.

    Can be used to find bridging waters or molecules in an interface.

    Similar to "*group* and (AROUND *A* *distance* and AROUND *B* *distance*)".

    Parameters
    ----------
    group : AtomGroup
        Find members of `group` that are between `A` and `B`
    A : AtomGroup
    B : AtomGroup
        `A` and `B` are the groups of atoms between which atoms in
        `group` are searched for.  The function works is more
        efficient if `group` is bigger than either `A` or `B`.
    distance : float
        maximum distance for an atom to be counted as in the vicinity of
        `A` or `B`

    Returns
    -------
    AtomGroup
        :class:`~MDAnalysis.core.groups.AtomGroup` of atoms that
        fulfill the criterion


    .. versionadded: 0.7.5

    """
    ns_group = AtomNeighborSearch(group)
    resA = ns_group.search(A, distance)
    resB = ns_group.search(B, distance)
    return resB.intersection(resA)
