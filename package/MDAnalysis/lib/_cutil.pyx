# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
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
#

import cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt

from MDAnalysis import NoDataError

from libcpp.set cimport set as cset
from libcpp.map cimport map as cmap

__all__ = ['unique_int_1d', 'make_whole']

cdef extern from "calc_distances.h":
    ctypedef float coordinate[3]
    void minimum_image(double *x, float *box, float *inverse_box)
    void minimum_image_triclinic(double *dx, coordinate *box)

ctypedef cset[int] intset
ctypedef cmap[int, intset] intmap


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def unique_int_1d(np.int64_t[:] values):
    """
    Find the unique elements of a 1D array of integers.

    This function is optimal on sorted arrays.

    Parameters
    ----------
    values: np.ndarray of type int64
        1D array of int in which to find the unique values.

    Returns
    -------
    np.ndarray

    .. versionadded:: 0.19.0
    """
    cdef bint is_monotonic = True
    cdef int i = 0
    cdef int j = 0
    cdef int n_values = values.shape[0]
    cdef np.int64_t[:] result = np.empty(n_values, dtype=np.int64)

    if n_values == 0:
        return result

    result[0] = values[0]
    for i in range(1, n_values):
        if values[i] != result[j]:
            j += 1
            result[j] = values[i]
        if values[i] < values[i - 1]:
            is_monotonic = False
    result = result[:j + 1]
    if not is_monotonic:
        result = unique_int_1d(np.sort(result))

    return np.array(result)


cdef intset difference(intset a, intset b):
    """a.difference(b)

    Returns set of values in a which are not in b
    """
    cdef intset output
    for val in a:
        if b.count(val) != 1:
            output.insert(val)
    return output


@cython.boundscheck(False)
@cython.wraparound(False)
def make_whole(atomgroup, reference_atom=None):
    """Move all atoms in a single molecule so that bonds don't split over images

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


    Parameters
    ----------
    atomgroup : AtomGroup
        The :class:`MDAnalysis.core.groups.AtomGroup` to work with.
        The positions of this are modified in place.  All these atoms
        must belong in the same molecule or fragment.
    reference_atom : :class:`~MDAnalysis.core.groups.Atom`
        The atom around which all other atoms will be moved.
        Defaults to atom 0 in the atomgroup.

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

        from MDAnalysis.lib.mdamath import make_whole

        # This algorithm requires bonds, these can be guessed!
        u = mda.Universe(......, guess_bonds=True)

        # MDAnalysis can split molecules into their fragments
        # based on bonding information.
        # Note that this function will only handle a single fragment
        # at a time, necessitating a loop.
        for frag in u.fragments:
          make_whole(frag)

    Alternatively, to keep a single atom in place as the anchor::

        # This will mean that atomgroup[10] will NOT get moved,
        # and all other atoms will move (if necessary).
        make_whole(atomgroup, reference_atom=atomgroup[10])


    .. versionadded:: 0.11.0
    """
    cdef intset refpoints, todo, done
    cdef int i, nloops, ref, atom, other, natoms
    cdef cmap[int, int] ix_to_rel
    cdef intmap bonding
    cdef int[:, :] bonds
    cdef float[:, :] oldpos, newpos
    cdef bint ortho
    cdef float[:] box
    cdef float tri_box[3][3]
    cdef float inverse_box[3]
    cdef double vec[3]
    cdef ssize_t[:] ix_view

    # map of global indices to local indices
    ix_view = atomgroup.ix[:]
    natoms = atomgroup.ix.shape[0]
    for i in range(natoms):
        ix_to_rel[ix_view[i]] = i

    if reference_atom is None:
        ref = 0
    else:
        # Sanity check
        if not reference_atom in atomgroup:
            raise ValueError("Reference atom not in atomgroup")
        ref = ix_to_rel[reference_atom.ix]

    box = atomgroup.dimensions

    for i in range(3):
        if box[i] == 0.0:
            raise ValueError("One or more dimensions was zero.  "
                             "You can set dimensions using 'atomgroup.dimensions='")

    ortho = True
    for i in range(3, 6):
        if box[i] != 90.0:
            ortho = False

    if ortho:
        for i in range(3):
            inverse_box[i] = 1.0 / box[i]
    else:
        from .mdamath import triclinic_vectors
        tri_box = triclinic_vectors(box)

    # C++ dict of bonds
    try:
        bonds = atomgroup.bonds.to_indices()
    except (AttributeError, NoDataError):
        raise NoDataError("The atomgroup is required to have bonds")
    for i in range(bonds.shape[0]):
        atom = bonds[i, 0]
        other = bonds[i, 1]
        # only add bonds if both atoms are in atoms set
        if ix_to_rel.count(atom) and ix_to_rel.count(other):
            atom = ix_to_rel[atom]
            other = ix_to_rel[other]

            bonding[atom].insert(other)
            bonding[other].insert(atom)

    oldpos = atomgroup.positions
    newpos = np.zeros((oldpos.shape[0], 3), dtype=np.float32)

    refpoints = intset()  # Who is safe to use as reference point?
    done = intset()  # Who have I already searched around?
    # initially we have one starting atom whose position is in correct image
    refpoints.insert(ref)
    for i in range(3):
        newpos[ref, i] = oldpos[ref, i]

    nloops = 0
    while refpoints.size() < natoms and nloops < natoms:
        # count iterations to prevent infinite loop here
        nloops += 1

        # We want to iterate over atoms that are good to use as reference
        # points, but haven't been searched yet.
        todo = difference(refpoints, done)
        for atom in todo:
            for other in bonding[atom]:
                # If other is already a refpoint, leave alone
                if refpoints.count(other):
                    continue
                # Draw vector from atom to other
                for i in range(3):
                    vec[i] = oldpos[other, i] - newpos[atom, i]
                # Apply periodic boundary conditions to this vector
                if ortho:
                    minimum_image(&vec[0], &box[0], &inverse_box[0])
                else:
                    minimum_image_triclinic(&vec[0], <coordinate*>&tri_box[0])
                # Then define position of other based on this vector
                for i in range(3):
                    newpos[other, i] = newpos[atom, i] + vec[i]

                # This other atom can now be used as a reference point
                refpoints.insert(other)
            done.insert(atom)

    if refpoints.size() < natoms:
        raise ValueError("AtomGroup was not contiguous from bonds, process failed")
    else:
        atomgroup.positions = newpos


@cython.boundscheck(False)
@cython.wraparound(False)
cdef float _dot(float * a, float * b):
    """Return dot product of two 3d vectors"""
    cdef ssize_t n
    cdef float sum1

    sum1 = 0.0
    for n in range(3):
        sum1 += a[n] * b[n]
    return sum1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _cross(float * a, float * b, float * result):
    """
    Calculates the cross product between 3d vectors

    Note
    ----
    Modifies the result array
    """

    result[0] = a[1]*b[2] - a[2]*b[1]
    result[1] = - a[0]*b[2] + a[2]*b[0]
    result[2] = a[0]*b[1] - a[1]*b[0]

cdef float _norm(float * a):
    """
    Calculates the magnitude of the vector
    """
    cdef float result
    cdef ssize_t n
    result = 0.0
    for n in range(3):
        result += a[n]*a[n]
    return sqrt(result)