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

import cython
import numpy as np
cimport numpy as cnp
from libc.math cimport sqrt, fabs

from MDAnalysis import NoDataError

from libcpp.set cimport set as cset
from libcpp.map cimport map as cmap
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cython.operator cimport dereference as deref

cnp.import_array()

__all__ = ['unique_int_1d', 'make_whole', 'find_fragments',
           '_sarrus_det_single', '_sarrus_det_multiple']

cdef extern from "calc_distances.h":
    ctypedef float coordinate[3]
    void minimum_image(double* x, float* box, float* inverse_box)
    void minimum_image_triclinic(double* dx, float* box)

ctypedef cset[int] intset
ctypedef cmap[int, intset] intmap


@cython.boundscheck(False)  # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def unique_int_1d(cnp.intp_t[:] values):
    """Find the unique elements of a 1D array of integers.

    This function is optimal on sorted arrays.

    Parameters
    ----------
    values: numpy.ndarray
        1D array of dtype ``numpy.int64`` in which to find the unique values.

    Returns
    -------
    numpy.ndarray
        A deduplicated copy of `values`.


    .. versionadded:: 0.19.0
    """
    cdef bint is_monotonic = True
    cdef int i = 0
    cdef int j = 0
    cdef int n_values = values.shape[0]
    cdef cnp.intp_t[:] result = np.empty(n_values, dtype=np.intp)

    if n_values == 0:
        return np.array(result)

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


@cython.boundscheck(False)
def _in2d(cnp.intp_t[:, :] arr1, cnp.intp_t[:, :] arr2):
    """Similar to np.in1d except works on 2d arrays

    Parameters
    ----------
    arr1, arr2 : numpy.ndarray, shape (n,2) and (m, 2)
       arrays of integers

    Returns
    -------
    in1d : bool array
      if an element of arr1 was in arr2

    .. versionadded:: 1.1.0
    """
    cdef object out
    cdef ssize_t i
    cdef cset[pair[cnp.intp_t, cnp.intp_t]] hits
    cdef pair[cnp.intp_t, cnp.intp_t] p

    """
    Construct a set from arr2 called hits
    then for each entry in arr1, check if there's a hit in this set

    python would look like:

    hits = {(i, j) for (i, j) in arr2}
    results = np.empty(arr1.shape[0])
    for i, (x, y) in enumerate(arr1):
        results[i] = (x, y) in hits

    return results
    """
    if not arr1.shape[1] == 2 or not arr2.shape[1] == 2:
        raise ValueError("Both arrays must be (n, 2) arrays")

    for i in range(arr2.shape[0]):
        p = pair[cnp.intp_t, cnp.intp_t](arr2[i, 0], arr2[i, 1])
        hits.insert(p)

    out = np.empty(arr1.shape[0], dtype=np.uint8)
    cdef unsigned char[::1] results = out
    for i in range(arr1.shape[0]):
        p = pair[cnp.intp_t, cnp.intp_t](arr1[i, 0], arr1[i, 1])

        if hits.count(p):
            results[i] = True
        else:
            results[i] = False

    return out.astype(bool)


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
def make_whole(atomgroup, reference_atom=None, inplace=True):
    """Move all atoms in a single molecule so that bonds don't split over
    images.

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
        must belong to the same molecule or fragment.
    reference_atom : :class:`~MDAnalysis.core.groups.Atom`
        The atom around which all other atoms will be moved.
        Defaults to atom 0 in the atomgroup.
    inplace : bool, optional
        If ``True``, coordinates are modified in place.

    Returns
    -------
    coords : numpy.ndarray
        The unwrapped atom coordinates.

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

        # MDAnalysis can split AtomGroups into their fragments
        # based on bonding information.
        # Note that this function will only handle a single fragment
        # at a time, necessitating a loop.
        for frag in u.atoms.fragments:
            make_whole(frag)

    Alternatively, to keep a single atom in place as the anchor::

        # This will mean that atomgroup[10] will NOT get moved,
        # and all other atoms will move (if necessary).
        make_whole(atomgroup, reference_atom=atomgroup[10])


    See Also
    --------
    :meth:`MDAnalysis.core.groups.AtomGroup.unwrap`


    .. versionadded:: 0.11.0
    .. versionchanged:: 0.20.0
        Inplace-modification of atom positions is now optional, and positions
        are returned as a numpy array.
    """
    cdef intset refpoints, todo, done
    cdef cnp.intp_t i, j, nloops, ref, atom, other, natoms
    cdef cmap[int, int] ix_to_rel
    cdef intmap bonding
    cdef int[:, :] bonds
    cdef float[:, :] oldpos, newpos
    cdef bint ortho
    cdef float[:] box
    cdef float[:, :] tri_box
    cdef float half_box[3]
    cdef float inverse_box[3]
    cdef double vec[3]
    cdef ssize_t[:] ix_view
    cdef bint is_unwrapped

    # map of global indices to local indices
    ix_view = atomgroup.ix[:]
    natoms = atomgroup.ix.shape[0]

    oldpos = atomgroup.positions

    # Nothing to do for less than 2 atoms
    if natoms < 2:
        return np.array(oldpos)

    for i in range(natoms):
        ix_to_rel[ix_view[i]] = i

    if reference_atom is None:
        ref = 0
    else:
        # Sanity check
        if reference_atom not in atomgroup:
            raise ValueError("Reference atom not in atomgroup")
        ref = ix_to_rel[reference_atom.ix]

    if atomgroup.dimensions is None:
        raise ValueError("No box information available. "
                         "You can set dimensions using 'atomgroup.dimensions='")
    box = atomgroup.dimensions
    for i in range(3):
        half_box[i] = 0.5 * box[i]

    ortho = True
    for i in range(3, 6):
        if box[i] != 90.0:
            ortho = False

    if ortho:
        # If atomgroup is already unwrapped, bail out
        is_unwrapped = True
        for i in range(1, natoms):
            for j in range(3):
                if fabs(oldpos[i, j] - oldpos[0, j]) >= half_box[j]:
                    is_unwrapped = False
                    break
            if not is_unwrapped:
                break
        if is_unwrapped:
            return np.array(oldpos)
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

    newpos = np.zeros((oldpos.shape[0], 3), dtype=np.float32)

    refpoints = intset()  # Who is safe to use as reference point?
    done = intset()  # Who have I already searched around?
    # initially we have one starting atom whose position is in correct image
    refpoints.insert(ref)
    for i in range(3):
        newpos[ref, i] = oldpos[ref, i]

    nloops = 0
    while <cnp.intp_t> refpoints.size() < natoms and nloops < natoms:
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
                    minimum_image_triclinic(&vec[0], &tri_box[0, 0])
                # Then define position of other based on this vector
                for i in range(3):
                    newpos[other, i] = newpos[atom, i] + vec[i]

                # This other atom can now be used as a reference point
                refpoints.insert(other)
            done.insert(atom)

    if <cnp.intp_t> refpoints.size() < natoms:
        raise ValueError("AtomGroup was not contiguous from bonds, process failed")
    if inplace:
        atomgroup.positions = newpos
    return np.array(newpos)


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


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.float64_t _sarrus_det_single(cnp.float64_t[:, ::1] m):
    """Computes the determinant of a 3x3 matrix."""
    cdef cnp.float64_t det
    det = m[0, 0] * m[1, 1] * m[2, 2]
    det -= m[0, 0] * m[1, 2] * m[2, 1]
    det += m[0, 1] * m[1, 2] * m[2, 0]
    det -= m[0, 1] * m[1, 0] * m[2, 2]
    det += m[0, 2] * m[1, 0] * m[2, 1]
    det -= m[0, 2] * m[1, 1] * m[2, 0]
    return det


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray _sarrus_det_multiple(cnp.float64_t[:, :, ::1] m):
    """Computes all determinants of an array of 3x3 matrices."""
    cdef cnp.intp_t n
    cdef cnp.intp_t i
    cdef cnp.float64_t[:] det
    n = m.shape[0]
    det = np.empty(n, dtype=np.float64)
    for i in range(n):
        det[i] = m[i, 0, 0] * m[i, 1, 1] * m[i, 2, 2]
        det[i] -= m[i, 0, 0] * m[i, 1, 2] * m[i, 2, 1]
        det[i] += m[i, 0, 1] * m[i, 1, 2] * m[i, 2, 0]
        det[i] -= m[i, 0, 1] * m[i, 1, 0] * m[i, 2, 2]
        det[i] += m[i, 0, 2] * m[i, 1, 0] * m[i, 2, 1]
        det[i] -= m[i, 0, 2] * m[i, 1, 1] * m[i, 2, 0]
    return np.array(det)


@cython.boundscheck(False)
@cython.wraparound(False)
def find_fragments(atoms, bondlist):
    """Calculate distinct fragments from nodes (atom indices) and edges (pairs
    of atom indices).

    Parameters
    ----------
    atoms : array_like
       1-D Array of atom indices (dtype will be converted to ``numpy.int64``
       internally)
    bonds : array_like
       2-D array of bonds (dtype will be converted to ``numpy.int32``
       internally), where ``bonds[i, 0]`` and ``bonds[i, 1]`` are the
       indices of atoms connected by the ``i``-th bond. Any bonds referring to
       atom indices not in `atoms` will be ignored.

    Returns
    -------
    fragments : list
       List of arrays, each containing the atom indices of a fragment.

    .. versionaddded:: 0.19.0
    """
    cdef intmap bondmap
    cdef intset todo, frag_todo, frag_done
    cdef vector[int] this_frag
    cdef int i, a, b
    cdef cnp.int64_t[:] atoms_view
    cdef cnp.int32_t[:, :] bonds_view

    atoms_view = np.asarray(atoms, dtype=np.int64)
    bonds_view = np.asarray(bondlist, dtype=np.int32)

    # grab record of which atoms I have to process
    # ie set of all nodes
    for i in range(atoms_view.shape[0]):
        todo.insert(atoms_view[i])
    # Process edges into map
    for i in range(bonds_view.shape[0]):
        a = bonds_view[i, 0]
        b = bonds_view[i, 1]
        # only include edges if both are known nodes
        if todo.count(a) and todo.count(b):
            bondmap[a].insert(b)
            bondmap[b].insert(a)

    frags = []

    while not todo.empty():  # While not all nodes have been done
        # Start a new fragment
        frag_todo.clear()
        frag_done.clear()
        this_frag.clear()
        # Grab a start point for next fragment
        frag_todo.insert(deref(todo.begin()))

        # Loop until fragment fully explored
        while not frag_todo.empty():
            # Pop next in this frag todo
            a = deref(frag_todo.begin())
            frag_todo.erase(a)
            if not frag_done.count(a):
                this_frag.push_back(a)
                frag_done.insert(a)
                todo.erase(a)
                for b in bondmap[a]:
                    if not frag_done.count(b):
                        frag_todo.insert(b)

        # Add fragment to output
        frags.append(np.asarray(this_frag))

    return frags
