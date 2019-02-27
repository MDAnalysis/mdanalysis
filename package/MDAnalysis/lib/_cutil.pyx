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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
#

# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: embedsignature=False
# Warning: Sphinx chokes if embedsignature is True

import numpy as np
cimport numpy as np
from libc.math cimport sqrt, fabs

from MDAnalysis import NoDataError

from libcpp.set cimport set as cset
from libcpp.map cimport map as cmap
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref


__all__ = ['unique_int_1d', 'isrange_int_1d', 'argwhere_int_1d',
           'make_whole', 'find_fragments']


cdef extern from "calc_distances.h":
    ctypedef float coordinate[3]
    void minimum_image(double* x, float* box, float* inverse_box)
    void minimum_image_triclinic(double* dx, float* box)

ctypedef cset[int] intset
ctypedef cmap[int, intset] intmap


def unique_int_1d(np.intp_t[:] values, bint return_counts=False):
    """Find the unique elements of a 1D array of integers.

    This function is optimal on sorted arrays.

    Parameters
    ----------
    values: numpy.ndarray
        1D array of dtype ``numpy.intp`` (or equivalent) in which to find the
        unique values.
    return_counts: bool, optional
        If ``True``, the number of occurrences of each unique value is returned
        as well.

    Returns
    -------
    unique: numpy.ndarray
        A deduplicated copy of `values`.
    counts: numpy.ndarray, optional
        An array of the same length as `unique` containing the number of
        occurrences of each unique value in the original `values` array. Only
        returned if `return_counts` is ``True``.

    Notes
    -----
    The dtype ``numpy.intp`` is usually equivalent to ``numpy.int32`` on a 32
    bit operating system, and, likewise, equivalent to ``numpy.int64`` on a 64
    bit operating system. The exact behavior is compiler-dependent and can be
    checked with ``print(numpy.intp)``.


    See Also
    --------
    :func:`numpy.unique`


    .. versionadded:: 0.19.0
    .. versionchanged:: 0.20.0
       Added optional  `return_counts` parameter and changed dtype from
       ``np.int64`` to ``np.intp`` (corresponds to atom indices).
    """
    cdef np.intp_t[::1] result
    cdef np.intp_t[::1] counts
    cdef np.intp_t n_values

    if return_counts:
        n_values = values.shape[0]
        counts = np.ones(n_values, dtype=np.intp)
        result = _unique_int_1d_with_counts(values, counts)
        n_values = result.shape[0]
        return np.asarray(result), np.asarray(counts[:n_values])

    result = _unique_int_1d(values)
    return np.asarray(result)


cdef inline np.intp_t[::1] _unique_int_1d(np.intp_t[:] values):
    """C/C++-level implementation of :func:`unique_int_1d` with
    ``return_counts=False``.


    .. versionadded:: 0.20.0
    """
    cdef bint is_monotonic = True
    cdef np.intp_t i = 0
    cdef np.intp_t j = 0
    cdef np.intp_t n_values = values.shape[0]
    cdef np.intp_t[::1] result = np.empty(n_values, dtype=np.intp)
    cdef np.intp_t[::1] srt
    if n_values > 0:
        result[0] = values[0]
        if n_values > 1:
            for i in range(1, n_values):
                if values[i] != result[j]:
                    j += 1
                    result[j] = values[i]
                if values[i] < values[i - 1]:
                    is_monotonic = False
            n_values = j + 1
            if is_monotonic:
                result = result[:n_values]
            else:
                # Here could also call the function recursively but that would
                # render inlining impossible
                j = 0
                srt = np.sort(result[:n_values])
                result[0] = srt[0]
                for i in range(1, n_values):
                    if srt[i] != result[j]:
                        j += 1
                        result[j] = srt[i]
                result = result[:j + 1]
    return result


cdef inline np.intp_t[::1] _unique_int_1d_with_counts(np.intp_t[:] values,
                                                      np.intp_t[::1] counts):
    """C/C++-level implementation of :func:`unique_int_1d` with
    ``return_counts=True``.


    .. versionadded:: 0.20.0
    """
    cdef bint is_monotonic = True
    cdef np.intp_t i = 0
    cdef np.intp_t j = 0
    cdef np.intp_t n_values = values.shape[0]
    cdef np.intp_t[::1] result = np.empty(n_values, dtype=np.intp)
    cdef np.intp_t[::1] srt
    if n_values > 0:
        result[0] = values[0]
        if n_values > 1:
            for i in range(1, n_values):
                if values[i] == result[j]:
                    counts[j] += 1
                else:
                    j += 1
                    result[j] = values[i]
                if values[i] < values[i - 1]:
                    is_monotonic = False
                    break
            if is_monotonic:
                result = result[:j + 1]
            else:
                # reset counts, sort the values and try again:
                counts[:j + 1] = 1
                j = 0
                srt = np.sort(values)
                result[0] = srt[0]
                for i in range(1, n_values):
                    if srt[i] == result[j]:
                        counts[j] += 1
                    else:
                        j += 1
                        result[j] = srt[i]
                result = result[:j + 1]
    return result


def isrange_int_1d(np.intp_t[:] values):
    """Checks if an integer array is a contiguous range.

    Checks if ``values[i+1] == values[i] + 1`` holds for all elements of
    `values`.

    Parameters
    ----------
    values: numpy.ndarray
        1D array of dtype ``numpy.intp``.

    Returns
    -------
    bool
        ``True`` if `values` is a contiguous range of numbers, ``False``
        otherwise or if `values` is empty.

    Notes
    -----
    The dtype ``numpy.intp`` is usually equivalent to ``numpy.int32`` on a 32
    bit operating system, and, likewise, equivalent to ``numpy.int64`` on a 64
    bit operating system. The exact behavior is compiler-dependent and can be
    checked with ``print(numpy.intp)``.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t n_values = values.shape[0]
    cdef np.intp_t i
    cdef bint isrange = True
    with nogil:
        if n_values > 1:
            if (values[n_values - 1] - values[0] + 1) == n_values:
                for i in range(n_values - 1):
                    if values[i] + 1 != values[i + 1]:
                        isrange = False
                        break
            else:
                isrange = False
        elif n_values == 0:
            isrange = False
    return isrange


def argwhere_int_1d(np.intp_t[:] arr, np.intp_t value):
    """Find the array indices where elements of `arr` are equal to `value`.

    This function is similar to calling `numpy.argwhere(arr == value).ravel()`
    but is a bit faster since it avoids creating an intermediate boolean array.

    Parameters
    ----------
    arr : numpy.ndarray
        The array to search. Must be one-dimensional and of dtype ``numpy.intp``
        or equivalent.
    value : int
        The value to search for. Must be of a type compatible with
        ``numpy.intp``.

    Returns
    -------
    numpy.ndarray
        A one-dimensional array of dtype ``numpy.intp`` containing all indices
        ``i`` which satisfy the condition ``arr[i] == value``.

    Notes
    -----
    The dtype ``numpy.intp`` is usually equivalent to ``numpy.int32`` on a 32
    bit operating system, and, likewise, equivalent to ``numpy.int64`` on a 64
    bit operating system. The exact behavior is compiler-dependent and can be
    checked with ``print(numpy.intp)``.

    See Also
    --------
    :func:`numpy.argwhere`


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t[::1] result = np.empty(arr.shape[0], dtype=np.intp)
    cdef np.intp_t nargs = _argwhere_int_1d(arr, value, result)
    return np.asarray(result[:nargs])


cdef inline np.intp_t _argwhere_int_1d(np.intp_t[:] arr, np.intp_t value,
                                       np.intp_t[::1] result) nogil:
    """C/C++-level implementation of :func:`argwhere_int_1d`.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t i = 0
    cdef np.intp_t nargs = 0
    cdef np.intp_t n = arr.shape[0]

    for i in range(n):
        if arr[i] == value:
            result[nargs] = i
            nargs += 1

    return nargs


cdef intset difference(intset a, intset b):
    """a.difference(b)

    Returns set of values in a which are not in b
    """
    cdef intset output
    for val in a:
        if b.count(val) != 1:
            output.insert(val)
    return output


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
        must belong in the same molecule or fragment.
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
    cdef np.intp_t i, j, nloops, ref, atom, other, natoms
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
        if not reference_atom in atomgroup:
            raise ValueError("Reference atom not in atomgroup")
        ref = ix_to_rel[reference_atom.ix]

    box = atomgroup.dimensions
    for i in range(3):
        half_box[i] = 0.5 * box[i]
        if box[i] == 0.0:
            raise ValueError("One or more dimensions was zero.  "
                             "You can set dimensions using 'atomgroup.dimensions='")

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
    while <np.intp_t> refpoints.size() < natoms and nloops < natoms:
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

    if <np.intp_t> refpoints.size() < natoms:
        raise ValueError("AtomGroup was not contiguous from bonds, process failed")
    if inplace:
        atomgroup.positions = newpos
    return np.array(newpos)


cdef float _dot(float * a, float * b):
    """Return dot product of two 3d vectors"""
    cdef ssize_t n
    cdef float sum1

    sum1 = 0.0
    for n in range(3):
        sum1 += a[n] * b[n]
    return sum1


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


cpdef np.float64_t _sarrus_det_single(np.float64_t[:, ::1] m):
    """Computes the determinant of a 3x3 matrix."""
    cdef np.float64_t det
    det = m[0, 0] * m[1, 1] * m[2, 2]
    det -= m[0, 0] * m[1, 2] * m[2, 1]
    det += m[0, 1] * m[1, 2] * m[2, 0]
    det -= m[0, 1] * m[1, 0] * m[2, 2]
    det += m[0, 2] * m[1, 0] * m[2, 1]
    det -= m[0, 2] * m[1, 1] * m[2, 0]
    return det


cpdef np.ndarray _sarrus_det_multiple(np.float64_t[:, :, ::1] m):
    """Computes all determinants of an array of 3x3 matrices."""
    cdef np.intp_t n
    cdef np.intp_t i
    cdef np.float64_t[:] det
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
    cdef np.int64_t[:] atoms_view
    cdef np.int32_t[:, :] bonds_view

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
