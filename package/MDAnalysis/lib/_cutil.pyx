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

from __future__ import division

import cython
cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt, fabs, INFINITY, NAN

from MDAnalysis import NoDataError

from libcpp.set cimport set as cset
from libcpp.map cimport map as cmap
from libcpp.list cimport list as clist
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref


__all__ = ['coords_add_vector', 'unique_int_1d', 'unique_masks_int_1d',
           'iscontiguous_int_1d', 'argwhere_int_1d', 'make_whole',
           'find_fragments']

cdef extern from "calc_distances.h":
    ctypedef float coordinate[3]
    void minimum_image(double* x, float* box, float* inverse_box)
    void minimum_image_triclinic(double* dx, float* box)

ctypedef cset[int] intset
ctypedef cmap[int, intset] intmap


def coords_add_vector(float[:, :] coordinates not None,
                      np.ndarray vector not None):
    """Add `vector` to each position in `coordinates`.

    Equivalent to ``coordinates += vector`` but faster for C-contiguous
    coordinate arrays. Coordinates are modified in place.

    Parameters
    ----------
    coordinates: numpy.ndarray
        Coordinate array of dtype ``numpy.float32`` and shape ``(n, 3)``.
    vector: numpy.ndarray
        Single coordinate vector of shape ``(3,)``.

    Raises
    ------
    ValueError
        If the shape of `coordinates` is not ``(n, 3)`` or if the shape of
        `vector` is not ``(3,)``.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t cshape1 = coordinates.shape[1]
    cdef np.intp_t vndim = vector.ndim
    cdef np.intp_t vshape0 = vector.shape[0]
    if cshape1 != 3:
        raise ValueError("Wrong shape: positions.shape != (n, 3)")
    if vndim != 1 or vshape0 != 3:
        raise ValueError("Wrong shape: vector.shape != (3,)")
    if vector.dtype == np.float32:
        _coords_add_vector32(coordinates, vector)
    else:
        _coords_add_vector64(coordinates, vector.astype(np.float64, copy=False))


cdef inline void _coords_add_vector32(float[:, :]& coordinates,
                                      float[:]& vector) nogil:
    """Low-level implementation of :func:`coords_add_vector` for float32
    vectors.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t i
    for i in range(coordinates.shape[0]):
        coordinates[i, 0] += vector[0]
        coordinates[i, 1] += vector[1]
        coordinates[i, 2] += vector[2]


cdef inline void _coords_add_vector64(float[:, :]& coordinates,
                                      double[:]& vector) nogil:
    """Low-level implementation of :func:`coords_add_vector` for non-float32
    vectors.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t i
    for i in range(coordinates.shape[0]):
        # Don't use += here!
        coordinates[i, 0] = coordinates[i, 0] + vector[0]
        coordinates[i, 1] = coordinates[i, 1] + vector[1]
        coordinates[i, 2] = coordinates[i, 2] + vector[2]


def unique_int_1d(np.intp_t[:] values not None, bint return_counts=False,
                  bint return_masks=False, bint assume_unsorted=False):
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
    return_masks: bool, optional
        If ``True``, an array of masks (with one mask per unique value) will be
        returned as well.
    assume_unsorted: bool, optional
        If `values` is known to be unsorted (i.e., its values are not
        monotonically increasing), setting `assume_unsorted` to ``True`` can
        speed up the computation.

    Returns
    -------
    unique: numpy.ndarray
        A deduplicated copy of `values`.
    counts: numpy.ndarray, optional
        An array of the same length as `unique` containing the number of
        occurrences of each unique value in the original `values` array. Only
        returned if `return_counts` is ``True``.
    masks: numpy.ndarray, optional
        An array of dtype ``object`` containing a mask for each value in
        `unique`. Each of the masks allows accessing all occurrences of its
        corresponding unique value in `values` such that
        ``numpy.all(values[masks[i]] == unique[i]) == True``. Thus, the masks
        array is roughly equivalent to
        ``[numpy.where(values == i) for i in numpy.unique(values)]``. Only
        returned if `return_masks` is ``True``.

    Notes
    -----
    The dtype ``numpy.intp`` is usually equivalent to ``numpy.int32`` on a 32
    bit operating system, and, likewise, equivalent to ``numpy.int64`` on a 64
    bit operating system. The exact behavior is compiler-dependent and can be
    checked with ``print(numpy.intp)``.


    See Also
    --------
    :func:`numpy.unique`
    :func:`unique_masks_int_1d`


    .. versionadded:: 0.19.0
    .. versionchanged:: 0.20.0
       Added optional  `return_counts`, `return_masks`, and `assume_unsorted`
       parameters and changed dtype from ``np.int64`` to ``np.intp``
       (corresponds to atom indices).
    """
    cdef np.intp_t n_values = values.shape[0]
    cdef np.intp_t n_unique
    cdef np.ndarray[np.intp_t, ndim=1] counts
    cdef np.ndarray[object, ndim=1] masks
    cdef np.ndarray[np.intp_t, ndim=1] unique = np.empty(n_values,
                                                         dtype=np.intp)
    if return_counts:
        counts = np.empty(n_values, dtype=np.intp)
        if return_masks:
            masks = np.empty(n_values, dtype=object)
            n_unique = _unique_int_1d_counts_masks(values, counts, masks,
                                                   assume_unsorted, unique)
            return unique[:n_unique], counts[:n_unique], masks[:n_unique]
        else:
            n_unique = _unique_int_1d_counts(values, counts, assume_unsorted,
                                             unique)
            return unique[:n_unique], counts[:n_unique]
    if return_masks:
        masks = np.empty(n_values, dtype=object)
        n_unique = _unique_int_1d_masks(values, masks, assume_unsorted, unique)
        return unique[:n_unique], masks[:n_unique]
    n_unique = _unique_int_1d(values, assume_unsorted, unique)
    return unique[:n_unique]


cdef inline np.intp_t _unique_int_1d(np.intp_t[:]& values, bint assume_unsorted,
                                     np.intp_t[::1] unique):
    """Low-level implementation of :func:`unique_int_1d` with
    ``return_counts=False``.


    .. versionadded:: 0.20.0
    """
    cdef bint monotonic = True
    cdef np.intp_t i = 0
    cdef np.intp_t n = 0
    cdef np.intp_t n_values = values.shape[0]
    cdef np.intp_t[::1] sorted_values
    if n_values > 0:
        if assume_unsorted and n_values > 1:
            sorted_values = np.sort(values)
            unique[0] = sorted_values[0]
            for i in range(1, n_values):
                if sorted_values[i] != unique[n]:
                    n += 1
                    unique[n] = sorted_values[i]
        else:
            unique[0] = values[0]
            if n_values > 1:
                for i in range(1, n_values):
                    if values[i] != unique[n]:
                        if monotonic and values[i] < unique[n]:
                            monotonic = False
                        n += 1
                        unique[n] = values[i]
                if not monotonic:
                    n_values = n + 1
                    n = 0
                    sorted_values = np.sort(unique[:n_values])
                    unique[0] = sorted_values[0]
                    for i in range(1, n_values):
                        if sorted_values[i] != unique[n]:
                            n += 1
                            unique[n] = sorted_values[i]
        n += 1
    return n


cdef inline np.intp_t _unique_int_1d_counts(np.intp_t[:]& values,
                                            np.intp_t[::1]& counts,
                                            bint assume_unsorted,
                                            np.intp_t[::1] unique):
    """Low-level implementation of :func:`unique_int_1d` with
    ``return_counts=True`` and ``return_masks=False``.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t i = 0
    cdef np.intp_t n = 0
    cdef np.intp_t n_values = values.shape[0]
    cdef np.intp_t[::1] sorted_values
    if n_values > 1:
        if not assume_unsorted:
            unique[0] = values[0]
            counts[0] = 1
            for i in range(1, n_values):
                if values[i] == unique[n]:
                    counts[n] += 1
                elif values[i] < unique[n]:
                    n = 0
                    break
                else:
                    n += 1
                    unique[n] = values[i]
                    counts[n] = 1
            else:
                return n + 1
        # values are unsorted or assume_unsorted == True:
        sorted_values = np.sort(values)
        unique[0] = sorted_values[0]
        counts[0] = 1
        for i in range(1, n_values):
            if sorted_values[i] == unique[n]:
                counts[n] += 1
            else:
                n += 1
                unique[n] = sorted_values[i]
                counts[n] = 1
        n += 1
    elif n_values == 1:
        n = 1
        unique[0] = values[0]
        counts[0] = 1
    return n


cdef inline np.intp_t _unique_int_1d_masks(np.intp_t[:]& values,
                                           object[::1]& masks,
                                           bint assume_unsorted,
                                           np.intp_t[::1] unique):
    """Low-level implementation of :func:`unique_int_1d` with
    ``return_counts=False`` and ``return_masks=True``.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t i = 0
    cdef np.intp_t n = 0
    cdef np.intp_t n_values = values.shape[0]
    cdef np.intp_t[::1] sort_ix
    cdef np.intp_t[::1] slice_ix
    if n_values > 1:
        slice_ix = np.empty(n_values + 1, dtype=np.intp)
        slice_ix[0] = 0
        if not assume_unsorted:
            unique[0] = values[0]
            for i in range(1, n_values):
                if values[i] != unique[n]:
                    if values[i] < unique[n]:
                        n = 0
                        break
                    n += 1
                    unique[n] = values[i]
                    slice_ix[n] = i
            else:
                n += 1
                slice_ix[n] = n_values
                for i in range(n):
                    masks[i] = slice(slice_ix[i], slice_ix[i + 1], 1)
                return n
        # values are unsorted or assume_unsorted == True:
        sort_ix = np.argsort(values)
        unique[0] = values[sort_ix[0]]
        for i in range(1, n_values):
            if values[sort_ix[i]] != unique[n]:
                n += 1
                unique[n] = values[sort_ix[i]]
                slice_ix[n] = i
        n += 1
        slice_ix[n] = n_values
        for i in range(n):
            masks[i] = sort_ix[slice_ix[i]:slice_ix[i + 1]]
    elif n_values == 1:
        n = 1
        unique[0] = values[0]
        masks[0] = slice(0, 1, 1)
    return n


cdef inline np.intp_t _unique_int_1d_counts_masks(np.intp_t[:]& values,
                                                  np.intp_t[::1]& counts,
                                                  object[::1]& masks,
                                                  bint assume_unsorted,
                                                  np.intp_t[::1] unique):
    """Low-level implementation of :func:`unique_int_1d` with
    ``return_counts=True`` and ``return_masks=True``.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t i = 0
    cdef np.intp_t n = 0
    cdef np.intp_t n_values = values.shape[0]
    cdef np.intp_t[::1] sort_ix
    cdef np.intp_t[::1] slice_ix
    if n_values > 1:
        slice_ix = np.empty(n_values + 1, dtype=np.intp)
        slice_ix[0] = 0
        counts[0] = 1
        if not assume_unsorted:
            unique[0] = values[0]
            for i in range(1, n_values):
                if values[i] == unique[n]:
                    counts[n] += 1
                elif values[i] < unique[n]:
                    n = 0
                    break
                else:
                    n += 1
                    unique[n] = values[i]
                    counts[n] = 1
                    slice_ix[n] = i
            else:
                n += 1
                slice_ix[n] = n_values
                for i in range(n):
                    masks[i] = slice(slice_ix[i], slice_ix[i + 1], 1)
                return n
        # values are unsorted or assume_unsorted == True:
        sort_ix = np.argsort(values)
        unique[0] = values[sort_ix[0]]
        counts[0] = 1
        for i in range(1, n_values):
            if values[sort_ix[i]] == unique[n]:
                counts[n] += 1
            else:
                n += 1
                unique[n] = values[sort_ix[i]]
                slice_ix[n] = i
                counts[n] = 1
        n += 1
        slice_ix[n] = n_values
        for i in range(n):
            masks[i] = sort_ix[slice_ix[i]:slice_ix[i + 1]]
    elif n_values == 1:
        n = 1
        unique[0] = values[0]
        counts[0] = 1
        masks[0] = slice(0, 1, 1)
    return n


def unique_masks_int_1d(np.intp_t[:] values not None,
                        bint assume_unsorted=False):
    """Find the indices of each unique element in a 1D array of integers and
    return them as an array of index masks or equivalent slices, similar to
    ``[numpy.where(values == i) for i in numpy.unique(values)]``.

    This function is optimal on sorted arrays.

    Parameters
    ----------
    values: numpy.ndarray
        1D array of dtype ``numpy.intp`` (or equivalent) in which to find the
        unique values.
    assume_unsorted: bool, optional
        If `values` is known to be unsorted (i.e., its values are not
        monotonically increasing), setting `assume_unsorted` to ``True`` can
        speed up the computation.

    Returns
    -------
    numpy.ndarray
        An array of dtype ``object`` containing index masks for each unique
        element in `values`.

    See Also
    --------
    :func:`unique_int_1d`


    .. versionadded:: 0.20.0
    """
    cdef object[::1] masks = _unique_masks_int_1d(values, assume_unsorted)
    return np.asarray(masks)


cdef inline object[::1] _unique_masks_int_1d(np.intp_t[:]& values,
                                             bint assume_unsorted):
    """Low-level implementation of :func:`unique_masks_int_1d`.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t i = 0
    cdef np.intp_t n = 1
    cdef np.intp_t n_values = values.shape[0]
    cdef np.intp_t[::1] sort_ix
    cdef np.intp_t[::1] slice_ix
    cdef object[::1] masks = np.empty(n_values, dtype=object)
    if n_values > 1:
        slice_ix = np.empty(n_values + 1, dtype=np.intp)
        slice_ix[0] = 0
        if not assume_unsorted:
            for i in range(1, n_values):
                if values[i] != values[i - 1]:
                    if values[i] < values[i - 1]:
                        n = 1
                        break
                    slice_ix[n] = i
                    n += 1
            else:
                slice_ix[n] = n_values
                for i in range(n):
                    masks[i] = slice(slice_ix[i], slice_ix[i + 1], 1)
                return masks[:n]
        # values are unsorted or assume_unsorted == True:
        sort_ix = np.argsort(values)
        for i in range(1, n_values):
            if values[sort_ix[i]] != values[sort_ix[i - 1]]:
                slice_ix[n] = i
                n += 1
        slice_ix[n] = n_values
        for i in range(n):
            masks[i] = sort_ix[slice_ix[i]:slice_ix[i + 1]]
        masks = masks[:n]
    elif n_values == 1:
        masks[0] = slice(0, 1, 1)
    return masks


def iscontiguous_int_1d(np.intp_t[:] values):
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
    if n_values > 1:
        for i in range(n_values - 1):
            if values[i] + 1 != values[i + 1]:
                return False
    elif n_values == 0:
        return False
    return True


def indices_to_slice_1d(np.ndarray[np.intp_t, ndim=1] indices):
    """Converts an index array to a slice if possible.

    Slice conversion is only performed if all indices are non-negative.

    Parameters
    ----------
    indices: numpy.ndarray
        1D array of dtype ``numpy.intp``.

    Returns
    -------
    slice or np.ndarray
        If `indices` can be represented by a slice and all its elements are
        non-negative, a slice object is returned. Otherwise, the `indices` array
        is returned.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t[:] ix = indices
    cdef np.intp_t n = ix.shape[0]
    cdef np.intp_t i
    cdef np.intp_t stride
    if n > 0 and ix[0] >= 0:
        if n > 1:
            stride = ix[1] - ix[0]
            if stride > 0:
                for i in range(1, n - 1):
                    if ix[i] + stride != ix[i + 1]:
                        return indices  # multiple strides
                return slice(ix[0], ix[n - 1] + 1, stride)
            elif stride < 0 and ix[n - 1] >= 0:
                for i in range(1, n - 1):
                    if ix[i] + stride != ix[i + 1]:
                        return indices  # multiple strides
                if ix[n - 1] == 0:
                    return slice(ix[0], None, stride)  # last index is 0
                return slice(ix[0], ix[n - 1] - 1, stride)
            return indices  # stride == 0
        return slice(ix[0], ix[0] + 1, 1)  # single index
    return indices  # empty array


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
    cdef np.ndarray[np.intp_t, ndim=1] result = np.empty(arr.shape[0], dtype=np.intp)
    cdef np.intp_t nargs = _argwhere_int_1d(arr, value, result)
    return result[:nargs]


cdef inline np.intp_t _argwhere_int_1d(np.intp_t[:] arr, np.intp_t value,
                                       np.intp_t[::1] result) nogil:
    """Low-level implementation of :func:`argwhere_int_1d`.


    .. versionadded:: 0.20.0
    """
    cdef np.intp_t i
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
