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

from libcpp.set cimport set as cset
from libcpp.map cimport map as cmap

__all__ = ['unique_int_1d', '_is_contiguous']


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def unique_int_1d(long[:] values):
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
    cdef long[:] result = np.empty(n_values, dtype=np.int64)

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
        result.sort()
        result = unique_int_1d(result)
    return result


ctypedef cset[int] intset
ctypedef cmap[int, intset] intmap

@cython.boundscheck(False)
@cython.wraparound(False)
cdef intset difference(intset a, intset b):
    """a.difference(b)

    Returns set of values in a which are not in b
    """
    cdef intset output

    output = intset()

    for val in a:
        if b.count(val) != 1:
            output.insert(val)

    return output


@cython.boundscheck(False)
@cython.wraparound(False)
def _is_contiguous(int[:] atoms, int[:, :] bonds, int start):
    """

    Parameters
    ----------
    atoms : np.ndarray
        array of atom indices to consider
    bonds : np.ndarray
        array of bonds
    start : int
        where to start walking

    Returns
    -------
    result : bool
        If entire molecule can be traversed by bonds
    """
    cdef bint result
    cdef intset seen, done, todo, total
    cdef intmap bonding
    cdef int i, N, nloops
    cdef int x, y

    total = intset()
    bonding = intmap()

    # make set of which atoms exist
    N = atoms.shape[0]
    for i in range(N):
        total.insert(atoms[i])

    if not total.count(start):
        raise ValueError

    # build C++ dict of bonds
    N = bonds.shape[0]
    for i in range(N):
        x = bonds[i, 0]
        y = bonds[i, 1]
        # only add bonds if both atoms are in atoms set
        if total.count(x):
            if total.count(y):
                bonding[x].insert(y)
                bonding[y].insert(x)

    seen = intset()
    seen.insert(start)
    done = intset()

    N = total.size()

    nloops = 0
    while seen.size() < N:
        nloops += 1
        if nloops >= N:
            break
        # todo is set of start points
        # can start on anyone that has been seen, but not done yet
        todo = difference(seen, done)

        for x in todo:  # for each start point
            for y in bonding[x]:  # add all bonded atoms
                seen.insert(y)
            # mark as done
            done.insert(x)

    # if we saw all Atoms when walking, is_contiguous
    if seen.size() == N:
        result = True
    else:
        result = False

    return result
