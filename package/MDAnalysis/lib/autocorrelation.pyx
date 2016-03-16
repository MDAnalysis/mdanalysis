# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2014 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
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

"""Vector autocorrelation functions --- :mod:`MDAnalysis.lib.autocorrelation
============================================================================

.. autofunction:: vector_autocorrelation(vectors)
.. autofunction:: windowed_vector_autocorrelation(vectors, window)

"""
import numpy as np
cimport cython
cimport numpy as np


@cython.boundscheck(False)
@cython.wraparound(False)
def vector_autocorrelation(np.ndarray[np.float32_t, ndim=2] vectors):
    """Perform a 3D vector autocorrelation

    .. math:: C_{VV}(t) = \langle V(t_0 + t) \cdot V(t_0) \rangle

    Parameters
    ----------
    vectors : ndarray
      Shape (n, 3) array of vectors.  Dtype must be float32

    Returns
    -------
    results : ndarray
      Shape (n) array of average autocorrelation

    Note
    ----
    This algorithm scales n^2, for a faster approach for larger
    data, see windowed_vector_autocorrelation.

    .. versionadded:: 0.15.0
    """
    cdef Py_ssize_t i, j, k, pos
    cdef float res
    cdef int n

    n = vectors.shape[0]

    cdef np.ndarray[np.float32_t] results = np.zeros(n, dtype=np.float32)

    for i in range(n):
        pos = 0
        for j in range(i, n):
            res = 0.0
            for k in range(3):
                res += vectors[i, k] * vectors[j, k]
            results[pos] += res
            pos += 1

    for i in range(n):
        results[i] /= n - i

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
def windowed_vector_autocorrelation(np.ndarray[np.float32_t, ndim=2] vectors,
                                    int window):
    """Perform a 3D vector autocorrelation over a rolling window

    Parameters
    ----------
    vectors : ndarray
      Shape (n, 3) array of vectors
    window : int
      The size of the window to perform the autocorrelation over

    Returns
    -------
    Results : ndarray
      Shape(window) of the average vector autocorrelation

    .. versionadded:: 0.15.0
    """
    cdef Py_ssize_t i, j, k
    cdef float res
    cdef int n

    n = vectors.shape[0] - window + 1

    cdef np.ndarray[np.float32_t] results = np.zeros(window, dtype=np.float32)

    for i in range(n):
        for j in range(window):
            res = 0.0
            for k in range(3):
                res += vectors[i, k] * vectors[i + j, k]
            results[j] += res

    for i in range(window):
        results[i] /= n

    return results
