# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; encoding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
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
Fast parallel distance array computation --- :mod:`MDAnalysis.core.parallel.distances`
======================================================================================

:Author:  Jan Domański
:Year:    2012
:Licence: GPL

A fast, parallel :func:`distance_array` function as a substitute for
:func:`MDAnalysis.core.distances.distance_array`; implemented with
`Cython Parallelism`_. For development notes see the comments for
`Issue 80`_.

.. _Cython Parallelism: http://docs.cython.org/src/userguide/parallelism.html
.. _Issue 80: https://code.google.com/p/mdanalysis/issues/detail?id=80

Load the module with ::

  import MDAnalysis.core.parallel.distances


.. function:: distance_array(ref,conf)

   Calculate all distances d_ij between the coordinates ref[i] and
   conf[j] in the numpy arrays *ref* and *conf*.

   Parallel version that will automatically decide on how many threads
   to run.

   .. versionadded:: 0.8

.. function:: distance_array_serial(ref,conf)

   Calculate all distances d_ij between the coordinates ref[i] and
   conf[j] in the numpy arrays *ref* and *conf*.

   Serial version (to check the parallel version). This version is
   slightly slower than the regular serial (pure C)
   :func:`MDAnalysis.core.distances.distance_array` function.

   .. versionadded:: 0.8
"""

cimport numpy as np
cimport cython
from cython.parallel import parallel, prange

import numpy as np

# Register a np.float32 as data type 'DTYPE_t'
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

# Register a C math sqrt function
cdef extern from "math.h":
    float sqrt(double x) nogil


def distance_array_serial(np.ndarray[DTYPE_t, ndim=2] coordA, np.ndarray[DTYPE_t, ndim=2]  coordB):
    """distance_array_serial(ref,conf)

    Calculate all distances d_ij between the coordinates ref[i] and
    conf[j] in the numpy arrays *ref* and *conf*.

    Serial version (to check the parallel version). This version is
    slightly slower than the regular serial (pure C)
    :func:`MDAnalysis.core.distances.distance_array` function.
    """
    cdef DTYPE_t x, y, z, dist
    cdef Py_ssize_t i, j
    cdef np.ndarray[DTYPE_t, ndim=2] distance

    # FIXME assume that coorA and coorB are of tha same length
    rows = coordA.shape[0];
    cols = coordB.shape[0];
    distance = np.empty((rows, cols), dtype=DTYPE)

    # The two loops are independent, let's use
    for i in range(rows) :
        for j in range(cols) :

            #if i == j:
            #  distance[i,j] = 0.0;
            #  continue


            x = coordA[i,0] - coordB[j,0];
            y = coordA[i,1] - coordB[j,1];
            z = coordA[i,2] - coordB[j,2];

            dist = sqrt((x*x)+(y*y)+(z*z));

            distance[i,j] = dist;

    return distance


@cython.boundscheck(False)
def distance_array(np.ndarray[DTYPE_t, ndim=2] coordA, np.ndarray[DTYPE_t, ndim=2]  coordB):
    """distance_array(ref,conf)

    Calculate all distances d_ij between the coordinates ref[i] and
    conf[j] in the numpy arrays *ref* and *conf*.

    Parallel version that will automatically decide on how many threads
    to run.
    """
    cdef DTYPE_t x, y, z, dist
    cdef Py_ssize_t i, j
    cdef np.ndarray[DTYPE_t, ndim=2] distance

    # FIXME assume that coorA and coorB are of tha same length
    rows = coordA.shape[0];
    cols = coordB.shape[0];
    distance = np.empty((rows, cols), dtype=DTYPE)
    with nogil, parallel():
        # The two loops are independent, let's use
        for i in prange(rows, schedule="dynamic", chunksize=50) :
            for j in range(cols) :

                #if i == j:
                #  distance[i,j] = 0.0;
                #  continue

                x = coordA[i,0] - coordB[j,0];
                y = coordA[i,1] - coordB[j,1];
                z = coordA[i,2] - coordB[j,2];

                # FIXME this might not be the optimal thing to do
                dist = sqrt((x*x)+(y*y)+(z*z));

                distance[i,j] = dist;
    return distance

