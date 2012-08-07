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

# Register a np.float64 as data type 'DTYPE_t'
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

# Register a C math sqrt function
cdef extern from "math.h":
    float sqrt(double x) nogil
    float fabs(double x) nogil

def distance_array_serial(np.ndarray[DTYPE_t, ndim=2] coordA, \
                          np.ndarray[DTYPE_t, ndim=2] coordB, \
                          np.ndarray[DTYPE_t, ndim=1] box = None, \
                          np.ndarray[DTYPE_t, ndim=2] result = None):
    """distance_array_serial(ref,conf)

    Calculate all distances d_ij between the coordinates ref[i] and
    conf[j] in the numpy arrays *ref* and *conf*.

    Serial version (to check the parallel version). This version is
    slightly slower than the regular serial (pure C)
    :func:`MDAnalysis.core.distances.distance_array` function.
    """
    cdef DTYPE_t x, y, z, dist
    cdef Py_ssize_t i, j
    
    cdef np.ndarray[DTYPE_t, ndim=1] box_half

    cdef char has_box = 0
    cdef DTYPE_t box_x, box_y, box_z, box_half_x, box_half_y, box_half_z

    if box != None:
      has_box = 1
      box_half = box / 2
      
      box_x = box[0]
      box_y = box[1]
      box_z = box[2]

      box_half_x = box_half[0]
      box_half_y = box_half[1]
      box_half_z = box_half[2]

    
    rows = coordA.shape[0];
    cols = coordB.shape[0];
    assert rows == cols, """Coordinate arrays of the same length must be used. 
    Distance matrix must be square: number of rows (%d) must be the same as the number of columns (%d)""" % (rows, cols)

    if result == None:
      result = np.empty((rows, cols), dtype=DTYPE)
    else: 
      assert result.shape[0] == rows, "Results array should have %d rows, has %d" % (rows, result.shape[0])
      assert result.shape[1] == cols, "Results array should have %d columns, has %d" % (cols, result.shape[1])

    # The two loops are independent, let's use
    for i in range(rows) :
        for j in range(cols) :

            x = coordA[i,0] - coordB[j,0];
            y = coordA[i,1] - coordB[j,1];
            z = coordA[i,2] - coordB[j,2];

            # Python objects, including np.ndarrays (even when defined)
            # cannot be indexed. This has been changed in Cython 0.17-beta1
            #
            if has_box == 1:
                if fabs(x) > box_half_x:
                    if x < 0.0: x = x + box_x
                    else: x = x - box_x
                if fabs(y) > box_half_y:
                    if y < 0.0: y = y + box_y
                    else: y = y - box_y
                if fabs(z) > box_half_z:
                    if z < 0.0: z = z + box_z
                    else: z = z - box_z  

            dist = sqrt((x*x)+(y*y)+(z*z));

            result[i,j] = dist;

    return result

# Jan: minimum_image has been dopted from calc_distances.h - there it's using 
# doubles for positions and floats for box/box_half, while I've used float32 for
# both; I'll be happy to change, if there is a convention.

    

@cython.boundscheck(False)
def distance_array(np.ndarray[DTYPE_t, ndim=2] coordA, \
                   np.ndarray[DTYPE_t, ndim=2] coordB, \
                   np.ndarray[DTYPE_t, ndim=1] box = None, \
                   np.ndarray[DTYPE_t, ndim=2] result = None):
                   
    """distance_array(ref,conf,box=None,result=None)

    Calculate all distances d_ij between the coordinates ref[i] and
    conf[j] in the numpy arrays *ref* and *conf*.

    Parallel version that will automatically decide on how many threads
    to run.
    """
    cdef DTYPE_t x, y, z, dist
    cdef Py_ssize_t i, j
    
    cdef np.ndarray[DTYPE_t, ndim=1] box_half

    cdef char has_box = 0
    cdef DTYPE_t box_x, box_y, box_z, box_half_x, box_half_y, box_half_z

    if box != None:
      has_box = 1
      box_half = box / 2
      
      box_x = box[0]
      box_y = box[1]
      box_z = box[2]

      box_half_x = box_half[0]
      box_half_y = box_half[1]
      box_half_z = box_half[2]

    # FIXME assume that coorA and coorB are of tha same length
    rows = coordA.shape[0];
    cols = coordB.shape[0];

    if result == None:
      result = np.empty((rows, cols), dtype=DTYPE)
    else: 
      assert result.shape[0] == rows, "Results array should have %d rows, has %d" % (rows, result.shape[0])
      assert result.shape[1] == cols, "Results array should have %d columns, has %d" % (cols, result.shape[1])
    
    with nogil, parallel():
        # The two loops are independent, let's use
        for i in prange(rows, schedule="dynamic", chunksize=50) :
            for j in range(cols) :
                
                x = coordA[i,0] - coordB[j,0];
                y = coordA[i,1] - coordB[j,1];
                z = coordA[i,2] - coordB[j,2];
                
                # Python objects, including np.ndarrays (even when defined)
                # cannot be indexed. This has been changed in Cython 0.17-beta1
                #
                if has_box == 1:
                    if fabs(x) > box_half_x:
                        if x < 0.0: x = x + box_x
                        else: x = x - box_x
                    if fabs(y) > box_half_y:
                        if y < 0.0: y = y + box_y
                        else: y = y - box_y
                    if fabs(z) > box_half_z:
                        if z < 0.0: z = z + box_z
                        else: z = z - box_z                        
                
                # FIXME this might not be the optimal thing to do
                #
                dist = sqrt((x*x)+(y*y)+(z*z));

                result[i,j] = dist;
    return result

