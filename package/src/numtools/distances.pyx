# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
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
Fast distance array computation --- :mod:`MDAnalysis.core.distances`
====================================================================

Fast C-routines to calculate distance arrays from coordinate arrays.

Overview
--------

.. function:: distance_array(ref,conf,[box,[,result]])

   Calculate all distances d_ij between the coordinates ref[i] and
   conf[j] in the numpy arrays *ref* and *conf*. If an orthorhombic
   *box* is supplied then a minimum image convention is used before
   calculating distances.

   If a 2D numpy array of dtype ``numpy.float64`` with the shape ``(len(ref),
   len(conf))`` is provided in *result* then this preallocated array is
   filled. This can speed up calculations.

.. function:: self_distance_array(ref,[box[,result]])

   Calculate all distances d_ij between atoms i and j in the reference
   coordinates *ref* for all N coordinates. Other options as in
   :func:`distance_array`.

   If a 1D numpy array of dtype ``numpy.float64`` with ``N*(N-1)/2`` elements is
   provided in *result* then this preallocated array is filled. This can speed
   up calculations.

Functions
---------
"""

cimport c_numpy
c_numpy.import_array()

cdef extern from "string.h":
    void* memcpy(void *dst, void *src, int len)

cdef extern from "calc_distances.h":
    ctypedef float coordinate[3]

    void calc_distance_array(coordinate* ref, int numref, coordinate* conf, int numconf, float* box, double* distances)
    void calc_distance_array_noPBC(coordinate* ref, int numref, coordinate* conf, int numconf, double* distances)
    void calc_self_distance_array(coordinate* ref, int numref, float* box, double* distances, int distnum)
    void calc_self_distance_array_noPBC(coordinate* ref, int numref, double* distances, int distnum)
    void coord_transform(coordinate* coords, int numCoords, coordinate* box)
    void calc_distance_array_triclinic(coordinate* ref, int numref, coordinate* conf, int numconf, coordinate* box, double* distances)
    void calc_self_distance_array_triclinic(coordinate* ref, int numref, coordinate* box, double* distances, int distnum)

import numpy
from MDAnalysis.coordinates.core import triclinic_vectors

def boxCheck(box):
    """Take a box input and deduce what type of system it represents based
    on the shape of the array and whether all angles are 90.
    
    :Returns:
      'ortho' orthogonal box
      'tri_vecs' triclinic box vectors
      'tri_box' triclinic box lengths and angles
      'unknown' boxCheck default, indicates no match found
    """
    boxtype = 'unknown'
    if box.shape == (3,):
        boxtype = 'ortho'
    elif box.shape == (3,3):
        if numpy.all([box[0][1] == 0.0, #Checks that tri box is properly formatted
                      box[0][2] == 0.0,
                      box[1][2] == 0.0]):
            boxtype = 'tri_vecs'
        else:
            boxtype = 'tri_vecs_bad'
    elif box.shape == (6,):
        if numpy.all(box[3:] == 90.):
            boxtype = 'ortho'
        else:
            boxtype = 'tri_box'

    return boxtype

def distance_array(c_numpy.ndarray reference, c_numpy.ndarray configuration, c_numpy.ndarray box=None, c_numpy.ndarray result=None):
    """Calculate all distances between a reference set and another configuration.

    d = distance_array(ref,conf[,box[,result=d]])

    :Arguments:
                *ref*
                        reference coordinate array
                *conf*
                        configuration coordinate array
                *box*
                        cell dimensions (minimum image convention is applied) or None [None]
                *result*
                        optional preallocated result array which must have the shape (len(ref),len(conf)) and dtype=numpy.float64. Avoids creating the              
                        array which saves time when the function is called repeatedly. [None]

    :Returns:
                *d*
                        (len(ref),len(conf)) numpy array with the distances d[i,j] between ref coordinates i and conf coordinates j

    .. Note:: This method is slower than it could be because internally we need to
          make copies of the ref and conf arrays.
    """
    cdef c_numpy.ndarray ref, conf
    cdef c_numpy.ndarray distances
    cdef int confnum, refnum

    # Work-around for a severe bug: function produces wrong numbers if
    # input arrays are views (eg slices from other arrays): copy to force a
    # new contiguous array in memory (and just make sure its in C order)
    ref = reference.copy('C')
    conf = configuration.copy('C')

    if (conf.nd != 2 or conf.dimensions[1] != 3):
        raise ValueError("conf must be a sequence of 3 dimensional coordinates")
    if (ref.nd != 2 or ref.dimensions[1] != 3):
        raise ValueError("ref must be a sequence of 3 dimensional coordinates")
    if (conf.dtype!=numpy.dtype(numpy.float32) or ref.dtype!=numpy.dtype(numpy.float32)):
        raise TypeError("coordinate data must be of type float32")

    with_PBC = (box is not None)
    if with_PBC:
        boxtype = boxCheck(box)
        if (boxtype == 'unknown'):
            raise ValueError("box input not recognised, must be an array of box dimensions")
        if (boxtype == 'tri_box'): # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box))
        if (box.dtype!=numpy.dtype(numpy.float32)):
            raise TypeError("periodic boundaries must be of type float32")

    confnum = conf.dimensions[0]
    refnum = ref.dimensions[0]

    if not result is None:
        if (result.nd != 2 or result.dimensions[0] != refnum or result.dimensions[1] != confnum):
            raise ValueError("result array has incorrect size - should be (%dx%d)"%(refnum,confnum))
        if (result.dtype != numpy.dtype(numpy.float64)):
            raise TypeError("result array must be of type numpy.float64")
        distances = numpy.asarray(result)
    else:
        distances = numpy.zeros((refnum, confnum), numpy.float64)

    if with_PBC:
        if boxtype == 'ortho':
            calc_distance_array(<coordinate*>ref.data, refnum, <coordinate*>conf.data, confnum, 
                                <float*>box.data, <double*>distances.data)
        else: 
            calc_distance_array_triclinic(<coordinate*> ref.data, refnum, <coordinate*> conf.data, confnum, 
                                          <coordinate*>box.data, <double*> distances.data)
    else:
        calc_distance_array_noPBC(<coordinate*>ref.data, refnum, <coordinate*>conf.data, confnum, <double*>distances.data)

    return distances

def self_distance_array(c_numpy.ndarray reference, c_numpy.ndarray box=None, c_numpy.ndarray result=None):
    """Calculate all distances d_ij between atoms i and j within a configuration *ref*.

    d = self_distance_array(ref[,box[,result=d]])

    :Arguments:
                *ref*
                        reference coordinate array with N=len(ref) coordinates
                *box*
                        cell dimensions (minimum image convention is applied) or None [None]
                *result*
                        optional preallocated result array which must have the shape
                           (N*(N-1)/2,) and dtype ``numpy.float64``. Avoids creating
                           the array which saves time when the function is called repeatedly. [None]

    :Returns:
                *d*
                        N*(N-1)/2 numpy 1D array with the distances dist[i,j] between ref
                           coordinates i and j at position d[k]. Loop through d::

                             for i in xrange(N):
                                for j in xrange(i+1, N):
                                    k += 1
                                    dist[i,j] = d[k]

    .. Note:: This method is slower than it could be because internally we need to
          make copies of the coordinate arrays.
    """

    cdef c_numpy.ndarray ref
    cdef c_numpy.ndarray distances
    cdef int refnum, distnum

    # Work-around for a severe bug: function produces wrong numbers if
    # input arrays are views (eg slices from other arrays): copy to force a
    # new contiguous array in memory (and just make sure its in C order)
    ref = reference.copy('C')

    if (ref.nd != 2 or ref.dimensions[1] != 3):
        raise ValueError("ref must be a sequence of 3 dimensional coordinates")
    if (ref.dtype!=numpy.dtype(numpy.float32)):
        raise TypeError("coordinate data must be of type float32")

    with_PBC = (box is not None)
    if with_PBC:
        boxtype = boxCheck(box)
        if (boxtype == 'unknown'):
            raise ValueError("box input not recognised, must be an array of box dimensions")
        if (boxtype == 'tri_box'): # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box))
        if (box.dtype!=numpy.dtype(numpy.float32)):
            raise TypeError("periodic boundaries must be of type float32")

    refnum = ref.dimensions[0]
    distnum = (refnum*(refnum-1))/2

    if not result is None:
        if (result.nd != 1 or result.dimensions[0] != distnum):
            raise ValueError("result array has incorrect size or datatype - should be (%d)"%(distnum))
        if (result.dtype != numpy.dtype(numpy.float64)):
            raise TypeError("result array must be of type numpy.float64")
        distances = numpy.asarray(result)
    else:
        distances = numpy.zeros((distnum,), numpy.float64)

    if with_PBC:
        if boxtype == 'ortho':
            calc_self_distance_array(<coordinate*>ref.data,refnum,<float*>box.data,
                                     <double*>distances.data, distnum)
        else:
            calc_self_distance_array_triclinic(<coordinate*>ref.data,refnum,<coordinate*>box.data,
                                               <double*>distances.data, distnum)
    else:
        calc_self_distance_array_noPBC(<coordinate*>ref.data,refnum,<double*>distances.data,distnum)

    return distances

def transform_RtoS(c_numpy.ndarray inputcoords, c_numpy.ndarray box):
    """Transform an array of coordinates from real space to S space (aka lambda space)

    S space represents fractional space within the unit cell for this system

    Reciprocal operation to :meth:`transform_StoR

    :Arguments:
      *inputcoords*
                      An n x 3 array of coordinate data, of type np.float32
      *box*
                      The unitcell dimesions for this system

    :Returns:
       *outcoords*
                      An n x 3 array of fracional coordiantes

    """

    cdef c_numpy.ndarray coords, inv
    cdef int numcoords

    #Create contiguous array
    coords = inputcoords.copy('C')
    numcoords = coords.dimensions[0]

    boxtype = boxCheck(box)
    if (boxtype == 'unknown'):
        raise ValueError("box input not recognised, must be an array of box dimensions")

    if (boxtype == 'tri_box'): # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
        box = triclinic_vectors(box)
    if (boxtype == 'tri_vecs_bad'):
        box = triclinic_vectors(triclinic_box(box))
    elif (boxtype == 'ortho'):
        box = numpy.array([[box[0], 0.0, 0.0],
                           [0.0, box[1], 0.0],
                           [0.0, 0.0, box[2]]], dtype = numpy.float32)

    #Create inverse matrix of box
    inv = numpy.array(numpy.matrix(box).I, dtype = numpy.float32, order='C') # need order C here

    # Checks on input arrays
    if (box.dtype!=numpy.dtype(numpy.float32)):
        raise TypeError("periodic boundaries must be of type float32")
    if (box.dimensions[0] != 3 or box.dimensions[1] != 3):
        raise ValueError("Box format not recognised, please use system dimensions")
    if (coords.nd != 2 or coords.dimensions[1] != 3):
        raise ValueError("Coordinates must be a sequence of 3 dimensional coordinates")
    if (coords.dtype != numpy.dtype(numpy.float32)):
        raise TypeError("Coordinate data must be of type numpy float32")
    
    coord_transform(<coordinate*> coords.data, numcoords, <coordinate*> inv.data)

    return coords

def transform_StoR(c_numpy.ndarray inputcoords, c_numpy.ndarray box):
    """Transform an array of coordinates from S space into real space.

    S space represents fractional space within the unit cell for this system    

    Reciprocal operation to :meth:`transform_RtoS`

    :Arguments:
      *inputcoords*
                      An n x 3 array of coordinate data, of type np.float32
      *box*
                      The unitcell dimesions for this system

    :Returns:
       *outcoords*
                      An n x 3 array of fracional coordiantes
    """

    cdef c_numpy.ndarray coords
    cdef int numcoords

    #Create contiguous array
    coords = inputcoords.copy('C')
    numcoords = coords.dimensions[0]

    boxtype = boxCheck(box)
    if (boxtype == 'unknown'):
        raise ValueError("box input not recognised, must be an array of box dimensions")

    if (boxtype == 'tri_box'): # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
        box = triclinic_vectors(box)
    elif (boxtype == 'ortho'):
        box = numpy.array([[box[0], 0.0, 0.0],
                           [0.0, box[1], 0.0],
                           [0.0, 0.0, box[2]]], dtype = numpy.float32)

    # Checks on input arrays
    if (box.dtype!=numpy.dtype(numpy.float32)):
        raise TypeError("periodic boundaries must be of type float32")
    if (coords.nd != 2 or coords.dimensions[1] != 3):
        raise ValueError("Coordinates must be a sequence of 3 dimensional coordinates")
    if (coords.dtype != numpy.dtype(numpy.float32)):
        raise TypeError("Coordinate data must be of type numpy float32")
    
    coord_transform(<coordinate*> coords.data, numcoords, <coordinate*> box.data)

    return coords
