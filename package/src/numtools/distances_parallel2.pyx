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

.. autofunction:: boxCheck(box)
.. autofunction:: calc_bonds(atom1, atom2 [, box, [,result]])
.. autofunction:: calc_angles(atom1, atom2, atom3 [,box [, result]])
.. autofunction:: calc_torsions(atom1, atom2, atom3, atom4 [,box [, result]])
.. autofunction:: applyPBC(coordinates, box)
.. autofunction:: transform_RtoS(coordinates, box)
.. autofunction:: transform_StoR(coordinates, box)
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
    void calc_bond_distance(coordinate* atom1, coordinate* atom2, int numatom, float*box, double* distances)
    void calc_bond_distance_triclinic(coordinate* atom1, coordinate* atom2, int numatom, coordinate* box, double* distances)
    void calc_bond_distance_noPBC(coordinate* atom1, coordinate* atom2, int numatom, double* distances)
    void calc_angle(coordinate* atom1, coordinate* atom2, coordinate* atom3, int numatom, double* angles)
    void calc_angle_ortho(coordinate* atom1, coordinate* atom2, coordinate* atom3, int numatom, float* box, double* angles)
    void calc_angle_triclinic(coordinate* atom1, coordinate* atom2, coordinate* atom3, int numatom, coordinate* box, double* angles)
    void calc_torsion(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, int numatom, double* angles)
    void calc_torsion_ortho(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, int numatom, float* box, double* angles)
    void calc_torsion_triclinic(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, int numatom, coordinate* box, double* angles)
    void ortho_pbc(coordinate* coords, int numcoords, float* box, float* box_inverse)
    void triclinic_pbc(coordinate* coords, int numcoords, coordinate* box, float* box_inverse)

import numpy
from MDAnalysis.coordinates.core import triclinic_vectors, triclinic_box

def boxCheck(box):
    """Take a box input and deduce what type of system it represents based
    on the shape of the array and whether all angles are 90.
    
    :Arguments:
      *box*
          box information of unknown format

    :Returns:
      * ``ortho`` orthogonal box
      * ``tri_vecs`` triclinic box vectors
      * ``tri_box`` triclinic box lengths and angles
      * ``unknown`` boxCheck default, indicates no match found
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
                        optional preallocated result array which must have the shape (len(ref),
                        len(conf)) and dtype=numpy.float64. Avoids creating the              
                        array which saves time when the function is called repeatedly. [None]

    :Returns:
                *d*
                        (len(ref),len(conf)) numpy array with the distances d[i,j] 
                        between ref coordinates i and conf coordinates j

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
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))
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
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))
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

    Reciprocal operation to :meth:`transform_StoR`

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
        box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))
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

def calc_bonds(c_numpy.ndarray list1, c_numpy.ndarray list2, c_numpy.ndarray box=None, c_numpy.ndarray result=None):
    """
    Calculate all distances between a pair of atoms.  *atom1* and *atom2* are both
    arrays of coordinates, where atom1[i] and atom2[i] represent a bond.  

    In comparison to distance_array and self_distance_array which calculate distances
    between all combinations of coordinates, calc_bonds can be used to calculate distance
    between pairs of objects, similar to::
    
       numpy.linalg.norm(a - b) for a, b in zip(coords1, coords2)

    The optional argument *box* applies minimum image convention if supplied.
    *box* can be either orthogonal or triclinic
    
    If a 1D numpy array of dtype ``numpy.float64`` with ``len(atom1)`` elements is
    provided in *result* then this preallocated array is filled. This can speed
    up calculations.

    bondlengths = calc_bonds(coords1, coords2 [, box [,result=bondlengths]])

    :Arguments:
       *coords1*
          An array of coordinates for one half of the bond
       *coords2*
          An array of coordinates for the other half of bond
       *box*
          Unit cell information if periodic boundary conditions are required [None]
       *result*
          optional preallocated result array which must be same length as coord 
          arrays and dtype=numpy.float64. Avoids creating the              
          array which saves time when the function is called repeatedly. [None]

    :Returns:
       *bondlengths*
          Numpy array with the length between each pair in coords1 and coords2

    .. versionadded:: 0.8
    """
    cdef c_numpy.ndarray atom1, atom2
    cdef c_numpy.ndarray distances
    cdef int numatom

    atom1 = list1.copy('C')
    atom2 = list2.copy('C')

    if (atom1.nd != 2 or atom1.dimensions[1] != 3):
        raise ValueError("list1 must be a sequence of 3 dimensional coordinates")
    if (atom1.dtype!=numpy.dtype(numpy.float32)):
        raise TypeError("coordinate data must be of type float32")
    if (atom2.nd != 2 or atom2.dimensions[1] != 3):
        raise ValueError("list2 must be a sequence of 3 dimensional coordinates")
    if (atom2.dtype!=numpy.dtype(numpy.float32)):
        raise TypeError("coordinate data must be of type float32")

    if (atom1.dimensions[0] != atom2.dimensions[0]):
        raise ValueError("list1 and list2 of different size")

    with_PBC = (box is not None)
    if with_PBC:
        boxtype = boxCheck(box)
        if (boxtype == 'unknown'):
            raise ValueError("box input not recognised, must be an array of box dimensions")
        if (boxtype == 'tri_box'): # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))
        if (box.dtype!=numpy.dtype(numpy.float32)):
            raise TypeError("periodic boundaries must be of type float32")

    numatom = atom1.dimensions[0]

    if not result is None:
        if (result.nd != 1 or result.dimensions[0] != numatom):
            raise ValueError("result array has incorrect size - should be (%d)"%(numatom))
        if (result.dtype != numpy.dtype(numpy.float64)):
            raise TypeError("result array must be of type numpy.float64")
        distances = numpy.asarray(result)
    else:
        distances = numpy.zeros((numatom,), numpy.float64)

    if with_PBC:
        if boxtype == 'ortho':
            calc_bond_distance(<coordinate*>atom1.data,<coordinate*>atom2.data,numatom,<float*>box.data,<double*>distances.data)
        else:
            calc_bond_distance_triclinic(<coordinate*>atom1.data, <coordinate*>atom2.data, numatom, <coordinate*>box.data, <double*> distances.data)
    else:
        calc_bond_distance_noPBC(<coordinate*>atom1.data,<coordinate*>atom2.data,numatom,<double*>distances.data)

    return distances


def calc_angles(c_numpy.ndarray list1, c_numpy.ndarray list2, c_numpy.ndarray list3, 
                c_numpy.ndarray box=None, c_numpy.ndarray result=None):
    """
    Calculates the angle formed between three atoms, over a list of coordinates.
    All *atom* inputs are lists of coordinates of equal length, with *atom2* 
    representing the apex of the angle.

    If a 1D numpy array of dtype ``numpy.float64`` with ``len(atom1)`` elements is
    provided in *result* then this preallocated array is filled. This can speed
    up calculations.

    The optional argument ``box`` ensures that periodic boundaries are taken into account when
    constructing the connecting vectors between atoms, ie that the vector between atoms 1 & 2
    goes between coordinates in the same image.

    angles = calc_angles(coords1, coords2, coords3, [[box=None],result=angles])

    :Arguments:
        *coords1*
            coordinate array of one side of angles
        *coords2*
            coordinate array of apex of angles
        *coords3*
            coordinate array of other side of angles
        *box*
            optional unit cell information.  This ensures that the connecting vectors between
            atoms respect minimum image convention.  This is import when the angle might
            be between atoms in different images.  
        *result*
            optional preallocated results array which must have same length as coordinate 
            array and dtype=numpy.float64. 

    :Returns:
        *angles*
            A numpy.array of angles in radians
    
    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Added optional box argument to account for periodic boundaries in calculation
    """
    cdef c_numpy.ndarray atom1, atom2, atom3
    cdef c_numpy.ndarray angles
    cdef int numatom

    atom1 = list1.copy('C')
    atom2 = list2.copy('C')
    atom3 = list3.copy('C')
    numatom = atom1.dimensions[0]

    #checks on input arrays
    if (atom1.nd != 2 or atom1.dimensions[1] != 3):
        raise ValueError("list1 must be an array of 3 dimensional coordinates")
    if (atom2.nd != 2 or atom2.dimensions[1] != 3):
        raise ValueError("list2 must be an array of 3 dimensional coordinates")
    if (atom3.nd != 2 or atom3.dimensions[1] != 3):
        raise ValueError("list3 must be an array of 3 dimensional coordinates")
    if (atom2.dimensions[0] != numatom or atom3.dimensions[0] != numatom):
        raise ValueError("all lists must be the same length")
    #type check
    if (atom1.dtype != numpy.dtype(numpy.float32) or atom2.dtype != numpy.dtype(numpy.float32) 
        or atom3.dtype != numpy.dtype(numpy.float32) ):
        raise TypeError("all coordinates must be of type numpy.float32")

    with_PBC = (box is not None)
    if with_PBC:
        boxtype = boxCheck(box)
        if (boxtype == 'unknown'):
            raise ValueError("box input not recognised, must be an array of box dimensions")
        if (boxtype == 'tri_box'): # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))
        if (box.dtype!=numpy.dtype(numpy.float32)):
            raise TypeError("periodic boundaries must be of type float32")

    if not result is None:
        if (result.nd != 1 or result.dimensions[0] != numatom):
            raise ValueError("result array has incorrect size - should be (%d)"%(numatom))
        if (result.dtype != numpy.dtype(numpy.float64)):
            raise TypeError("result array must be of type numpy.float64")
        angles = numpy.asarray(result)
    else:
        angles = numpy.zeros((numatom,), numpy.float64)

    if with_PBC:
        if boxtype == 'ortho':
            calc_angle_ortho(<coordinate*>atom1.data,<coordinate*>atom2.data,<coordinate*>atom3.data,numatom,
                              <float*>box.data, <double*>angles.data)
        else:
            calc_angle_triclinic(<coordinate*>atom1.data,<coordinate*>atom2.data,<coordinate*>atom3.data,numatom,
                                  <coordinate*>box.data, <double*>angles.data)
    else:
        calc_angle(<coordinate*>atom1.data,<coordinate*>atom2.data,<coordinate*>atom3.data,numatom,<double*>angles.data)

    return angles

def calc_torsions(c_numpy.ndarray list1, c_numpy.ndarray list2, c_numpy.ndarray list3, c_numpy.ndarray list4, 
                  c_numpy.ndarray box=None, c_numpy.ndarray result=None):
    """
    Calculate the torsional angle formed by four atoms, over a list of coordinates.
    
    Torsional angle around axis connecting atoms 1 and 2 (i.e. the angle
    between the planes spanned by atoms (0,1,2) and (1,2,3))::

                  3
                  |
            1-----2
           /
          0

    If a 1D numpy array of dtype ``numpy.float64`` with ``len(atom1)`` elements is
    provided in *result* then this preallocated array is filled. This can speed
    up calculations.

    The optional argument ``box`` ensures that periodic boundaries are taken into account when
    constructing the connecting vectors between atoms, ie that the vector between atoms 1 & 2
    goes between coordinates in the same image.

    angles = calc_torsions(coords1, coords2, coords3, coords4 [,box=box, result=angles])

    :Arguments:
        *coords1*
            coordinate array of 1st atom in torsions
        *coords2*
            coordinate array of 2nd atom in torsions
        *coords3*
            coordinate array of 3rd atom in torsions
        *coords4*
            coordinate array of 4th atom in torsions
        *box*
            optional unit cell information.  This ensures that the connecting vectors between
            atoms respect minimum image convention.  This is import when the angle might
            be between atoms in different images.  
        *result*
            optional preallocated results array which must have same length as coordinate 
            array and dtype=numpy.float64. 

    :Returns:
        *angles*
            A numpy.array of angles in radians
    
    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Added optional box argument to account for periodic boundaries in calculation
    """
    cdef c_numpy.ndarray atom1, atom2, atom3, atom4
    cdef c_numpy.ndarray angles
    cdef int numatom

    atom1 = list1.copy('C')
    atom2 = list2.copy('C')
    atom3 = list3.copy('C')
    atom4 = list4.copy('C')
    numatom = atom1.dimensions[0]

    #checks on input arrays
    if (atom1.nd != 2 or atom1.dimensions[1] != 3):
        raise ValueError("list1 must be an array of 3 dimensional coordinates")
    if (atom2.nd != 2 or atom2.dimensions[1] != 3):
        raise ValueError("list2 must be an array of 3 dimensional coordinates")
    if (atom3.nd != 2 or atom3.dimensions[1] != 3):
        raise ValueError("list3 must be an array of 3 dimensional coordinates")
    if (atom4.nd != 2 or atom4.dimensions[1] != 3):
        raise ValueError("list3 must be an array of 3 dimensional coordinates")
    if (atom2.dimensions[0] != numatom or atom3.dimensions[0] != numatom or atom4.dimensions[0] != numatom):
        raise ValueError("all lists must be the same length")
   #type check
    if (atom1.dtype != numpy.dtype(numpy.float32) or atom2.dtype != numpy.dtype(numpy.float32) 
        or atom3.dtype != numpy.dtype(numpy.float32) or atom4.dtype != numpy.dtype(numpy.float32) ):
        raise TypeError("all coordinates must be of type numpy.float32")

    with_PBC = (box is not None)
    if with_PBC:
        boxtype = boxCheck(box)
        if (boxtype == 'unknown'):
            raise ValueError("box input not recognised, must be an array of box dimensions")
        if (boxtype == 'tri_box'): # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))
        if (box.dtype!=numpy.dtype(numpy.float32)):
            raise TypeError("periodic boundaries must be of type float32")

    if not result is None:
        if (result.nd != 1 or result.dimensions[0] != numatom):
            raise ValueError("result array has incorrect size - should be (%d)"%(numatom))
        if (result.dtype != numpy.dtype(numpy.float64)):
            raise TypeError("result array must be of type numpy.float64")
        angles = numpy.asarray(result)
    else:
        angles = numpy.zeros((numatom,), numpy.float64)

    if with_PBC:
        if boxtype == 'ortho':
            calc_torsion_ortho(<coordinate*>atom1.data,<coordinate*>atom2.data,<coordinate*>atom3.data,<coordinate*>atom4.data,
                                numatom, <float*>box.data, <double*>angles.data)
        else:
            calc_torsion_triclinic(<coordinate*>atom1.data,<coordinate*>atom2.data,<coordinate*>atom3.data,<coordinate*>atom4.data,
                                    numatom, <coordinate*> box.data, <double*>angles.data)
    else:
        calc_torsion(<coordinate*>atom1.data,<coordinate*>atom2.data,<coordinate*>atom3.data,<coordinate*>atom4.data,
                      numatom,<double*>angles.data)

    return angles

def applyPBC(c_numpy.ndarray incoords, c_numpy.ndarray box):
    """Moves a set of coordinates to all be within the primary unit cell

    newcoords = applyPBC(coords, box)

    :Arguments:
                *coords*
                          coordinate array (of type numpy.float32)
                *box*
                          box dimensions, can be either orthogonal or triclinic information

    :Returns:
                *newcoords*
                          coordinates that are now all within the primary unit cell, as defined by box

    .. versionadded:: 0.8
    """
    cdef c_numpy.ndarray coords
    cdef c_numpy.ndarray box_inv
    cdef int coordnum

    coords = incoords.copy('C')

    # checks in input array
    if (coords.nd != 2 or coords.dimensions[1] != 3):
        raise ValueError("conf must be a sequence of 3 dimensional coordinates")
    if (coords.dtype != numpy.dtype(numpy.float32)):
        raise TypeError("coordinate data must be of type float32")

    coordnum = coords.dimensions[0]

    # determine boxtype
    boxtype = boxCheck(box)
    if (boxtype == 'unknown'):
        raise ValueError("box input not recognised, must be an array of box dimensions")
    if (boxtype == 'tri_box'): # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
        box = triclinic_vectors(box)
    if (boxtype == 'tri_vecs_bad'):
        box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))
    if (box.dtype!=numpy.dtype(numpy.float32)):
        raise TypeError("periodic boundaries must be of type float32")

    box_inv = numpy.zeros((3), dtype=numpy.float32)
    if boxtype == 'ortho':
        box_inv[0] = 1.0 / box[0]
        box_inv[1] = 1.0 / box[1]
        box_inv[2] = 1.0 / box[2]
        ortho_pbc(<coordinate*> coords.data, coordnum, <float*> box.data, <float*> box_inv.data)
    else:
        box_inv[0] = 1.0 / box[0][0]
        box_inv[1] = 1.0 / box[1][1]
        box_inv[2] = 1.0 / box[2][2]
        triclinic_pbc(<coordinate*> coords.data, coordnum, <coordinate*> box.data, <float*> box_inv.data)

    return coords
