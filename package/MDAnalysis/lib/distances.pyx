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
Distance calculation library --- :mod:`MDAnalysis.lib._distances`
=================================================================

Serial versions of all distance calculations
"""

cimport cython
cimport c_numpy
c_numpy.import_array()

cdef extern from "string.h":
    void* memcpy(void *dst, void *src, int len)

cdef extern from "calc_distances.h":
    ctypedef float coordinate[3]

    void _calc_distance_array(coordinate* ref, int numref, coordinate* conf, int numconf, double* distances)
    void _calc_distance_array_ortho(coordinate* ref, int numref, coordinate* conf, int numconf, float* box, double* distances)
    void _calc_distance_array_triclinic(coordinate* ref, int numref, coordinate* conf, int numconf, coordinate* box, double* distances)
    void _calc_self_distance_array(coordinate* ref, int numref, double* distances, int distnum)
    void _calc_self_distance_array_ortho(coordinate* ref, int numref, float* box, double* distances, int distnum)
    void _calc_self_distance_array_triclinic(coordinate* ref, int numref, coordinate* box, double* distances, int distnum)
    void _coord_transform(coordinate* coords, int numCoords, coordinate* box)
    void _calc_bond_distance(coordinate* atom1, coordinate* atom2, int numatom, double* distances)
    void _calc_bond_distance_ortho(coordinate* atom1, coordinate* atom2, int numatom, float*box, double* distances)
    void _calc_bond_distance_triclinic(coordinate* atom1, coordinate* atom2, int numatom, coordinate* box, double* distances)
    void _calc_angle(coordinate* atom1, coordinate* atom2, coordinate* atom3, int numatom, double* angles)
    void _calc_angle_ortho(coordinate* atom1, coordinate* atom2, coordinate* atom3, int numatom, float* box, double* angles)
    void _calc_angle_triclinic(coordinate* atom1, coordinate* atom2, coordinate* atom3, int numatom, coordinate* box, double* angles)
    void _calc_dihedral(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, int numatom, double* angles)
    void _calc_dihedral_ortho(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, int numatom, float* box, double* angles)
    void _calc_dihedral_triclinic(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, int numatom, coordinate* box, double* angles)
    void _ortho_pbc(coordinate* coords, int numcoords, float* box, float* box_inverse)
    void _triclinic_pbc(coordinate* coords, int numcoords, coordinate* box, float* box_inverse)
    void minimum_image(double *x, float *box, float *inverse_box)

def calc_distance_array(c_numpy.ndarray ref, c_numpy.ndarray conf,
                        c_numpy.ndarray result):
    cdef int confnum, refnum
    confnum = conf.dimensions[0]
    refnum = ref.dimensions[0]

    _calc_distance_array(<coordinate*>ref.data, refnum,
                         <coordinate*>conf.data, confnum,
                         <double*>result.data)

def calc_distance_array_ortho(c_numpy.ndarray ref, c_numpy.ndarray conf,
                              c_numpy.ndarray box,
                              c_numpy.ndarray result):
    cdef int confnum, refnum
    confnum = conf.dimensions[0]
    refnum = ref.dimensions[0]

    _calc_distance_array_ortho(<coordinate*>ref.data, refnum,
                               <coordinate*>conf.data, confnum,
                               <float*>box.data,
                               <double*>result.data)

def calc_distance_array_triclinic(c_numpy.ndarray ref, c_numpy.ndarray conf,
                                  c_numpy.ndarray box,
                                  c_numpy.ndarray result):
    cdef int confnum, refnum
    confnum = conf.dimensions[0]
    refnum = ref.dimensions[0]

    _calc_distance_array_triclinic(<coordinate*>ref.data, refnum,
                                   <coordinate*>conf.data, confnum,
                                   <coordinate*>box.data,
                                   <double*>result.data)

def calc_self_distance_array(c_numpy.ndarray ref,
                             c_numpy.ndarray result):
    cdef int refnum, distnum
    refnum = ref.dimensions[0]
    distnum = (refnum*(refnum-1))/2

    _calc_self_distance_array(<coordinate*>ref.data, refnum,
                              <double*>result.data, distnum)

def calc_self_distance_array_ortho(c_numpy.ndarray ref,
                                   c_numpy.ndarray box,
                                   c_numpy.ndarray result):
    cdef int refnum, distnum
    refnum = ref.dimensions[0]
    distnum = (refnum*(refnum-1))/2

    _calc_self_distance_array_ortho(<coordinate*>ref.data, refnum,
                                    <float*>box.data,
                                    <double*>result.data, distnum)

def calc_self_distance_array_triclinic(c_numpy.ndarray ref,
                                       c_numpy.ndarray box,
                                       c_numpy.ndarray result):
    cdef int refnum, distnum
    refnum = ref.dimensions[0]
    distnum = (refnum*(refnum-1))/2

    _calc_self_distance_array_triclinic(<coordinate*>ref.data, refnum,
                                        <coordinate*>box.data,
                                        <double*>result.data, distnum)

def coord_transform(c_numpy.ndarray coords,
                    c_numpy.ndarray box):
    cdef int numcoords
    numcoords = coords.dimensions[0]

    _coord_transform(<coordinate*> coords.data, numcoords,
                     <coordinate*> box.data)

def calc_bond_distance(c_numpy.ndarray coords1,
                       c_numpy.ndarray coords2,
                       c_numpy.ndarray results):
    cdef int numcoords
    numcoords = coords1.dimensions[0]

    _calc_bond_distance(<coordinate*> coords1.data, <coordinate*> coords2.data,
                        numcoords,
                        <double*>results.data)

def calc_bond_distance_ortho(c_numpy.ndarray coords1,
                             c_numpy.ndarray coords2,
                             c_numpy.ndarray box,
                             c_numpy.ndarray results):
    cdef int numcoords
    numcoords = coords1.dimensions[0]

    _calc_bond_distance_ortho(<coordinate*> coords1.data, <coordinate*> coords2.data,
                              numcoords,
                              <float*>box.data,
                              <double*>results.data)

def calc_bond_distance_triclinic(c_numpy.ndarray coords1,
                                 c_numpy.ndarray coords2,
                                 c_numpy.ndarray box,
                                 c_numpy.ndarray results):
    cdef int numcoords
    numcoords = coords1.dimensions[0]

    _calc_bond_distance_triclinic(<coordinate*> coords1.data, <coordinate*> coords2.data,
                                  numcoords,
                                  <coordinate*>box.data,
                                  <double*>results.data)

def calc_angle(c_numpy.ndarray coords1,
               c_numpy.ndarray coords2,
               c_numpy.ndarray coords3,
               c_numpy.ndarray results):
    cdef int numcoords
    numcoords = coords1.dimensions[0]

    _calc_angle(<coordinate*> coords1.data, <coordinate*> coords2.data,
                <coordinate*> coords3.data,
                numcoords,
                <double*>results.data)

def calc_angle_ortho(c_numpy.ndarray coords1,
                     c_numpy.ndarray coords2,
                     c_numpy.ndarray coords3,
                     c_numpy.ndarray box,
                     c_numpy.ndarray results):
    cdef int numcoords
    numcoords = coords1.dimensions[0]

    _calc_angle_ortho(<coordinate*> coords1.data, <coordinate*> coords2.data,
                      <coordinate*> coords3.data,
                      numcoords,
                      <float*>box.data,
                      <double*>results.data)

def calc_angle_triclinic(c_numpy.ndarray coords1,
                         c_numpy.ndarray coords2,
                         c_numpy.ndarray coords3,
                         c_numpy.ndarray box,
                         c_numpy.ndarray results):
    cdef int numcoords
    numcoords = coords1.dimensions[0]

    _calc_angle_triclinic(<coordinate*> coords1.data, <coordinate*> coords2.data,
                          <coordinate*> coords3.data,
                          numcoords,
                          <coordinate*>box.data,
                          <double*>results.data)

def calc_dihedral(c_numpy.ndarray coords1,
                 c_numpy.ndarray coords2,
                 c_numpy.ndarray coords3,
                 c_numpy.ndarray coords4,
                 c_numpy.ndarray results):
    cdef int numcoords
    numcoords = coords1.dimensions[0]

    _calc_dihedral(<coordinate*> coords1.data, <coordinate*> coords2.data,
                  <coordinate*> coords3.data, <coordinate*> coords4.data,
                  numcoords,
                  <double*>results.data)

def calc_dihedral_ortho(c_numpy.ndarray coords1,
                       c_numpy.ndarray coords2,
                       c_numpy.ndarray coords3,
                       c_numpy.ndarray coords4,
                       c_numpy.ndarray box,
                       c_numpy.ndarray results):
    cdef int numcoords
    numcoords = coords1.dimensions[0]

    _calc_dihedral_ortho(<coordinate*> coords1.data, <coordinate*> coords2.data,
                        <coordinate*> coords3.data, <coordinate*> coords4.data,
                        numcoords,
                        <float*>box.data,
                        <double*>results.data)

def calc_dihedral_triclinic(c_numpy.ndarray coords1,
                           c_numpy.ndarray coords2,
                           c_numpy.ndarray coords3,
                           c_numpy.ndarray coords4,
                           c_numpy.ndarray box,
                           c_numpy.ndarray results):
    cdef int numcoords
    numcoords = coords1.dimensions[0]

    _calc_dihedral_triclinic(<coordinate*> coords1.data, <coordinate*> coords2.data,
                            <coordinate*> coords3.data, <coordinate*> coords4.data,
                            numcoords,
                            <coordinate*>box.data,
                            <double*>results.data)

def ortho_pbc(c_numpy.ndarray coords,
              c_numpy.ndarray box, c_numpy.ndarray box_inverse):
    cdef int numcoords
    numcoords = coords.dimensions[0]

    _ortho_pbc(<coordinate*> coords.data, numcoords,
               <float*>box.data, <float*>box_inverse.data)

def triclinic_pbc(c_numpy.ndarray coords,
                  c_numpy.ndarray box, c_numpy.ndarray box_inverse):
    cdef int numcoords
    numcoords = coords.dimensions[0]

    _triclinic_pbc(<coordinate*> coords.data, numcoords,
                   <coordinate*> box.data, <float*>box_inverse.data)


@cython.boundscheck(False)
def contact_matrix_no_pbc(coord, sparse_contacts, cutoff):
    cdef int rows = len(coord)
    cdef double cutoff2 = cutoff ** 2
    cdef float[:, ::1] coord_view = coord

    cdef int i, j
    cdef double[3] rr;
    cdef double dist
    for i in range(rows):
        sparse_contacts[i, i] = True
        for j in range(i+1, rows):
            rr[0] = coord_view[i, 0] - coord_view[j, 0]
            rr[1] = coord_view[i, 1] - coord_view[j, 1]
            rr[2] = coord_view[i, 2] - coord_view[j, 2]
            dist = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2]
            if dist < cutoff2:
                sparse_contacts[i, j] = True
                sparse_contacts[j, i] = True


@cython.boundscheck(False)
def contact_matrix_pbc(coord, sparse_contacts, box, cutoff):
    cdef int rows = len(coord)
    cdef double cutoff2 = cutoff ** 2
    cdef float[:, ::1] coord_view = coord
    cdef float[::1] box_view = box
    cdef float[::1] box_inv = 1. / box

    cdef int i, j
    cdef double[3] rr;
    cdef double dist
    for i in range(rows):
        sparse_contacts[i, i] = True
        for j in range(i+1, rows):
            rr[0] = coord_view[i, 0] - coord_view[j, 0]
            rr[1] = coord_view[i, 1] - coord_view[j, 1]
            rr[2] = coord_view[i, 2] - coord_view[j, 2]

            minimum_image(rr, &box_view[0], &box_inv[0])

            dist = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2]

            if dist < cutoff2:
                sparse_contacts[i, j] = True
                sparse_contacts[j, i] = True
