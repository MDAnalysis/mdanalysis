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

"""
Distance calculation library --- :mod:`MDAnalysis.lib.c_distances`
==================================================================

Serial versions of all distance calculations
"""

cimport cython
from libc.stdint cimport uint64_t, UINT64_MAX
import numpy
cimport numpy
numpy.import_array()

from libc.math cimport round as cround
from libc.float cimport FLT_MAX, DBL_MAX

# make UINT64_MAX visible at the python layer
_UINT64_MAX = UINT64_MAX


OPENMP_ENABLED = True if USED_OPENMP else False


def calc_distance_array(numpy.ndarray ref, numpy.ndarray conf,
                        numpy.ndarray result):
    cdef uint64_t confnum, refnum
    confnum = conf.shape[0]
    refnum = ref.shape[0]

    _calc_distance_array(<coordinate*> ref.data, refnum,
                         <coordinate*> conf.data, confnum,
                         <double*> result.data)


def calc_distance_array_ortho(numpy.ndarray ref, numpy.ndarray conf,
                              numpy.ndarray box, numpy.ndarray result):
    cdef uint64_t confnum, refnum
    confnum = conf.shape[0]
    refnum = ref.shape[0]

    _calc_distance_array_ortho(<coordinate*> ref.data, refnum,
                               <coordinate*> conf.data, confnum,
                               <float*> box.data, <double*> result.data)


def calc_distance_array_triclinic(numpy.ndarray ref, numpy.ndarray conf,
                                  numpy.ndarray box, numpy.ndarray result):
    cdef uint64_t confnum, refnum
    confnum = conf.shape[0]
    refnum = ref.shape[0]

    _calc_distance_array_triclinic(<coordinate*> ref.data, refnum,
                                   <coordinate*> conf.data, confnum,
                                   <float*> box.data, <double*> result.data)


def calc_self_distance_array(numpy.ndarray ref, numpy.ndarray result):
    cdef uint64_t refnum
    refnum = ref.shape[0]

    _calc_self_distance_array(<coordinate*> ref.data, refnum,
                              <double*> result.data)


def calc_self_distance_array_ortho(numpy.ndarray ref, numpy.ndarray box,
                                   numpy.ndarray result):
    cdef uint64_t refnum
    refnum = ref.shape[0]

    _calc_self_distance_array_ortho(<coordinate*> ref.data, refnum,
                                    <float*> box.data, <double*> result.data)


def calc_self_distance_array_triclinic(numpy.ndarray ref, numpy.ndarray box,
                                       numpy.ndarray result):
    cdef uint64_t refnum
    refnum = ref.shape[0]

    _calc_self_distance_array_triclinic(<coordinate*> ref.data, refnum,
                                        <float*> box.data,
                                        <double*> result.data)


def coord_transform(numpy.ndarray coords, numpy.ndarray box):
    cdef uint64_t numcoords
    numcoords = coords.shape[0]

    _coord_transform(<coordinate*> coords.data, numcoords, <double*> box.data)


def calc_bond_distance(numpy.ndarray coords1, numpy.ndarray coords2,
                       numpy.ndarray results):
    cdef uint64_t numcoords
    numcoords = coords1.shape[0]

    _calc_bond_distance(<coordinate*> coords1.data, <coordinate*> coords2.data,
                        numcoords, <double*> results.data)


def calc_bond_distance_ortho(numpy.ndarray coords1, numpy.ndarray coords2,
                             numpy.ndarray box, numpy.ndarray results):
    cdef uint64_t numcoords
    numcoords = coords1.shape[0]

    _calc_bond_distance_ortho(<coordinate*> coords1.data,
                              <coordinate*> coords2.data, numcoords,
                              <float*> box.data, <double*> results.data)


def calc_bond_distance_triclinic(numpy.ndarray coords1, numpy.ndarray coords2,
                                 numpy.ndarray box, numpy.ndarray results):
    cdef uint64_t numcoords
    numcoords = coords1.shape[0]

    _calc_bond_distance_triclinic(<coordinate*> coords1.data,
                                  <coordinate*> coords2.data, numcoords,
                                  <float*> box.data, <double*> results.data)


def calc_angle(numpy.ndarray coords1, numpy.ndarray coords2,
               numpy.ndarray coords3, numpy.ndarray results):
    cdef uint64_t numcoords
    numcoords = coords1.shape[0]

    _calc_angle(<coordinate*> coords1.data, <coordinate*> coords2.data,
                <coordinate*> coords3.data, numcoords, <double*> results.data)


def calc_angle_ortho(numpy.ndarray coords1, numpy.ndarray coords2,
                     numpy.ndarray coords3, numpy.ndarray box,
                     numpy.ndarray results):
    cdef uint64_t numcoords
    numcoords = coords1.shape[0]

    _calc_angle_ortho(<coordinate*> coords1.data, <coordinate*> coords2.data,
                      <coordinate*> coords3.data, numcoords, <float*> box.data,
                      <double*> results.data)


def calc_angle_triclinic(numpy.ndarray coords1, numpy.ndarray coords2,
                         numpy.ndarray coords3, numpy.ndarray box,
                         numpy.ndarray results):
    cdef uint64_t numcoords
    numcoords = coords1.shape[0]

    _calc_angle_triclinic(<coordinate*> coords1.data,
                          <coordinate*> coords2.data,
                          <coordinate*> coords3.data, numcoords,
                          <float*> box.data, <double*> results.data)


def calc_dihedral(numpy.ndarray coords1, numpy.ndarray coords2,
                  numpy.ndarray coords3, numpy.ndarray coords4,
                  numpy.ndarray results):
    cdef uint64_t numcoords
    numcoords = coords1.shape[0]

    _calc_dihedral(<coordinate*> coords1.data, <coordinate*> coords2.data,
                   <coordinate*> coords3.data, <coordinate*> coords4.data,
                   numcoords, <double*> results.data)


def calc_dihedral_ortho(numpy.ndarray coords1, numpy.ndarray coords2,
                        numpy.ndarray coords3, numpy.ndarray coords4,
                        numpy.ndarray box, numpy.ndarray results):
    cdef uint64_t numcoords
    numcoords = coords1.shape[0]

    _calc_dihedral_ortho(<coordinate*> coords1.data, <coordinate*> coords2.data,
                         <coordinate*> coords3.data, <coordinate*> coords4.data,
                         numcoords, <float*> box.data, <double*> results.data)


def calc_dihedral_triclinic(numpy.ndarray coords1, numpy.ndarray coords2,
                            numpy.ndarray coords3, numpy.ndarray coords4,
                            numpy.ndarray box, numpy.ndarray results):
    cdef uint64_t numcoords
    numcoords = coords1.shape[0]

    _calc_dihedral_triclinic(<coordinate*> coords1.data,
                             <coordinate*> coords2.data,
                             <coordinate*> coords3.data,
                             <coordinate*> coords4.data, numcoords,
                             <float*> box.data, <double*> results.data)


def ortho_pbc(numpy.ndarray coords, numpy.ndarray box):
    cdef uint64_t numcoords
    numcoords = coords.shape[0]

    _ortho_pbc(<coordinate*> coords.data, numcoords, <float*> box.data)


def triclinic_pbc(numpy.ndarray coords, numpy.ndarray box):
    cdef uint64_t numcoords
    numcoords = coords.shape[0]

    _triclinic_pbc(<coordinate*> coords.data, numcoords, <float*> box.data)


@cython.boundscheck(False)
def contact_matrix_no_pbc(coord, sparse_contacts, cutoff):
    cdef int rows = len(coord)
    cdef double cutoff2 = cutoff ** 2
    cdef float[:, ::1] coord_view = coord

    cdef int i, j
    cdef double[3] rr
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
    cdef double[3] rr
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


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void _minimum_image_orthogonal(cython.floating[:] dx,
                                           cython.floating[:] box,
                                           cython.floating[:] inverse_box) nogil:
    """Minimize dx to be the shortest vector

    Parameters
    ----------
    dx : numpy.array, shape (3,)
      vector to minimize
    box : numpy.array, shape (3,)
      box length in each dimension
    inverse_box : numpy.array, shape (3,)
      inverse of box

    Operates in-place on dx!
    """
    cdef int i
    cdef cython.floating s

    for i in range(3):
        if box[i] > 0:
            s = inverse_box[i] * dx[i]
            dx[i] = box[i] * (s - cround(s))


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void _minimum_image_triclinic(cython.floating[:] dx,
                                          cython.floating[:] box,
                                          cython.floating[:] inverse_box) nogil:
    """Minimise dx to be the shortest vector

    Parameters
    ----------
    dx : numpy.array, shape (3,)
      the vector to apply minimum image convention to
    box : numpy.array, shape (9,)
      flattened 3x3 representation of the unit cell
    inverse_box : numpy.array, shape (3,)
      inverse of the **diagonal** of the 3x3 representation

    This version is near identical to the version in calc_distances.h, with the
    difference being that this version does not assume that coordinates are at
    most a single box length apart.

    Operates in-place on dx!
    """
    cdef cython.floating dx_min[3]
    cdef cython.floating s, dsq, dsq_min, rx
    cdef cython.floating ry[2]
    cdef cython.floating rz[3]
    cdef int ix, iy, iz

    # first make shift only 1 cell in any direction
    s = cround(inverse_box[2] * dx[2])
    dx[0] -= s * box[6]
    dx[1] -= s * box[7]
    dx[2] -= s * box[8]
    s = cround(inverse_box[1] * dx[1])
    dx[0] -= s * box[3]
    dx[1] -= s * box[4]
    s = cround(inverse_box[0] * dx[0])
    dx[0] -= s * box[0]

    if cython.floating is float:
        dsq_min = FLT_MAX
    else:
        dsq_min = DBL_MAX
    dx_min[0] = 0.0
    dx_min[1] = 0.0
    dx_min[2] = 0.0

    # then check all images to see which combination of 1 cell shifts gives the best shift
    for ix in range(-1, 2):
        rx = dx[0] + box[0] * ix
        for iy in range(-1, 2):
            ry[0] = rx + box[3] * iy
            ry[1] = dx[1] + box[4] * iy
            for iz in range(-1, 2):
                rz[0] = ry[0] + box[6] * iz
                rz[1] = ry[1] + box[7] * iz
                rz[2] = dx[2] + box[8] * iz
                dsq = rz[0] * rz[0] + rz[1] * rz[1] + rz[2] * rz[2]
                if (dsq < dsq_min):
                    dsq_min = dsq
                    dx_min[0] = rz[0]
                    dx_min[1] = rz[1]
                    dx_min[2] = rz[2]

    dx[0] = dx_min[0]
    dx[1] = dx_min[1]
    dx[2] = dx_min[2]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _minimize_vectors_ortho(cython.floating[:, :] vectors not None, cython.floating[:] box not None,
                            cython.floating[:, :] output not None):
    cdef int i, n
    cdef cython.floating box_inverse[3]
    cdef cython.floating[:] box_inverse_view

    box_inverse[0] = 1.0 / box[0]
    box_inverse[1] = 1.0 / box[1]
    box_inverse[2] = 1.0 / box[2]

    box_inverse_view = box_inverse

    n = len(vectors)
    with nogil:
        for i in range(n):
            output[i, 0] = vectors[i, 0]
            output[i, 1] = vectors[i, 1]
            output[i, 2] = vectors[i, 2]
            _minimum_image_orthogonal(output[i, :], box, box_inverse_view)


@cython.boundscheck(False)
@cython.wraparound(False)
def _minimize_vectors_triclinic(cython.floating[:, :] vectors not None, cython.floating[:] box not None,
                                cython.floating[:, :] output not None):
    cdef int i, n
    cdef cython.floating box_inverse[3]
    cdef cython.floating[:] box_inverse_view

    # this is only inverse of diagonal, used for initial shift to ensure vector is within a single shift away
    box_inverse[0] = 1.0 / box[0]
    box_inverse[1] = 1.0 / box[4]
    box_inverse[2] = 1.0 / box[8]

    box_inverse_view = box_inverse

    n = len(vectors)
    with nogil:
        for i in range(n):
            output[i, 0] = vectors[i, 0]
            output[i, 1] = vectors[i, 1]
            output[i, 2] = vectors[i, 2]
            _minimum_image_triclinic(output[i, :], box, box_inverse_view)
