# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
Parallel distance calculation library --- :mod:`MDAnalysis.lib.c_distances_openmp`
==================================================================================


Contains OpenMP versions of the contents of "calc_distances.h"
"""

from libc.stdint cimport uint64_t
import numpy
cimport numpy
numpy.import_array()

cdef extern from "string.h":
    void* memcpy(void* dst, void* src, int len)

cdef extern from "calc_distances.h":
    ctypedef float coordinate[3]
    cdef bint USED_OPENMP
    void _calc_distance_array(coordinate* ref, uint64_t numref, coordinate* conf, uint64_t numconf, double* distances)
    void _calc_distance_array_ortho(coordinate* ref, uint64_t numref, coordinate* conf, uint64_t numconf, float* box, double* distances)
    void _calc_distance_array_triclinic(coordinate* ref, uint64_t numref, coordinate* conf, uint64_t numconf, float* box, double* distances)
    void _calc_self_distance_array(coordinate* ref, uint64_t numref, double* distances)
    void _calc_self_distance_array_ortho(coordinate* ref, uint64_t numref, float* box, double* distances)
    void _calc_self_distance_array_triclinic(coordinate* ref, uint64_t numref, float* box, double* distances)
    void _coord_transform(coordinate* coords, uint64_t numCoords, double* box)
    void _calc_bond_distance(coordinate* atom1, coordinate* atom2, uint64_t numatom, double* distances)
    void _calc_bond_distance_ortho(coordinate* atom1, coordinate* atom2, uint64_t numatom, float* box, double* distances)
    void _calc_bond_distance_triclinic(coordinate* atom1, coordinate* atom2, uint64_t numatom, float* box, double* distances)
    void _calc_angle(coordinate* atom1, coordinate* atom2, coordinate* atom3, uint64_t numatom, double* angles)
    void _calc_angle_ortho(coordinate* atom1, coordinate* atom2, coordinate* atom3, uint64_t numatom, float* box, double* angles)
    void _calc_angle_triclinic(coordinate* atom1, coordinate* atom2, coordinate* atom3, uint64_t numatom, float* box, double* angles)
    void _calc_dihedral(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, uint64_t numatom, double* angles)
    void _calc_dihedral_ortho(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, uint64_t numatom, float* box, double* angles)
    void _calc_dihedral_triclinic(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4, uint64_t numatom, float* box, double* angles)
    void _ortho_pbc(coordinate* coords, uint64_t numcoords, float* box)
    void _triclinic_pbc(coordinate* coords, uint64_t numcoords, float* box)


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


def calc_bond_distance_ortho(numpy.ndarray coords1,
                             numpy.ndarray coords2,
                             numpy.ndarray box,
                             numpy.ndarray results):
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
