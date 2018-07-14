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
import numpy
cimport numpy

from libcpp.vector cimport vector
from libc.math cimport sqrt

__all__ = ['augment', 'undo_augment']

@cython.boundscheck(False)
@cython.wraparound(False)
def augment(float[:, ::1] coordinates, float[:, ::1] dm, float r):
    """Calculate augmented coordinate set

    Parameters
    ----------
    coordinates : array
      Input coordinate array to generate duplicate images
    dm : array
      Real space box vectors in a matrix with shape (3, 3)
      Vectors `[[a], [b], [c]]` are box vectors
    r : float
      thickness of cutoff region for duplicate image generation

    Returns
    -------
    output : array
      coordinates of duplicates generated due to periodic boundary conditions
    indices : array
      original indices of the augmented coordinates
    """
    cdef bint lo_x, hi_x, lo_y, hi_y, lo_z, hi_z
    cdef int i, j, p, N
    cdef float norm
    cdef float shiftX[3]
    cdef float shiftY[3]
    cdef float shiftZ[3]
    cdef float coord[3]
    cdef float end[3]
    cdef float other[3]
    cdef float[:, ::1] reciprocal = numpy.zeros((3, 3), dtype=numpy.float32)
    for i in range(3):
        shiftX[i] = dm[0, i]
        shiftY[i] = dm[1, i]
        shiftZ[i] = dm[2, i]
        end[i] = dm[0, i] + dm[1, i] + dm[2, i]
    # Calculate reciprocal vectors
    _cross(&dm[1, 0], &dm[2, 0], &reciprocal[0, 0])
    _cross(&dm[2, 0], &dm[0, 0], &reciprocal[1, 0])
    _cross(&dm[0, 0], &dm[1, 0], &reciprocal[2, 0])
    # Normalize
    for i in range(3):
        norm = _norm(&reciprocal[i, 0])
        for j in range(3):
            reciprocal[i, j] = reciprocal[i, j]/norm

    N = coordinates.shape[0]

    cdef vector[float] output
    cdef vector[int] indices

    for i in range(0, N):
        for j in range(3):
            coord[j] = coordinates[i, j]
            other[j] = end[j] - coordinates[i, j]
        # identify the condition
        lo_x = _dot(&coord[0], &reciprocal[0, 0]) <= r
        hi_x = _dot(&other[0], &reciprocal[0, 0]) <= r
        lo_y = _dot(&coord[0], &reciprocal[1, 0]) <= r
        hi_y = _dot(&other[0], &reciprocal[1, 0]) <= r
        lo_z = _dot(&coord[0], &reciprocal[2, 0]) <= r
        hi_z = _dot(&other[0], &reciprocal[2, 0]) <= r

        if lo_x:
            # if X, face piece
            for j in range(3):
                # add to output
                output.push_back(coord[j] + shiftX[j])
            # keep record of which index this augmented
            # position was created from
            indices.push_back(i)

            if lo_y:
                # if X&Y, edge piece
                for j in range(3):
                    output.push_back(coord[j] + shiftX[j] + shiftY[j])
                indices.push_back(i)

                if lo_z:
                    # if X&Y&Z, corner piece
                    for j in range(3):
                        output.push_back(coord[j] + shiftX[j] + shiftY[j] + shiftZ[j])
                    indices.push_back(i)

                elif hi_z:
                    for j in range(3):
                        output.push_back(coord[j] + shiftX[j] + shiftY[j] - shiftZ[j])
                    indices.push_back(i)

            elif hi_y:
                for j in range(3):
                    output.push_back(coord[j] + shiftX[j] - shiftY[j])
                indices.push_back(i)

                if lo_z:
                    for j in range(3):
                        output.push_back(coord[j] + shiftX[j] - shiftY[j] + shiftZ[j])
                    indices.push_back(i)

                elif hi_z:
                    for j in range(3):
                        output.push_back(coord[j] + shiftX[j] - shiftY[j] - shiftZ[j])
                    indices.push_back(i)

            if lo_z:
                for j in range(3):
                    output.push_back(coord[j] + shiftX[j] + shiftZ[j])
                indices.push_back(i)

            elif hi_z:
                for j in range(3):
                    output.push_back(coord[j] + shiftX[j] - shiftZ[j])
                indices.push_back(i)

        elif hi_x:
            for j in range(3):
                output.push_back(coord[j] - shiftX[j])
            indices.push_back(i)

            if lo_y:
                for j in range(3):
                    output.push_back(coord[j] - shiftX[j] + shiftY[j])
                indices.push_back(i)

                if lo_z:
                    for j in range(3):
                        output.push_back(coord[j] - shiftX[j] + shiftY[j] + shiftZ[j])
                    indices.push_back(i)

                elif hi_z:
                    for j in range(3):
                        output.push_back(coord[j] - shiftX[j] + shiftY[j] - shiftZ[j])
                    indices.push_back(i)

            elif hi_y:
                for j in range(3):
                    output.push_back(coord[j] - shiftX[j] - shiftY[j])
                indices.push_back(i)

                if lo_z:
                    for j in range(3):
                        output.push_back(coord[j] - shiftX[j] - shiftY[j] + shiftZ[j])
                    indices.push_back(i)

                elif hi_z:
                    for j in range(3):
                        output.push_back(coord[j] - shiftX[j] - shiftY[j] - shiftZ[j])
                    indices.push_back(i)

            if lo_z:
                for j in range(3):
                    output.push_back(coord[j] - shiftX[j] + shiftZ[j])
                indices.push_back(i)

            elif hi_z:
                for j in range(3):
                    output.push_back(coord[j] - shiftX[j] - shiftZ[j])
                indices.push_back(i)

        if lo_y:
            for j in range(3):
                output.push_back(coord[j] + shiftY[j])
            indices.push_back(i)

            if lo_z:
                for j in range(3):
                    output.push_back(coord[j] + shiftY[j] + shiftZ[j])
                indices.push_back(i)

            elif hi_z:
                for j in range(3):
                    output.push_back(coord[j] + shiftY[j] - shiftZ[j])
                indices.push_back(i)

        elif hi_y:
            for j in range(3):
                output.push_back(coord[j] - shiftY[j])
            indices.push_back(i)

            if lo_z:
                for j in range(3):
                    output.push_back(coord[j] - shiftY[j] + shiftZ[j])
                indices.push_back(i)

            elif hi_z:
                for j in range(3):
                    output.push_back(coord[j] - shiftY[j] - shiftZ[j])
                indices.push_back(i)

        if lo_z:
            for j in range(3):
                output.push_back(coord[j] + shiftZ[j])
            indices.push_back(i)

        elif hi_z:
            for j in range(3):
                output.push_back(coord[j] - shiftZ[j])
            indices.push_back(i)
    n = indices.size()
    return numpy.asarray(output, dtype=numpy.float32).reshape(n, 3), numpy.asarray(indices, dtype=numpy.int32)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef float _dot(float * a, float * b):
    """Return dot product of two sequences in range."""
    cdef ssize_t n
    cdef float sum1

    sum1 = 0.0
    for n in range(3):
        sum1 += a[n] * b[n]
    return sum1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _cross(float * a, float * b, float * result):
    """
    Calculates the cross product between vectors
    given by pointers a and b

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


@cython.boundscheck(False)
@cython.wraparound(False)
def undo_augment(int[:] results, int[:] translation, int nreal):
    """Translate augmented indices back to originals

    Parameters
    ----------
    results : ndarray of ints
      indices of coordinates, including "augmented" indices
    translation : ndarray of ints
      original indices of augmented coordinates
    nreal : int
      number of real coordinates, ie values in results equal or larger
      than this need to be translated to their real counterpart

    Returns
    -------
    results : ndarray of ints

    Note
    ----
    Modifies the results array

    """
    cdef int N
    N = results.shape[0]

    for i in range(N):
        if results[i] >= nreal:
            results[i] = translation[results[i] - nreal]
    return results
