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
from .mdamath import triclinic_vectors
cimport numpy as np
cimport _cutil
from _cutil cimport _dot ,_norm, _cross

from libcpp.vector cimport vector


__all__ = ['augment_coordinates', 'undo_augment']


@cython.boundscheck(False)
@cython.wraparound(False)
def augment_coordinates(float[:, ::1] coordinates, float[:] box, float r):
    r"""Calculates the relevant images of particles which are within a
    distance 'r' from the box walls

    The algorithm works by generating explicit periodic images of
    interior atoms residing close to any of the six box walls. 
    The steps involved in generating images involves
    evaluation of reciprocal vectors for the given box vectors 
    followed by calculation of projection distance of atom along the 
    reciprocal vectors. If the distance is less than a
    specified cutoff distance, relevant periodic images are generated
    using box translation vectors i.e. ``l[a] + m[b] + n[c]``, where 
    ``[l, m, n]`` are the neighbouring cell indices relative to the central cell, 
    and ``[a, b, c]`` are the box vectors. For instance, an atom close to
    ``XY`` plane containing origin will generate a periodic image
    outside the central cell and close to the opposite `XY` plane
    of the box i.e. at ``0[a] + 0[b] + 1[c]``. 
    Similarly, if the particle is close to more than
    one box walls, images along the diagonals are also generated ::


                           |  x               x
        +---------------+  |    +---------------+
        |               |  |    |               |
        |               |  |    |               |
        |               |  |    |               |
        |             o |  |  x |             o |
        +---------------+  |    +---------------+
                           |



    Parameters
    ----------
    coordinates : array
      Input coordinate array to generate duplicate images
      in the vicinity of the central cell. All the coordinates
      must be within the primary unit cell. (dtype = numpy.float32)
    box : array
      Box dimension of shape (6, ). The dimensions must be
      provided in the same format as returned
      by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
      ``[lx, ly, lz, alpha, beta, gamma]`` (dtype = numpy.float32)
    r : float
      thickness of cutoff region for duplicate image generation

    Returns
    -------
    output : array
      coordinates of duplicate(augmented) particles (dtype = numpy.float32)
    indices : array
      original indices of the augmented coordinates (dtype = numpy.int64)
      A map which translates the indices of augmented particles
      to their original particle index such that
      ``indices[augmentedindex] = originalindex``

    Note
    ----
    Output doesnot return coordinates from the initial array.
    To merge the particles with their respective images, following operation
    needs to be superseded after generating the images:
    
    .. code-block:: python

            images, mapping = augment_coordinates(coordinates, box, max_cutoff)
            all_coords = np.concatenate([coordinates, images])


    See Also
    --------
    MDAnalysis.lib._augment.undo_augment


    .. versionadded:: 0.19.0
    """
    cdef bint lo_x, hi_x, lo_y, hi_y, lo_z, hi_z
    cdef int i, j, N
    cdef float norm
    cdef float shiftX[3]
    cdef float shiftY[3]
    cdef float shiftZ[3]
    cdef float coord[3]
    cdef float end[3]
    cdef float other[3]
    cdef float dm[3][3]
    cdef float reciprocal[3][3]

    dm = triclinic_vectors(box)

    for i in range(3):
        shiftX[i] = dm[0][i]
        shiftY[i] = dm[1][i]
        shiftZ[i] = dm[2][i]
        end[i] = dm[0][i] + dm[1][i] + dm[2][i]
    # Calculate reciprocal vectors
    _cross(&dm[1][0], &dm[2][0], &reciprocal[0][0])
    _cross(&dm[2][0], &dm[0][0], &reciprocal[1][0])
    _cross(&dm[0][0], &dm[1][0], &reciprocal[2][0])
    # Normalize
    for i in range(3):
        norm = _norm(&reciprocal[i][0])
        for j in range(3):
            reciprocal[i][j] = reciprocal[i][j]/norm

    N = coordinates.shape[0]

    cdef vector[float] output
    cdef vector[int] indices

    for i in range(N):
        for j in range(3):
            coord[j] = coordinates[i, j]
            other[j] = end[j] - coordinates[i, j]
        # identify the condition
        lo_x = _dot(&coord[0], &reciprocal[0][0]) <= r
        hi_x = _dot(&other[0], &reciprocal[0][0]) <= r
        lo_y = _dot(&coord[0], &reciprocal[1][0]) <= r
        hi_y = _dot(&other[0], &reciprocal[1][0]) <= r
        lo_z = _dot(&coord[0], &reciprocal[2][0]) <= r
        hi_z = _dot(&other[0], &reciprocal[2][0]) <= r

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
    return np.asarray(output, dtype=np.float32).reshape(n, 3), np.asarray(indices, dtype=np.int64)


@cython.boundscheck(False)
@cython.wraparound(False)
def undo_augment(np.int64_t[:] results, np.int64_t[:] translation, int nreal):
    """Translate augmented indices back to original indices

    Parameters
    ----------
    results : numpy.ndarray
      indices of coordinates, including "augmented" indices (dtype = numpy.int64)
    translation : numpy.ndarray
      Map to link the augmented indices to the original particle indices
      such that ``translation[augmentedindex] = originalindex``
      (dtype = numpy.int64)
    nreal : int
      number of real coordinates, i.e. values in results equal or larger
      than this need to be translated to their real counterpart


    Returns
    -------
    results : numpy.ndarray
      modified input results with all the augmented indices
      translated to their corresponding initial original indices
      (dtype = numpy.int64)

    Note
    ----
    Modifies the results array in place

    See Also
    --------
    'MDAnalysis.lib._augment.augment_coordinates'


    .. versionadded:: 0.19.0
    """
    cdef int N
    cdef ssize_t i
    N = results.shape[0]

    for i in range(N):
        if results[i] >= nreal:
            results[i] = translation[results[i] - nreal]
    return np.asarray(results, dtype=np.int64)
