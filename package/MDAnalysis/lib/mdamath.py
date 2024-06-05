# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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

"""
Mathematical helper functions --- :mod:`MDAnalysis.lib.mdamath`
===============================================================


Helper functions for common mathematical operations. Some of these functions
are written in C/cython for higher performance.

Linear algebra
--------------

.. autofunction:: normal
.. autofunction:: norm
.. autofunction:: pdot
.. autofunction:: pnorm
.. autofunction:: angle
.. autofunction:: dihedral
.. autofunction:: stp
.. autofunction:: sarrus_det
.. autofunction:: triclinic_box
.. autofunction:: triclinic_vectors
.. autofunction:: box_volume


Connectivity
------------

.. autofunction:: make_whole
.. autofunction:: find_fragments


.. versionadded:: 0.11.0
.. versionchanged: 1.0.0
   Unused function :func:`_angle()` has now been removed.

"""
import numpy as np

from ..exceptions import NoDataError
from . import util
from ._cutil import (make_whole, find_fragments, _sarrus_det_single,
                     _sarrus_det_multiple)
import numpy.typing as npt
from typing import Union

# geometric functions


def norm(v: npt.ArrayLike) -> float:
    r"""Calculate the norm of a vector v.

    .. math:: v = \sqrt{\mathbf{v}\cdot\mathbf{v}}

    This version is faster then numpy.linalg.norm because it only works for a
    single vector and therefore can skip a lot of the additional fuss
    linalg.norm does.

    Parameters
    ----------
    v : array_like
        1D array of shape (N) for a vector of length N

    Returns
    -------
    float
        norm of the vector

    """
    return np.sqrt(np.dot(v, v))


def normal(vec1: npt.ArrayLike, vec2: npt.ArrayLike) -> npt.NDArray:
    r"""Returns the unit vector normal to two vectors.

    .. math::

       \hat{\mathbf{n}} = \frac{\mathbf{v}_1 \times \mathbf{v}_2}{|\mathbf{v}_1 \times \mathbf{v}_2|}

    If the two vectors are collinear, the vector :math:`\mathbf{0}` is returned.

    .. versionchanged:: 0.11.0
       Moved into lib.mdamath
    """
    # TODO: enable typing when https://github.com/python/mypy/issues/11347 done
    normal: npt.NDArray = np.cross(vec1, vec2)  # type: ignore
    n = norm(normal)
    if n == 0.0:
        return normal  # returns [0,0,0] instead of [nan,nan,nan]
    # ... could also use numpy.nan_to_num(normal/norm(normal))
    return normal / n


def pdot(a: npt.NDArray, b: npt.NDArray) -> npt.NDArray:
    """Pairwise dot product.

    ``a`` must be the same shape as ``b``.

    Parameters
    ----------
    a: :class:`numpy.ndarray` of shape (N, M)
    b: :class:`numpy.ndarray` of shape (N, M)

    Returns
    -------
    :class:`numpy.ndarray` of shape (N,)
    """
    return np.einsum('ij,ij->i', a, b)


def pnorm(a: npt.NDArray) -> npt.NDArray:
    """Euclidean norm of each vector in a matrix

    Parameters
    ----------
    a: :class:`numpy.ndarray` of shape (N, M)

    Returns
    -------
    :class:`numpy.ndarray` of shape (N,)
    """
    return pdot(a, a)**0.5


def angle(a: npt.ArrayLike, b: npt.ArrayLike) -> float:
    """Returns the angle between two vectors in radians

    .. versionchanged:: 0.11.0
       Moved into lib.mdamath
    .. versionchanged:: 1.0.0
       Changed rounding-off code to use `np.clip()`. Values lower than
       -1.0 now return `np.pi` instead of `-np.pi`
    """
    x = np.dot(a, b) / (norm(a) * norm(b))
    # catch roundoffs that lead to nan otherwise
    x = np.clip(x, -1.0, 1.0)
    return np.arccos(x)


def stp(vec1: npt.ArrayLike, vec2: npt.ArrayLike, vec3: npt.ArrayLike) -> float:
    r"""Takes the scalar triple product of three vectors.

    Returns the volume *V* of the parallel epiped spanned by the three
    vectors

    .. math::

        V = \mathbf{v}_3 \cdot (\mathbf{v}_1 \times \mathbf{v}_2)

    .. versionchanged:: 0.11.0
       Moved into lib.mdamath
    """
    # TODO: enable typing when https://github.com/python/mypy/issues/11347 done
    return np.dot(vec3, np.cross(vec1, vec2))  # type: ignore


def dihedral(ab: npt.ArrayLike, bc: npt.ArrayLike, cd: npt.ArrayLike) -> float:
    r"""Returns the dihedral angle in radians between vectors connecting A,B,C,D.

    The dihedral measures the rotation around bc::

         ab
       A---->B
              \ bc
              _\'
                C---->D
                  cd

    The dihedral angle is restricted to the range -π <= x <= π.

    .. versionadded:: 0.8
    .. versionchanged:: 0.11.0
       Moved into lib.mdamath
    """
    x = angle(normal(ab, bc), normal(bc, cd))
    return (x if stp(ab, bc, cd) <= 0.0 else -x)


def sarrus_det(matrix: npt.NDArray) -> Union[float, npt.NDArray]:
    """Computes the determinant of a 3x3 matrix according to the
    `rule of Sarrus`_.

    If an array of 3x3 matrices is supplied, determinants are computed per
    matrix and returned as an appropriately shaped numpy array.

    .. _rule of Sarrus:
       https://en.wikipedia.org/wiki/Rule_of_Sarrus

    Parameters
    ----------
    matrix : numpy.ndarray
        An array of shape ``(..., 3, 3)`` with the 3x3 matrices residing in the
        last two dimensions.

    Returns
    -------
    det : float or numpy.ndarray
        The determinant(s) of `matrix`.
        If ``matrix.shape == (3, 3)``, the determinant will be returned as a
        scalar. If ``matrix.shape == (..., 3, 3)``, the determinants will be
        returned as a :class:`numpy.ndarray` of shape ``(...,)`` and dtype
        ``numpy.float64``.

    Raises
    ------
    ValueError:
        If `matrix` has less than two dimensions or its last two dimensions
        are not of shape ``(3, 3)``.


    .. versionadded:: 0.20.0
    """
    m = matrix.astype(np.float64)
    shape = m.shape
    ndim = m.ndim
    if ndim < 2 or shape[-2:] != (3, 3):
        raise ValueError("Invalid matrix shape: must be (3, 3) or (..., 3, 3), "
                         "got {}.".format(shape))
    if ndim == 2:
        return _sarrus_det_single(m)
    return _sarrus_det_multiple(m.reshape((-1, 3, 3))).reshape(shape[:-2])


def triclinic_box(x: npt.ArrayLike, y: npt.ArrayLike, z: npt.ArrayLike) -> npt.NDArray:
    """Convert the three triclinic box vectors to
    ``[lx, ly, lz, alpha, beta, gamma]``.

    If the resulting box is invalid, i.e., any box length is zero or negative,
    or any angle is outside the open interval ``(0, 180)``, a zero vector will
    be returned.

    All angles are in degrees and defined as follows:

    * ``alpha = angle(y,z)``
    * ``beta  = angle(x,z)``
    * ``gamma = angle(x,y)``

    Parameters
    ----------
    x : array_like
        Array of shape ``(3,)`` representing the first box vector
    y : array_like
        Array of shape ``(3,)`` representing the second box vector
    z : array_like
        Array of shape ``(3,)`` representing the third box vector

    Returns
    -------
    numpy.ndarray
        A numpy array of shape ``(6,)`` and dtype ``np.float32`` providing the
        unitcell dimensions in the same format as returned by
        :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:\n
        ``[lx, ly, lz, alpha, beta, gamma]``.\n
        Invalid boxes are returned as a zero vector.

    Note
    ----
    Definition of angles: http://en.wikipedia.org/wiki/Lattice_constant

    See Also
    --------
    :func:`~MDAnalysis.lib.mdamath.triclinic_vectors`


    .. versionchanged:: 0.20.0
       Calculations are performed in double precision and invalid box vectors
       result in an all-zero box.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    lx = norm(x)
    ly = norm(y)
    lz = norm(z)
    with np.errstate(invalid="ignore"):
        alpha = np.rad2deg(np.arccos(np.dot(y, z) / (ly * lz)))
        beta = np.rad2deg(np.arccos(np.dot(x, z) / (lx * lz)))
        gamma = np.rad2deg(np.arccos(np.dot(x, y) / (lx * ly)))
    box = np.array([lx, ly, lz, alpha, beta, gamma], dtype=np.float32)
    # Only positive edge lengths and angles in (0, 180) are allowed:
    if np.all(box > 0.0) and alpha < 180.0 and beta < 180.0 and gamma < 180.0:
        return box
    # invalid box, return zero vector:
    return np.zeros(6, dtype=np.float32)


def triclinic_vectors(dimensions: npt.ArrayLike,
                      dtype: npt.DTypeLike = np.float32) -> npt.NDArray:
    """Convert ``[lx, ly, lz, alpha, beta, gamma]`` to a triclinic matrix
    representation.

    Original `code by Tsjerk Wassenaar`_ posted on the Gromacs mailinglist.

    If `dimensions` indicates a non-periodic system (i.e., all lengths are
    zero), zero vectors are returned. The same applies for invalid `dimensions`,
    i.e., any box length is zero or negative, or any angle is outside the open
    interval ``(0, 180)``.

    .. _code by Tsjerk Wassenaar:
       http://www.mail-archive.com/gmx-users@gromacs.org/msg28032.html

    Parameters
    ----------
    dimensions : array_like
        Unitcell dimensions provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:\n
        ``[lx, ly, lz, alpha, beta, gamma]``.
    dtype: numpy.dtype
        The data type of the returned box matrix.

    Returns
    -------
    box_matrix : numpy.ndarray
        A numpy array of shape ``(3, 3)`` and dtype `dtype`,
        with ``box_matrix[0]`` containing the first, ``box_matrix[1]`` the
        second, and ``box_matrix[2]`` the third box vector.

    Notes
    -----
    * The first vector is guaranteed to point along the x-axis, i.e., it has the
      form ``(lx, 0, 0)``.
    * The second vector is guaranteed to lie in the x/y-plane, i.e., its
      z-component is guaranteed to be zero.
    * If any box length is negative or zero, or if any box angle is zero, the
      box is treated as invalid and an all-zero-matrix is returned.


    .. versionchanged:: 0.7.6
       Null-vectors are returned for non-periodic (or missing) unit cell.
    .. versionchanged:: 0.20.0
       * Calculations are performed in double precision and zero vectors are
         also returned for invalid boxes.
       * Added optional output dtype parameter.
    """
    dim = np.asarray(dimensions, dtype=np.float64)
    lx, ly, lz, alpha, beta, gamma = dim
    # Only positive edge lengths and angles in (0, 180) are allowed:
    if not (np.all(dim > 0.0) and
            alpha < 180.0 and beta < 180.0 and gamma < 180.0):
        # invalid box, return zero vectors:
        box_matrix = np.zeros((3, 3), dtype=dtype)
    # detect orthogonal boxes:
    elif alpha == beta == gamma == 90.0:
        # box is orthogonal, return a diagonal matrix:
        box_matrix = np.diag(dim[:3].astype(dtype, copy=False))
    # we have a triclinic box:
    else:
        box_matrix = np.zeros((3, 3), dtype=np.float64)
        box_matrix[0, 0] = lx
        # Use exact trigonometric values for right angles:
        if alpha == 90.0:
            cos_alpha = 0.0
        else:
            cos_alpha = np.cos(np.deg2rad(alpha))
        if beta == 90.0:
            cos_beta = 0.0
        else:
            cos_beta = np.cos(np.deg2rad(beta))
        if gamma == 90.0:
            cos_gamma = 0.0
            sin_gamma = 1.0
        else:
            gamma = np.deg2rad(gamma)
            cos_gamma = np.cos(gamma)
            sin_gamma = np.sin(gamma)
        box_matrix[1, 0] = ly * cos_gamma
        box_matrix[1, 1] = ly * sin_gamma
        box_matrix[2, 0] = lz * cos_beta
        box_matrix[2, 1] = lz * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        box_matrix[2, 2] = np.sqrt(lz * lz - box_matrix[2, 0] ** 2 -
                                   box_matrix[2, 1] ** 2)
        # The discriminant of the above square root is only negative or zero for
        # triplets of box angles that lead to an invalid box (i.e., the sum of
        # any two angles is less than or equal to the third).
        # We don't need to explicitly test for np.nan here since checking for a
        # positive value already covers that.
        if box_matrix[2, 2] > 0.0:
            # all good, convert to correct dtype:
            box_matrix = box_matrix.astype(dtype, copy=False)
        else:
            # invalid box, return zero vectors:
            box_matrix = np.zeros((3, 3), dtype=dtype)
    return box_matrix


def box_volume(dimensions: npt.ArrayLike) -> float:
    """Return the volume of the unitcell described by `dimensions`.

    The volume is computed as the product of the box matrix trace, with the
    matrix obtained from :func:`triclinic_vectors`.
    If the box is invalid, i.e., any box length is zero or negative, or any
    angle is outside the open interval ``(0, 180)``, the resulting volume will
    be zero.

    Parameters
    ----------
    dimensions : array_like
        Unitcell dimensions provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:\n
        ``[lx, ly, lz, alpha, beta, gamma]``.

    Returns
    -------
    volume : float
        The volume of the unitcell. Will be zero for invalid boxes.


    .. versionchanged:: 0.20.0
        Calculations are performed in double precision and zero is returned
        for invalid dimensions.
    """
    dim = np.asarray(dimensions, dtype=np.float64)
    lx, ly, lz, alpha, beta, gamma = dim
    if alpha == beta == gamma == 90.0 and lx > 0 and ly > 0 and lz > 0:
        # valid orthogonal box, volume is the product of edge lengths:
        volume = lx * ly * lz
    else:
        # triclinic or invalid box, volume is trace product of box matrix
        # (invalid boxes are set to all zeros by triclinic_vectors):
        tri_vecs = triclinic_vectors(dim, dtype=np.float64)
        volume = tri_vecs[0, 0] * tri_vecs[1, 1] * tri_vecs[2, 2]
    return volume
