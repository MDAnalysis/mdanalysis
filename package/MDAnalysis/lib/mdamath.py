# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Mathematical helper functions --- :mod:`MDAnalysis.lib.mdamath`
===============================================================

Helper functions for common mathematical operations

.. autofunction:: triclinic_box
.. autofunction:: triclinic_vectors
.. autofunction:: box_volume

.. versionadded:: 0.11.0
"""
import numpy as np


def _veclength(v):
    """Length of vector *v*."""
    # note: this is 3 times faster than numpy.linalg.norm
    return np.sqrt(np.dot(v, v))


def _angle(a, b):
    """Angle between two vectors *a* and *b* in degrees.

    If one of the lengths is 0 then the angle is returned as 0
    (instead of `nan`).
    """
    angle = np.arccos(np.dot(a, b) / (_veclength(a) * _veclength(b)))
    if np.isnan(angle):
        return 0.0
    return np.rad2deg(angle)


def triclinic_box(x, y, z):
    """Convert the three triclinic box vectors to [A,B,C,alpha,beta,gamma].

    Angles are in degrees.

    * alpha  = angle(y,z)
    * beta   = angle(x,z)
    * gamma  = angle(x,y)

    .. SeeAlso:: Definition of angles: http://en.wikipedia.org/wiki/Lattice_constant
    """
    A, B, C = [_veclength(v) for v in x, y, z]
    alpha = _angle(y, z)
    beta = _angle(x, z)
    gamma = _angle(x, y)
    return np.array([A, B, C, alpha, beta, gamma], dtype=np.float32)


def triclinic_vectors(dimensions):
    """Convert `[A,B,C,alpha,beta,gamma]` to a triclinic box representation.

    Original `code by Tsjerk Wassenaar`_ posted on the Gromacs mailinglist.

    If *dimensions* indicates a non-periodic system (i.e. all lengths
    0) then null-vectors are returned.

    .. _code by Tsjerk Wassenaar:
       http://www.mail-archive.com/gmx-users@gromacs.org/msg28032.html

    :Arguments:
      *dimensions*
        list of box lengths and angles (in degrees) such as
        [A,B,C,alpha,beta,gamma]

    :Returns: numpy 3x3 array B, with B[0] = first box vector,
              B[1] = second vector, B[2] third box vector.

    .. note::

       The first vector is always pointing along the X-axis
       i.e. parallel to (1,0,0).

    .. versionchanged:: 0.7.6
       Null-vectors are returned for non-periodic (or missing) unit cell.

    """
    B = np.zeros((3, 3), dtype=np.float32)
    x, y, z, a, b, c = dimensions[:6]

    if np.all(dimensions[:3] == 0):
        return B

    B[0][0] = x
    if a == 90. and b == 90. and c == 90.:
        B[1][1] = y
        B[2][2] = z
    else:
        a = np.deg2rad(a)
        b = np.deg2rad(b)
        c = np.deg2rad(c)
        B[1][0] = y * np.cos(c)
        B[1][1] = y * np.sin(c)
        B[2][0] = z * np.cos(b)
        B[2][1] = z * (np.cos(a) - np.cos(b) * np.cos(c)) / np.sin(c)
        B[2][2] = np.sqrt(z * z - B[2][0] ** 2 - B[2][1] ** 2)
    return B


def box_volume(dimensions):
    """Return the volume of the unitcell described by *dimensions*.

    The volume is computed as `det(x1,x2,x2)` where the xi are the
    triclinic box vectors from :func:`triclinic_vectors`.

    :Arguments:
       *dimensions*
          list of box lengths and angles (in degrees) such as
          [A,B,C,alpha,beta,gamma]

    :Returns: numpy 3x3 array B, with B[0] = first box vector,
              B[1] = second vector, B[2] third box vector.
    """
    return np.linalg.det(triclinic_vectors(dimensions))
