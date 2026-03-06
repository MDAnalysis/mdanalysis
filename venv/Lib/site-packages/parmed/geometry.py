"""
This module contains the functionality for carrying out geometrical calculations
for molecules and molecular systems

Author: Jason Swails
Contributors:

Copyright (C) 2014 - 2015  Jason Swails

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
59 Temple Place - Suite 330
Boston, MA 02111-1307, USA.
"""
import warnings
from collections import defaultdict
from math import pi, cos, sin, sqrt, acos

import numpy as np

from . import unit as u
from .constants import TINY, DEG_TO_RAD, RAD_TO_DEG
from .vec3 import Vec3

def box_lengths_and_angles_to_vectors(a, b, c, alpha, beta, gamma):
    """
    This function takes the lengths of the unit cell vectors and the angles
    between them and returns 3 unit cell vectors satisfying those dimensions

    Parameters
    ----------
    a : double (or length Quantity)
        Length of the first unit cell vector
    b : double (or length Quantity)
        Length of the second unit cell vector
    c : double (or length Quantity)
        Length of the third unit cell vector
    alpha : double (or angle Quantity)
        Angle between vectors b and c
    beta : double (or angle Quantity)
        Angle between vectors a and c
    gamma : double (or angle Quantity)
        Angle between vectors a and b

    Returns
    -------
    list Quantity, list Quantity, list Quantity
        The 3, 3-element vectors as quantities with dimension length

    Notes
    -----
    The unit cell lengths are assumed to be Angstroms if no explicit unit is
    given. The angles are assumed to be degrees
    """
    if u.is_quantity(a): a = a.value_in_unit(u.angstroms)
    if u.is_quantity(b): b = b.value_in_unit(u.angstroms)
    if u.is_quantity(c): c = c.value_in_unit(u.angstroms)
    if u.is_quantity(alpha): alpha = alpha.value_in_unit(u.degrees)
    if u.is_quantity(beta): beta = beta.value_in_unit(u.degrees)
    if u.is_quantity(gamma): gamma = gamma.value_in_unit(u.degrees)

    if alpha <= 2*pi and beta <= 2*pi and gamma <= 2*pi:
        warnings.warn('Strange unit cell vector angles detected. They '
                      'appear to be in radians...')

    alpha *= DEG_TO_RAD
    beta *= DEG_TO_RAD
    gamma *= DEG_TO_RAD

    av = [a, 0.0, 0.0]
    bx = b * cos(gamma)
    by = b * sin(gamma)
    bv = [bx, by, 0.0]
    cx = c * cos(beta)
    cy = c * (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma)
    cz = sqrt(c*c - cx*cx - cy*cy)
    cv = [cx, cy, cz]

    # Make sure that any tiny components are exactly 0
    if abs(bx) < TINY: bv[0] = 0
    if abs(by) < TINY: bv[1] = 0
    if abs(cx) < TINY: cv[0] = 0
    if abs(cy) < TINY: cv[1] = 0
    if abs(cz) < TINY: cv[2] = 0

    return (av, bv, cv) * u.angstroms

def box_vectors_to_lengths_and_angles(a, b, c):
    """
    This function takes the lengths of the unit cell vectors and the angles
    between them and returns 3 unit cell vectors satisfying those dimensions

    Parameters
    ----------
    a : collection of 3 floats (or length Quantity)
        The first unit cell vector
    b : collection of 3 floats (or length Quantity)
        The second unit cell vector
    c : collection of 3 floats (or length Quantity)
        The third unit cell vector

    Returns
    -------
    (a, b, c), (alpha, beta, gamma)
        Two tuples, the first is the 3 unit cell vector lengths as
        length-dimension Quantity objects and the second is the set of angles
        between the unit cell vectors as angle-dimension Quantity objects

    Notes
    -----
    The unit cell lengths are assumed to be Angstroms if no explicit unit is
    given.
    """
    if u.is_quantity(a): a = a.value_in_unit(u.angstroms)
    if u.is_quantity(b): b = b.value_in_unit(u.angstroms)
    if u.is_quantity(c): c = c.value_in_unit(u.angstroms)
    # Get the lengths
    la = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
    lb = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2])
    lc = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2])
    # Angles
    alpha = acos((b[0]*c[0] + b[1]*c[1] + b[2]*c[2]) / (lb*lc))
    beta = acos((a[0]*c[0] + a[1]*c[1] + a[2]*c[2]) / (la*lc))
    gamma = acos((b[0]*a[0] + b[1]*a[1] + b[2]*a[2]) / (lb*la))
    # Convert to degrees
    alpha *= RAD_TO_DEG
    beta *= RAD_TO_DEG
    gamma *= RAD_TO_DEG

    return (la, lb, lc) * u.angstroms, (alpha, beta, gamma) * u.degrees

def reduce_box_vectors(a, b, c):
    """
    This function puts three unit cell vectors in a reduced form where a is
    "mostly" in x, b is "mostly" in y, and c is "mostly" in z. This form is
    necessary for some programs (notably OpenMM and Gromacs)

    Parameters
    ----------
    a : 3-element collection of float
        First unit cell vector
    b : 3-element collection of float
        Second unit cell vector
    c : 3-element collection of float
        Third unit cell vector

    Returns
    -------
    red_a, red_b, red_c : Vec3, Vec3, Vec3
        The reduced unit cell vectors in units of angstroms

    Notes
    -----
    The implementation here is taken from the OpenMM Python application layer
    written by Peter Eastman
    """
    if u.is_quantity(a):
        a = a.value_in_unit(u.angstroms)
    if u.is_quantity(b):
        b = b.value_in_unit(u.angstroms)
    if u.is_quantity(c):
        c = c.value_in_unit(u.angstroms)

    a = Vec3(*a)
    b = Vec3(*b)
    c = Vec3(*c)

    c = c - b*round(c[1]/b[1])
    c = c - a*round(c[0]/a[0])
    b = b - a*round(b[0]/a[0])

    return a, b, c


def center_of_mass(coordinates, masses):
    """ Compute the center of mass of a group of coordinates.

    Parameters
    ----------
    coordinates : numpy.ndarray
        Coordinate array
    masses : numpy.ndarray
        Array of masses

    Returns
    -------
    COM
        np.ndarray of shape (3,) identifying the cartesian center of mass

    Notes
    -----
    This method *requires* that the parameters be passed in as numpy arrays.
    AttributeError's will ensue if this is not the case. Also, coordinates must
    be able to be reshaped to (len(masses), 3), or ValueError's will ensue
    """
    masses = masses.flatten()
    coordinates = coordinates.reshape((masses.shape[0], 3))
    return np.average(coordinates, weights=masses, axis=0)

def distance2(a1, a2):
    """ Computes the cartesian distance between two atoms. Ignores periodic boundary conditions.

    Parameters
    ----------
    a1, a2 : Atom or collection of 3 coordinates
        The two atoms between whom the distance should be calculated

    Returns
    -------
    d2 : float
        The square of the distance between the two atoms

    Notes
    -----
    This is done in pure Python, so it should not be used for large numbers of
    distance calculations. For that, use numpy-vectorized routines and the numpy
    coordinate arrays

    Raises
    ------
    TypeError if a1 or a2 are not Atom or iterable
    ValueError if a1 or a2 are iterable, but do not have exactly 3 items
    """
    x1, y1, z1 = _get_coords_from_atom_or_tuple(a1)
    x2, y2, z2 = _get_coords_from_atom_or_tuple(a2)
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2
    return dx*dx + dy*dy + dz*dz

def angle(a1, a2, a3):
    """ Computes the cartesian angle between three atoms. Ignores periodic boundary conditions.

    Parameters
    ----------
    a1, a2, a3 : Atom or collection of 3 coordinates
        The three atoms between whom the angle should be calculated (with a2
        being the central atoms)

    Returns
    -------
    ang : float
        The angle between the vectors a1-a2 and a2-a3 in degrees

    Notes
    -----
    This is done in pure Python, so it should not be used for large numbers of
    distance calculations. For that, use numpy-vectorized routines and the numpy
    coordinate arrays

    Raises
    ------
    TypeError if a1, a2, or a3 are not Atom or iterable
    ValueError if a1, a2, or a3 are iterable, but do not have exactly 3 items
    """
    x1, y1, z1 = _get_coords_from_atom_or_tuple(a1)
    x2, y2, z2 = _get_coords_from_atom_or_tuple(a2)
    x3, y3, z3 = _get_coords_from_atom_or_tuple(a3)
    v1 = np.array([x2 - x1, y2 - y1, z2 - z1])
    v2 = np.array([x2 - x3, y2 - y3, z2 - z3])
    l1 = np.sqrt(np.dot(v1, v1))
    l2 = np.sqrt(np.dot(v2, v2))
    cosa = np.dot(v1, v2) / (l1 * l2)
    return np.degrees(np.arccos(cosa))

def dihedral(a1, a2, a3, a4):
    """
    Computes the angle between three vectors made up of four points (all three
    vectors share one point with one other vector)

    Parameters
    ----------
    a1, a2, a3, a4 : Atom or collection of 4 coordinates
        The four atoms between whom the torsion angle should be calculated (with
        a1 and a4 being the two end-point atoms not shared between two vectors)

    Returns
    -------
    dihed : float
        The measured dihedral between the 4 points in degrees
    """
    p = np.array([_get_coords_from_atom_or_tuple(a1),
                  _get_coords_from_atom_or_tuple(a2),
                  _get_coords_from_atom_or_tuple(a3),
                  _get_coords_from_atom_or_tuple(a4)])
    v1 = p[1] - p[0]
    v2 = p[1] - p[2]
    v3 = p[3] - p[2]
    # Take the cross product between v1-v2 and v2-v3
    v1xv2 = _cross(v1, v2)
    v2xv3 = _cross(v2, v3)
    # Now find the angle between these cross-products
    l1 = np.sqrt(np.dot(v1xv2, v1xv2))
    l2 = np.sqrt(np.dot(v2xv3, v2xv3))
    cosa = np.dot(v1xv2, v2xv3) / (l1 * l2)
    if np.dot(v3, v1xv2) <= 0.0 :
        return np.degrees(np.arccos(cosa))
    else :
        return -np.degrees(np.arccos(cosa))

def _cross(v1, v2):
    """ Computes the cross-product """
    # Can't use np.cross for pypy, since it's not yet implemented
    return np.array([v1[1]*v2[2] - v1[2]*v2[1],
                     v1[2]*v2[0] - v1[0]*v2[2],
                     v1[0]*v2[1] - v1[1]*v2[0]])

def _unit_diff(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Unit vector going from a to b"""
    u = (b - a).astype(np.float64)
    u /= sqrt(np.dot(u, u))
    return u

def _unit_cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Unit vector along cross product of a and b"""
    u = _cross(a, b)
    u /= sqrt(np.dot(u, u))
    return u

# Skipping type annotations for Atom bc it causes a circular import
def local_axes(atom1, atom2, atom3) -> np.ndarray:
    """Return an orthonormal basis defined by three atoms.

    This can be used to bootstrap a cartesian coordinate system from internal
    coordinates by choosing the origin and orientation based on a choice of
    three atoms. In the present convention the first atom is placed at the
    origin, the second atom is placed on the positive side of the z-axis and
    the third atom is placed on the positive side of the x-axis.
    """
    x1 = np.array([atom1.xx, atom1.xy, atom1.xz])
    x2 = np.array([atom2.xx, atom2.xy, atom2.xz])
    x3 = np.array([atom3.xx, atom3.xy, atom3.xz])
    z = _unit_diff(x1, x2)
    u23 = _unit_diff(x2, x3)
    if np.isclose(abs(np.dot(z, u23)), 1.0):
        raise ValueError('Atoms defining local_frame are co-linear')
    x = _unit_cross(z, _unit_cross(u23, z))
    y = _unit_cross(z, x)
    return np.array([x, y, z])

# Skipping type annotations for Atom bc it causes a circular import
def bat2xyz(atom1, atom2, atom3, atom4) -> np.ndarray:
    """
    Compute cartesian coordinates for atom4 using its bond/angle/torsion
    coordinates relative to atoms 3, 2 and 1.
    """
    r = sqrt(distance2(atom3, atom4))
    th = np.radians(angle(atom2, atom3, atom4))
    ph = np.radians(dihedral(atom1, atom2, atom3, atom4))
    sin_th = sin(th)
    x_over_r = sin_th * cos(ph) # This matches the definition in local_axes()
    y_over_r = sin_th * sin(ph)
    z_over_r = cos(th)
    return r * np.array([x_over_r, y_over_r, z_over_r])

def _get_coords_from_atom_or_tuple(a):
    from .topologyobjects import Atom
    if isinstance(a, Atom):
        return a.xx, a.xy, a.xz
    return a

# tuples are pairs of atomic numbers followed by the distance cutoff below which
# they are considered "bonded". This is taken from Atom.cpp in cpptraj
# (Atom::GetBondLength)
_OFFSET = 0.20 # offset for what is considered a bond
_DEFAULT_CUTOFF = (1.6 + _OFFSET)**2
STANDARD_BOND_LENGTHS_SQUARED = defaultdict(lambda: _DEFAULT_CUTOFF)
STANDARD_BOND_LENGTHS_SQUARED.update({
    # Self-bonds
    (1, 1) : (0.74 + _OFFSET)**2,
    (6, 6) : (1.54 + _OFFSET)**2,
    (7, 7) : (1.45 + _OFFSET)**2,
    (8, 8) : (1.48 + _OFFSET)**2,
    (15, 15) : (2.21 + _OFFSET)**2,
    (16, 16) : (2.05 + _OFFSET)**2,
    # H-
    (1, 6) : (1.09 + _OFFSET)**2, (6, 1) : (1.09 + _OFFSET)**2,
    (1, 7) : (1.01 + _OFFSET)**2, (7, 1) : (1.01 + _OFFSET)**2,
    (1, 8) : (0.96 + _OFFSET)**2, (8, 1) : (0.96 + _OFFSET)**2,
    (1, 15) : (1.44 + _OFFSET)**2, (15, 1) : (1.44 + _OFFSET)**2,
    (1, 16) : (1.34 + _OFFSET)**2, (16, 1) : (1.34 + _OFFSET)**2,
    # C-
    (6, 7) : (1.47 + _OFFSET)**2, (7, 6) : (1.47 + _OFFSET)**2,
    (6, 8) : (1.43 + _OFFSET)**2, (8, 6) : (1.43 + _OFFSET)**2,
    (6, 9) : (1.35 + _OFFSET)**2, (9, 6) : (1.35 + _OFFSET)**2,
    (6, 15) : (1.84 + _OFFSET)**2, (15, 6) : (1.84 + _OFFSET)**2,
    (6, 16) : (1.82 + _OFFSET)**2, (16, 6) : (1.82 + _OFFSET)**2,
    (6, 17) : (1.77 + _OFFSET)**2, (17, 6) : (1.77 + _OFFSET)**2,
    (6, 35) : (1.94 + _OFFSET)**2, (35, 6) : (1.94 + _OFFSET)**2,
    # N-
    (7, 8) : (1.40 + _OFFSET)**2, (8, 7) : (1.40 + _OFFSET)**2,
    (7, 9) : (1.36 + _OFFSET)**2, (9, 7) : (1.36 + _OFFSET)**2,
    (7, 15) : (1.71 + _OFFSET)**2, (15, 7) : (1.71 + _OFFSET)**2,
    (7, 16) : (1.68 + _OFFSET)**2, (16, 7) : (1.68 + _OFFSET)**2,
    (7, 17) : (1.75 + _OFFSET)**2, (17, 7) : (1.75 + _OFFSET)**2,
    # O-
    (8, 9) : (1.42 + _OFFSET)**2, (9, 8) : (1.42 + _OFFSET)**2,
    (8, 15) : (1.63 + _OFFSET)**2, (15, 8) : (1.63 + _OFFSET)**2,
    (8, 16) : (1.48 + _OFFSET)**2, (16, 8) : (1.48 + _OFFSET)**2,
    # F-
    (9, 15) : (1.54 + _OFFSET)**2, (15, 9) : (1.54 + _OFFSET)**2,
    (9, 16) : (1.56 + _OFFSET)**2, (16, 9) : (1.56 + _OFFSET)**2,
    # P-
    (15, 16) : (1.86 + _OFFSET)**2, (16, 15) : (1.86 + _OFFSET)**2,
    (15, 17) : (2.03 + _OFFSET)**2, (17, 15) : (2.03 + _OFFSET)**2,
    # S-
    (16, 17) : (2.07 + _OFFSET)**2, (17, 16) : (2.07 + _OFFSET)**2,
})
