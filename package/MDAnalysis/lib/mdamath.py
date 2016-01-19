# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
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

.. autofunction:: normal
.. autofunction:: norm
.. autofunction:: angle
.. autofunction:: dihedral
.. autofunction:: stp
.. autofunction:: triclinic_box
.. autofunction:: triclinic_vectors
.. autofunction:: box_volume
.. autofunction:: make_whole

.. versionadded:: 0.11.0
"""
import numpy as np

from ..exceptions import NoDataError


# geometric functions
def norm(v):
    r"""Calculate the norm of a vector v.

    .. math:: v = \sqrt{\mathbf{v}\cdot\mathbf{v}}

    This version is faster then numpy.linalg.norm because it only works for a
    single vector and therefore can skip a lot of the additional fuss
    linalg.norm does.

    Parameters
    ----------
    v: array_like
        1D array of shape (N) for a vector of length N

    Returns
    -------
    float
        norm of the vector

    """
    return np.sqrt(np.dot(v, v))


def normal(vec1, vec2):
    r"""Returns the unit vector normal to two vectors.

    .. math::

       \hat{\mathbf{n}} = \frac{\mathbf{v}_1 \times \mathbf{v}_2}{|\mathbf{v}_1 \times \mathbf{v}_2|}

    If the two vectors are collinear, the vector :math:`\mathbf{0}` is returned.

    .. versionchanged:: 0.11.0
       Moved into lib.mdamath
    """
    normal = np.cross(vec1, vec2)
    n = norm(normal)
    if n == 0.0:
        return normal  # returns [0,0,0] instead of [nan,nan,nan]
    return normal / n  # ... could also use numpy.nan_to_num(normal/norm(normal))


def angle(a, b):
    """Returns the angle between two vectors in radians

    .. versionchanged:: 0.11.0
       Moved into lib.mdamath
    """
    x = np.dot(a, b) / (norm(a) * norm(b))
    # catch roundoffs that lead to nan otherwise
    if x > 1.0:
        return 0.0
    elif x < -1.0:
        return -np.pi
    return np.arccos(x)


def stp(vec1, vec2, vec3):
    r"""Takes the scalar triple product of three vectors.

    Returns the volume *V* of the parallel epiped spanned by the three
    vectors

    .. math::

        V = \mathbf{v}_3 \cdot (\mathbf{v}_1 \times \mathbf{v}_2)

    .. versionchanged:: 0.11.0
       Moved into lib.mdamath
    """
    return np.dot(vec3, np.cross(vec1, vec2))


def dihedral(ab, bc, cd):
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


def _angle(a, b):
    """Angle between two vectors *a* and *b* in degrees.

    If one of the lengths is 0 then the angle is returned as 0
    (instead of `nan`).
    """
    # This function has different limits than angle?

    angle = np.arccos(np.dot(a, b) / (norm(a) * norm(b)))
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
    A, B, C = [norm(v) for v in (x, y, z)]
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


def _is_contiguous(atomgroup, atom):
    """Walk through atomgroup, starting with atom.

    :Returns:
       ``True`` if all of *atomgroup* is accessible through walking
        along bonds.
       ``False`` otherwise.
    """
    seen = set([atom])
    walked = set()
    ag_set = set(atomgroup)

    nloops = 0
    while len(seen) < len(atomgroup):
        nloops += 1
        if nloops > len(atomgroup):
            return False

        todo = seen.difference(walked)
        for atom in todo:
            for other in atom.bonded_atoms:
                if other in ag_set:
                    seen.add(other)
            walked.add(atom)

    return True


def make_whole(atomgroup, reference_atom=None):
    """Move all atoms in a single molecule so that bonds don't split over images

    :Arguments:
      *atomgroup*
        The :class:`MDAnalysis.core.AtomGroup.AtomGroup` to work with.
        The positions of this are modified in place.  All these atoms
        must belong in the same molecule or fragment.

    :Keywords:
      *reference_atom*
        The atom around which all other atoms will be moved.
        Defaults to atom 0 in the atomgroup.

    :Returns:
       ``None``.  Atom positions are modified in place

    :Raises:
      `NoDataError`
         if there are no bonds present.
         (See :func:`~MDAnalysis.topology.core.guess_bonds`)

      `ValueError`
         if the algorithm fails to work.  This is usually
         caused by the atomgroup not being a single fragment.
         (ie the molecule can't be traversed by following bonds)

    This function is most useful when atoms have been packed into the primary
    unit cell, causing breaks mid molecule, with the molecule then appearing
    on either side of the unit cell. This is problematic for operations
    such as calculating the center of mass of the molecule. ::

       +-----------+     +-----------+
       |           |     |           |
       | 6       3 |     |         3 | 6
       | !       ! |     |         ! | !
       |-5-8   1-2-| ->  |       1-2-|-5-8
       | !       ! |     |         ! | !
       | 7       4 |     |         4 | 7
       |           |     |           |
       +-----------+     +-----------+

    Usage::

      from MDAnalysis.util.mdamath import make_whole

      # This algorithm requires bonds, these can be guessed!
      u = mda.Universe(......, guess_bonds=True)

      # MDAnalysis can split molecules into their fragments
      # based on bonding information.
      # Note that this function will only handle a single fragment
      # at a time, necessitating a loop.
      for frag in u.fragments:
          make_whole(frag)

    Alternatively, to keep a single atom in place as the anchor::

      # This will mean that atomgroup[10] will NOT get moved,
      # and all other atoms will move (if necessary).
      make_whole(atomgroup, reference_atom=atomgroup[10])

    .. Note::

       Only orthogonal boxes are currently supported

    .. versionadded:: 0.11.0
    """
    if not atomgroup.bonds:
        raise NoDataError("The atomgroup is required to have bonds")

    if reference_atom is None:
        ref = atomgroup[0]
    else:
        ref = reference_atom
        # Sanity check
        if not ref in atomgroup:
            raise ValueError("Reference atom not in atomgroup")

    # Check all of atomgroup is accessible from ref
    if not _is_contiguous(atomgroup, ref):
        raise ValueError("atomgroup not contiguous from bonds")

    # Not sure if this is actually a requirement...
    # I think application of pbc would need to be changed for triclinic boxes
    # but that's all?  How does minimum bond length criteria change?
    if not all(atomgroup.dimensions[3:] == 90.0):
        raise ValueError("Non orthogonal boxes are not supported")
    box = atomgroup.dimensions[:3]

    if all(box == 0.0):
        raise ValueError("Supplied box had zero size")

    box_length = box.min() / 2.0

    bondlengths = atomgroup.bonds.bonds(pbc=True)
    if bondlengths.min() * 1.4 > box_length:
        raise ValueError("Box lengths are too small relative to bond lengths")

    # All checks done, let's continue
    # If bond lengths don't change after pbc applied, then no bonds
    # straddle the box boundaries
    if np.allclose(atomgroup.bonds.bonds(), bondlengths):
        return
    # Can't reuse this calculation of bond lengths as we're changing
    # stuff as we go.

    processed = set()  # Who have I already done?
    ref_points = set([ref])  # Who is safe to use as reference point?

    ag_set = set(atomgroup)
    nres = len(atomgroup)  # total size of the problem
    nloops = 0
    while len(ref_points) < nres:  # While all atoms aren't correct
        nloops += 1
        if nloops > nres:  # To prevent infinite loop
            # This point probable isn't reachable with the above _is_contiguous
            # check, but better safe than sorry.
            raise ValueError("Algorithm couldn't traverse atomgroup. "
                             "Perhaps the atomgroup isn't fully connected")

        # We want to iterate over atoms that are good to use as reference
        # points, but haven't been processed yet.
        todo = ref_points.difference(processed)
        for atom in todo:
            for b in atom.bonds:
                other = b.partner(atom)
                # Avoid atoms not in our scope
                if other not in ag_set:
                    continue
                if other in ref_points:
                    continue
                if b.length() > box_length:
                    # Vector from ref atom to other
                    vec = other.position - ref.position
                    # Apply pbc to this vector
                    vec -= np.rint(vec/box) * box
                    # Define the position of other based on this vector
                    other.position = ref.position + vec
                # This atom can now be used as a reference point
                ref_points.add(other)

            processed.add(atom)
