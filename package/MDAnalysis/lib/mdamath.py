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
from __future__ import division, absolute_import
from six.moves import zip
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
    v : array_like
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

    Note
    ----
    Definition of angles: http://en.wikipedia.org/wiki/Lattice_constant
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

    Parameters
    ----------
    dimensions : [A, B, C, alpha, beta, gamma]
        list of box lengths and angles (in degrees) such as
        ``[A,B,C,alpha,beta,gamma]``

    Returns
    -------
    numpy.array
        numpy 3x3 array B, with ``B[0]`` as first box vector,
        ``B[1]``as second vector, ``B[2]`` as third box vector.

    Notes
    -----
    The first vector is always pointing along the X-axis,
    i.e., parallel to (1, 0, 0).


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

    Parameters
    ----------
    dimensions : [A, B, C, alpha, beta, gamma]
        list of box lengths and angles (in degrees) such as
        [A,B,C,alpha,beta,gamma]

    Returns
    -------
    volume : float
    """
    return np.linalg.det(triclinic_vectors(dimensions))


def _is_contiguous(atomgroup, atom):
    """Walk through atomgroup, starting with atom.

    Returns
    -------
    bool
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


def make_whole(atomgroup):
    """Move all atoms in a single molecule so that bonds don't split over images

    Atoms are modified in place.

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


    Parameters
    ----------
    atomgroup : AtomGroup
        The :class:`MDAnalysis.core.groups.AtomGroup` to work with.
        The positions of this are modified in place.  All these atoms
        must belong in the same molecule or fragment.

    Raises
    ------
    NoDataError
        There are no bonds present.
        (See :func:`~MDAnalysis.topology.core.guess_bonds`)

    ValueError
        The algorithm fails to work.  This is usually
        caused by the atomgroup not being a single fragment.
        (ie the molecule can't be traversed by following bonds)


    Example
    -------
    Make fragments whole::

        from MDAnalysis.transformations.pbc import make_whole

        # This algorithm requires bonds, these can be guessed!
        u = mda.Universe(......, guess_bonds=True)

        # MDAnalysis can split molecules into their fragments
        # based on bonding information.
        # Note that this function will only handle a single fragment
        # at a time, necessitating a loop.
        for frag in u.fragments:
          make_whole(frag)


    """
    from ..lib.distances import apply_PBC
    
    # check if box is valid
    box = atomgroup.dimensions[:3]
    if all(box == 0.0):
        raise ValueError("Supplied box had zero size")
    # this will fail if the atomgroup has no bonds
    try:
        b = atomgroup.bonds
    except (AttributeError, NoDataError):
        raise NoDataError("The atomgroup is required to have bonds")
    # check if box is too small
    box_length = box.min() / 2.0
    bondlengths = atomgroup.bonds.bonds(pbc=True)
    if bondlengths.min() * 1.4 > box_length:
        raise ValueError("Box lengths are too small relative to bond lengths")
    
    tric = triclinic_vectors(atomgroup.dimensions)
    halfbox = tric.sum(axis=0) / 2.0
    atgp = atomgroup
    diffs = None
    count = 0
    idx = None
    prev = {}
    repeats = 0
    while True:
        # this method is iterative, as correcting a bond causes other to become distorted
        if diffs is not None:
            if idx is not None:
                prev[count] = idx
            idx = np.unique(np.where(diffs)[0])
            atgp = atgp.bonds.atom2[idx]
        # get the bond vectors and translate them to the center of the unit cell
        # all the vectors that are bigger than half the unit cell will be outside
        # these are the ones we want to "correct"
        vecs = atgp.bonds.atom2.positions - atgp.bonds.atom1.positions
        vecs += halfbox
        # apply PBC on the vectors
        newvecs = apply_PBC(vecs, atgp.dimensions)
        # diffs is used to check which the atoms are already in their correct positions
        # if all of them are correct then break the loop
        diffs = np.any(~np.isclose(vecs, newvecs, 1e-6), axis=1)
        if not diffs.sum():
            break
        # the new vectors are translated back to the origin and positions are updated
        newvecs -= halfbox
        atgp.bonds.atom2.positions = atgp.bonds.atom1.positions + newvecs
        
        # the following check is for cases when an atom in the atom2 group
        # is connected to two atoms in atom1. Since where changing the positions
        # of the atom2 group, suchs cases lead to the while loop being stuck:
        # this atom is moved back and forth between the two atoms it is conected
        # to. Since were tracking the atoms that need changing, if this group stays
        # the same over multiple iterations, we change atom1.positions instead.
        if diffs is not None and idx is not None:
            if prev.has_key(count) and np.array_equal(prev[count], idx):
                repeats +=1 
                if repeats > len(idx):
                    atgp.bonds.atom1.positions = atgp.bonds.atom2.positions - newvecs
            else:
                repeats = 0
        count += 1
        # for debugging purposes. Also makes a nice animation
        #atomgroup.write(filename='updated_%d.pdb'%count)


def one_to_many_pointers(Ni, Nj, i2j):
    """Based on a many to one mapping of i to j, create the reverse mapping

    Arguments
    ---------
    Ni : int
       number of i  components
    Nj : int
       number of j components
    i2j : numpy.ndarray
       relates i to parent js

    Returns
    -------
    ordered : list
       an ordered list of i indices [size (i,)]
    ptrs : list
       the start and end index for each j [size (Nj, 2)]

    Example
    -------
    .. code:: python

        # Residx - the resid of each Atom
        ordered, ptrs = one_to_many_pointers(Natoms, Nres, Residx)

        # Returns an array of the atom indices that are in resid 7
        atoms = ordered[ptrs[7,0]:ptrs[7,1]]

    """
    ordered = i2j.argsort()
    sorted_idx = i2j[ordered]
    borders = np.concatenate([[0],
                              np.where(np.diff(sorted_idx))[0] + 1,
                              [Ni]])

    ptrs = np.zeros((Nj, 2), dtype=np.int32)
    for x, y in zip(borders[:-1], borders[1:]):
        i = sorted_idx[x]
        ptrs[i] = x, y

    return ordered, ptrs
