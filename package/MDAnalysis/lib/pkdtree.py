# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
PeriodicKDTree --- :mod:`MDAnalysis.lib.pkdtree`
===============================================================================

This module contains a class to allow searches on a KDTree involving periodic
boundary conditions.
"""

from __future__ import absolute_import
from six.moves import range

import numpy as np
from Bio.KDTree import _CKDTree

from MDAnalysis.lib.distances import _box_check, _check_array, apply_PBC
from MDAnalysis.lib.mdamath import norm, triclinic_vectors, triclinic_box

__all__ = ['PeriodicKDTree', ]


class PeriodicKDTree(object):
    """
    Wrapper around Bio.KDTree._CKDTree to enable search with periodic boundary
     conditions.

    A tree is first constructed with the coordinates wrapped onto the central
    cell. A query for neighbors around a center point is performed first by
    wrapping the center point coordinates to the central cell, then generating
    images of this wrapped center point (thanks to
    https://github.com/patvarilly/periodic_kdtree for the idea) and finally
    searching for neighbors close to the images.

    Only the necessary number of center point images is generated for each
    case. For instance, if the wrapped center point lies well within the cell
    and far away from the cell boundaries, there may be no need to generate
    any center point image.
    """

    def __init__(self, box, bucket_size=10):
        """
        Parameters
        ----------
        box : array-like or ``None``, optional, default ``None``
          Simulation cell dimensions in the form of
          :attr:`MDAnalysis.trajectory.base.Timestep.dimensions` when
          periodic boundary conditions should be taken into account for
          the calculation of contacts.
        bucket_size : int
          Number of entries in leafs of the KDTree. If you suffer poor
          performance you can play around with this number. Increasing the
          `bucket_size` will speed up the construction of the KDTree but
          slow down the search.
        """
        self.dim = 3  # 3D systems
        self.box = None
        self._dm = None  # matrix of central-cell vectors
        self._rm = None  # matrix of normalized reciprocal vectors
        self.initialize_bm(box)
        self.kdt = _CKDTree.KDTree(self.dim, bucket_size)
        self._indices = list()

    def initialize_bm(self, box):
        """
        Store box information and define direct and reciprocal box matrices.
        Rows of the direct matrix are the components of the central cell
        vectors. Rows of the reciprocal matrix are the components of the
        normalized reciprocal lattice vectors. Each row of the reciprocal
        matrix represents the vector normal to the unit cell face
        associated to each axis.
        For instance, in an orthorhombic cell the YZ-plane is associated to
        the X-axis and its normal vector is (1, 0, 0). In a triclinic cell,
        the plane associated to vector ``\vec{a}`` is perpendicular to the
        normalized cross product of ``\vec{b}`` and ``\vec{c}``.

        Parameters
        ----------
        box : array-like or ``None``, optional, default ``None``
          Simulation cell dimensions in the form of
          :attr:`MDAnalysis.trajectory.base.Timestep.dimensions` when
          periodic boundary conditions should be taken into account for
          the calculation of contacts.
        """
        box_type = _box_check(box)
        if box_type == 'ortho':
            a, b, c = box[:3]
            dm = np.array([[a, 0, 0],
                           [0, b, 0],
                           [0, 0, c]], dtype=np.float32)
            rm = np.array([[1, 0, 0],
                           [0, 1, 0],
                           [0, 0, 1]], dtype=np.float32)
        elif box_type in 'tri_box tri_vecs tri_vecs_bad':
            if box_type == 'tri_box':
                dm = triclinic_vectors(box)
            elif box_type == 'tri_vecs':
                dm = box.copy('C')
            else:  # case 'tri_vecs_bad'
                dm = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))
            rm = np.zeros(9, dtype=np.float32).reshape(3, 3)
            rm[0] = np.cross(dm[1], dm[2])
            rm[1] = np.cross(dm[2], dm[0])
            rm[2] = np.cross(dm[0], dm[1])
            for i in range(self.dim):
                rm[i] /= norm(rm[i])  # normalize
        else:
            raise ValueError('Failed to initialize direct/reciprocal matrices')
        self.box = box
        self._dm = dm
        self._rm = rm

    def set_coords(self, coords):
        """
        Add coordinates of the points. Wrapping of coordinates to the central
        cell is enforced along the periodic axes.

        Parameters
        ----------
        coords: NumPy.array
          Positions of points, shape=(N, 3) for N atoms.
        """
        _check_array(coords, 'coords')
        wrapped_data = apply_PBC(coords, self.box)
        self.kdt.set_data(wrapped_data)
        self.built = 1

    def find_centers(self, center_point, radius):
        """
        Find relevant images of a center point, inspired by
        https://github.com/patvarilly/periodic_kdtree

        Parameters
        ----------
        center_point: NumPy.array
          Coordinates of the query center point
        radius: float
          Maximum distance from center in search for neighbors

        Returns
        ------
        :class:`List`
          wrapped center point and its relevant images
        """
        wrapped_c = apply_PBC(center_point.reshape(1, 3), self.box)[0]
        # extents marks the max distance to plane just to consider an image
        extents = np.array([norm(self._dm[i])/2.0 for i in range(self.dim)])
        extents = np.where(extents < radius, extents, radius)
        # displacements are vectors that we add to wrapped_c to
        # generate images "up" or "down" the central cell along
        # the axis that we happen to be looking.
        displacements = list()
        for i in range(self.dim):
            # distance to the plane containing the origin
            if np.dot(wrapped_c, self._rm[i]) < extents[i]:
                displacements.append(self._dm[i])  # "upper" image
            # distance to the plane containing point self._dm[i]
            elif np.dot(self._dm[i] - wrapped_c, self._rm[i]) < extents[i]:
                displacements.append(-self._dm[i])  # "lower" image
        # If we have displacements along more than one axis, we have
        # to combine them. This happens when wrapped_center is close
        # to any edge or vertex of the central cell.
        # face, n_displacements==1; no combination
        # edge, n_displacements==2; combinations produce one extra displacement
        # vertex, n_displacements==3; five extra displacements
        n_displacements = len(displacements)
        if n_displacements > 1:
            for start in range(n_displacements - 1, -1, -1):
                for i in range(start+1, len(displacements)):
                    displacements.append(displacements[start]+displacements[i])
        return [wrapped_c, ] + [wrapped_c + d for d in displacements]

    def search(self, center, radius):
        """Search all points within radius of center and its periodic images.
        Wrapping of center coordinates is enforced to enable comparison to
        wrapped coordinates of points in the tree.

        Parameter
        ---------
        center: NumPy.array
          origin around which to search for neighbors
        radius: float
          maximum distance around which to search for neighbors. The search
          radius is half the smallest periodicity if radius exceeds this value
        """
        if not self.built:
            raise RuntimeError('Unbuilt tree. Run tree.set_coords first')
        if center.shape != (self.dim,):
            raise ValueError('Expected a ({},) NumPy array'.format(self.dim))
        self._indices = set()  # clear previous search
        for c in self.find_centers(center, radius):
            self.kdt.search_center_radius(c, radius)
            new_indices = self.kdt.get_indices()  # returns None or np.array
            if new_indices is not None:
                self._indices.update(new_indices)
        self._indices = sorted(list(self._indices))

    def get_indices(self):
        return self._indices
