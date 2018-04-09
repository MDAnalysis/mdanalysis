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
PeriodicKDTree --- :mod:`MDAnalysis.lib.pkdtree`
===============================================================================

This module contains a class to allow searches on a KDTree involving periodic
boundary conditions.
"""

from __future__ import absolute_import
from six.moves import range

import itertools
import numpy as np
from Bio.KDTree import _CKDTree

from MDAnalysis.lib.distances import _box_check, _check_array, apply_PBC
from MDAnalysis.lib.mdamath import norm, triclinic_vectors, triclinic_box

__all__ = [
    'PeriodicKDTree',
]


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
            dm = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]], dtype=np.float32)
            rm = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float32)
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
        coords: array_like
          Positions of points, shape=(N, 3) for N atoms.
        """
        _check_array(np.asanyarray(coords), 'coords')
        wrapped_data = apply_PBC(coords, self.box)
        self.kdt.set_data(wrapped_data)
        self.built = 1

    def find_images(self, center_points, radius):
        """
        Find relevant images of a a set of center points, inspired by
        https://github.com/patvarilly/periodic_kdtree

        Parameters
        ----------
        center_points: array_like (N,3)
          Coordinates of the query center points. Must be in centrall cell.
        radius: float
          Maximum distance from center in search for neighbors

        Returns
        ------
        :class:`List`
          images of the center points close to any cell boundary plane
        """
        images = list()

        center_points = np.asanyarray(center_points)
        if center_points.shape == (3, ):
            center_points = center_points.reshape((1, 3))

        # Calculating displacement vectors for images of `center_point`
        # possible point for parallel loop version (Benchmark before!)
        for center_point in center_points:
            # distances to cell boundary planes passing through the origin
            distances = np.dot(self._rm, center_point)
            displacements = list(self._dm[np.where(distances < radius)[0]])
            # distances to the remaining three cell boundary planes
            distances = np.einsum('ij,ij->i', self._rm,
                                  self._dm - center_point)
            displacements.extend(
                list(-self._dm[np.where(distances < radius)[0]]))
            # If we have displacements along more than one axis, we have
            # to combine them. This happens when center_point is close
            # to any edge or vertex of the central cell.
            # face case: n_displacements==1; no combination
            # edge case: n_displacements==2; one extra displacement
            # vertex case: n_displacements==3; five extra displacements
            n_displacements = len(displacements)
            if n_displacements > 1:
                for start in range(n_displacements - 1, -1, -1):
                    for i in range(start+1, len(displacements)):
                        displacements.append(displacements[start]+displacements[i])
            images.extend([center_point + d for d in displacements])
        return images

    def search(self, centers, radius):
        """Search all points within radius of centers and its periodic images.
        Wrapping of center coordinates is enforced to enable comparison to
        wrapped coordinates of points in the tree.

        Parameter
        ---------
        centers: array_like (N,3)
          origins around which to search for neighbors
        radius: float
          maximum distance around which to search for neighbors. The search
          radius is half the smallest periodicity if radius exceeds this value
        """
        if not self.built:
            raise RuntimeError('Unbuilt tree. Run tree.set_coords(...) first')
        centers = np.asarray(centers)
        if centers.shape == (self.dim, ):
            centers = centers.reshape((1, self.dim))
        wrapped_centers = apply_PBC(centers, self.box)
        self._indices = set()  # clear previous search
        # possible loop for parallel execution (benchmark before!)
        for c in itertools.chain(wrapped_centers,
                                 self.find_images(wrapped_centers, radius)):
            self.kdt.search_center_radius(c, radius)
            n = self.kdt.get_count()
            if n:
                new_indices = np.empty(n, int)
                self.kdt.get_indices(new_indices)
                self._indices.update(new_indices)
        self._indices = sorted(list(self._indices))

    def get_indices(self):
        """
        Returns
        ------
        indices : list
          neighbors for the last query points and search radius
        """
        return self._indices
