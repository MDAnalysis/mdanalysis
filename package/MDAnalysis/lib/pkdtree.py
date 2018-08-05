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

import itertools
import numpy as np
from scipy.spatial import cKDTree

from ._cutil import unique_int_1d
from ._augment import augment_coordinates, undo_augment
from .util import unique_rows

from MDAnalysis.lib.distances import apply_PBC

__all__ = [
    'PeriodicKDTree'
]


class PeriodicKDTree(object):
    """Wrapper around scipy.spatial.cKDtree

    Creates an object which can handle periodic as well as
    non periodic boundary condtions depending on the parameters
    provided while constructing the tree.

    To enable periodic boundary conditions, box dimensions must be
    provided. Periodic Boundary conditions are implemented by creating
    duplicates of the particles which are within the specified cutoff
    distance from the boundary. These duplicates along with the
    original particle coordinates are used with the cKDTree
    without any special treatment due to PBC beyond this point.
    The final results after any operation with duplicate particle indices
    can be traced back to the original particle using undo augment function.
    """
    def __init__(self, box=None, leafsize=10):
        """

        Parameters
        ----------
        box : array-like or ``None``, optional, default ``None``
          Simulation cell dimensions in the form of
          :attr:`MDAnalysis.trajectory.base.Timestep.dimensions` when
          periodic boundary conditions should be taken into account for
          the calculation of contacts.
        leafsize : int (optional)
          Number of entries in leafs of the KDTree. If you suffer poor
          performance you can play around with this number. Increasing the
          `leafsize` will speed up the construction of the KDTree but
          slow down the search.

        """
        self.leafsize = leafsize
        self.dim = 3  # 3D systems
        self.box = box
        self._built = False
        self.cutoff = None

    @property
    def pbc(self):
        return self.box is not None

    def set_coords(self, coords, cutoff=None):
        """Constructs KDTree from the coordinates

        Wrapping of coordinates to the primary unit cell is enforced
        before any distance evaluations. If periodic boundary conditions
        are enabled, then duplicate particles are generated in the
        vicinity of the box. An additional array `mapping` is also
        generated which can be later used to trace the origin of
        duplicate particle coordinates.

        For non-periodic calculations, cutoff should not be provided
        the parameter is only required for periodic calculations.

        Parameters
        ----------
        coords: array_like
          Coordinate array of shape ``(N, 3)`` for N atoms.
        cutoff: float
          Specified cutoff distance to create duplicate images
          Typically equivalent to the desired search radius
          or the maximum of the desired cutoff radius. Relevant images
          corresponding to every atom which lies
          within ``cutoff`` distance from either of the box boundary
          will be generated.


        See Also
        --------
        MDAnalysis.lib._augment.augment_coordinates


        """
        # If no cutoff distance is provided but PBC aware
        if self.pbc and (cutoff is None):
            raise RuntimeError('Provide a cutoff distance'
                               ' with tree.set_coords(...)')

        # set coords dtype to float32
        # augment coordinates will work only with float32
        coords = np.asarray(coords, dtype=np.float32)

        if self.pbc:
            self.cutoff = cutoff
            # Bring the coordinates in the central cell
            self.coords = apply_PBC(coords, self.box)
            # generate duplicate images
            self.aug, self.mapping = augment_coordinates(self.coords,
                                                         self.box,
                                                         self.cutoff)
            # Images + coords
            self.all_coords = np.concatenate([self.coords, self.aug])
            self.ckdt = cKDTree(self.all_coords, leafsize=self.leafsize)
        else:
            # if cutoff distance is provided for non PBC calculations
            if cutoff is not None:
                raise RuntimeError('Donot provide cutoff distance for'
                                   ' non PBC aware calculations')
            self.coords = coords
            self.ckdt = cKDTree(self.coords, self.leafsize)
        self._built = True

    def search(self, centers, radius):
        """Search all points within radius from centers and their periodic images.

        All the centers coordinates are wrapped around the central cell
        to enable distance evaluations from points in the tree
        and their images.

        Parameter
        ---------
        centers: array_like (N,3)
          coordinate array to search for neighbors
        radius: float
          maximum distance to search for neighbors.
        """

        if not self._built:
            raise RuntimeError('Unbuilt tree. Run tree.set_coords(...)')

        centers = np.asarray(centers)
        if centers.shape == (self.dim, ):
            centers = centers.reshape((1, self.dim))

        self._indices = set()  # clear previous search
        # Sanity check
        if self.pbc:
            if self.cutoff < radius:
                raise RuntimeError('Set cutoff greater or equal to the radius.')
            # Bring all query points to the central cell
            wrapped_centers = apply_PBC(centers, self.box)
            indices = list(self.ckdt.query_ball_point(wrapped_centers,
                                                      radius))
            self._indices = np.array(list(
                                     itertools.chain.from_iterable(indices)),
                                     dtype=np.int64)
            if self._indices.size > 0:
                self._indices = undo_augment(self._indices,
                                             self.mapping,
                                             len(self.coords))
        else:
            wrapped_centers = np.asarray(centers)
            indices = list(self.ckdt.query_ball_point(wrapped_centers,
                                                      radius))
            self._indices = np.array(list(
                                     itertools.chain.from_iterable(indices)),
                                     dtype=np.int64)
        self._indices = np.asarray(unique_int_1d(self._indices))
        return self._indices

    def get_indices(self):
        """
        Returns
        ------
        indices : list
          neighbors for the last query points and search radius
        """
        return self._indices

    def search_pairs(self, radius):
        """Search all the pairs within a specified radius

        Parameters
        ----------
        radius : float
          Maximum distance between pairs of coordinates

        Returns
        -------
        pairs : array
          Indices of all the pairs which are within the specified radius
        """
        if not self._built:
            raise RuntimeError(' Unbuilt Tree. Run tree.set_coords(...)')

        if self.pbc:
            if self.cutoff < radius:
                raise RuntimeError('Set cutoff greater or equal to the radius.')

        pairs = np.array(list(self.ckdt.query_pairs(radius)), dtype=np.int64)
        if self.pbc:
            if len(pairs) > 1:
                pairs[:, 0] = undo_augment(pairs[:, 0], self.mapping,
                                           len(self.coords))
                pairs[:, 1] = undo_augment(pairs[:, 1], self.mapping,
                                           len(self.coords))
        if pairs.size > 0:
            # First sort the pairs then pick the unique pairs
            pairs = np.sort(pairs, axis=1)
            pairs = unique_rows(pairs)
        return pairs
