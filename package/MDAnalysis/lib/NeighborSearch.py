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
Neighbor Search wrapper for MDAnalysis --- :mod:`MDAnalysis.lib.NeighborSearch`
===============================================================================

This module contains classes that allow neighbor searches directly with
`AtomGroup` objects from `MDAnalysis`.
"""
from __future__ import absolute_import

import numpy as np
from MDAnalysis.lib.pkdtree import PeriodicKDTree
from MDAnalysis.lib.util import unique_int_1d

from MDAnalysis.core.groups import AtomGroup, Atom


class AtomNeighborSearch(object):
    """This class can be used to find all atoms/residues/segments within the
    radius of a given query position.

    For the neighbor search, this class uses the BioPython KDTree and its
    wrapper PeriodicKDTree for non-periodic and periodic systems, respectively.
    """

    def __init__(self, atom_group, box=None, bucket_size=10):
        """

        Parameters
        ----------
        atom_list : AtomGroup
          list of atoms
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
        self.atom_group = atom_group
        self._u = atom_group.universe
        self._box = box
        self.kdtree = PeriodicKDTree(box=box, leafsize=bucket_size)

    def search(self, atoms, radius, level='A'):
        """
        Return all atoms/residues/segments that are within *radius* of the
        atoms in *atoms*.

        Parameters
        ----------
        atoms : AtomGroup, MDAnalysis.core.groups.Atom
          list of atoms
        radius : float
          Radius for search in Angstrom.
        level : str
          char (A, R, S). Return atoms(A), residues(R) or segments(S) within
          *radius* of *atoms*.
        """
        if isinstance(atoms, Atom):
            positions = atoms.position.reshape(1, 3)
        else:
            positions = atoms.positions

        # check if already built
        
        cutoff = radius if self._box is not None else None
        self.kdtree.set_coords(self.atom_group.positions, cutoff=cutoff)

        indices = []
        for pos in positions:
            self.kdtree.search(pos, radius)
            indices.append(self.kdtree.get_indices())
        unique_idx = unique_int_1d(
            np.array([i for l in indices for i in l], dtype=np.int64)
        )
        return self._index2level(unique_idx, level)

    def _index2level(self, indices, level):
        """Convert list of atom_indices in a AtomGroup to either the
        Atoms or segments/residues containing these atoms.

        Parameters
        ----------
        indices
           list of atom indices
        level : str
          char (A, R, S). Return atoms(A), residues(R) or segments(S) within
          *radius* of *atoms*.
        """
        n_atom_list = self.atom_group[indices]
        if level == 'A':
            if not n_atom_list:
                return []
            else:
                return n_atom_list
        elif level == 'R':
            return list({a.residue for a in n_atom_list})
        elif level == 'S':
            return list(set([a.segment for a in n_atom_list]))
        else:
            raise NotImplementedError('{0}: level not implemented'.format(level))
