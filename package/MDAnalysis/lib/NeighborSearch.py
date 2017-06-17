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
Neighbor Search wrapper for MDAnalysis --- :mod:`MDAnalysis.lib.NeighborSearch`
===============================================================================

This module contains classes that allow neighbor searches directly with
`AtomGroup` objects from `MDAnalysis`.
"""
from __future__ import absolute_import

import numpy as np
from Bio.KDTree import KDTree

from MDAnalysis.core.groups import AtomGroup, Atom

class AtomNeighborSearch(object):
    """This class can be used to find all atoms/residues/segements within the
    radius of a given query position.

    This class is using the BioPython KDTree for the neighborsearch. This class
    also does not apply PBC to the distance calculattions. So you have to ensure
    yourself that the trajectory has been corrected for PBC artifacts.
    """

    def __init__(self, atom_group, bucket_size=10):
        """

        Parameters
        ----------
        atom_list : AtomGroup
          list of atoms
        bucket_size : int
          Number of entries in leafs of the KDTree. If you suffer poor
          performance you can play around with this number. Increasing the
          `bucket_size` will speed up the construction of the KDTree but
          slow down the search.
        """
        self.atom_group = atom_group
        self._u = atom_group.universe
        self.kdtree = KDTree(dim=3, bucket_size=bucket_size)
        self.kdtree.set_coords(atom_group.positions)

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

        indices = []
        for pos in positions:
            self.kdtree.search(pos, radius)
            indices.append(self.kdtree.get_indices())
        unique_idx = np.unique([i for l in indices for i in l]).astype(np.int64)
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
