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
Neighbor Search wrapper for MDAnalysis --- :mod:`MDAnalysis.lib.NeighborSearch`
===============================================================================

This module contains classes that allow neighbor searches directly with
`AtomGroup` objects from `MDAnalysis`.
"""

import numpy as np
from Bio.KDTree import KDTree

from MDAnalysis.core.AtomGroup import AtomGroup

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
        self.kdtree = KDTree(dim=3, bucket_size=bucket_size)
        self.kdtree.set_coords(atom_group.positions)

    def search(self, atoms, radius, level='A'):
        """
        Return all atoms/residues/segments that are within *radius* of the
        atoms in *atoms*.

        Parameters
        ----------
        atoms : AtomGroup
          list of atoms
        radius : float
          Radius for search in Angstrom.
        level : str
          char (A, R, S). Return atoms(A), residues(R) or segments(S) within
          *radius* of *atoms*.
        """
        indices = []
        for atom in atoms.coordinates():
            self.kdtree.search(atom, radius)
            indices.append(self.kdtree.get_indices())
        unique_idx = np.unique([i for l in indices for i in l])
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
        n_atom_list = [self.atom_group[i] for i in indices]
        if level == 'A':
            if len(n_atom_list) == 0:
                return []
            else:
                return AtomGroup(n_atom_list)
        elif level == 'R':
            return list({a.residue for a in n_atom_list})
        elif level == 'S':
            return list(set([a.segment for a in n_atom_list]))
        else:
            raise NotImplementedError('{0}: level not implemented'.format(level))
