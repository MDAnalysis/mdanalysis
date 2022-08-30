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
# doi: 10.25080/majora-629e541a-00e
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
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.lib.util import unique_int_1d
from MDAnalysis.core.groups import AtomGroup, SegmentGroup, ResidueGroup
import numpy.typing as npt
from typing import Optional, Union, List


class AtomNeighborSearch(object):
    """This class can be used to find all atoms/residues/segments within the
    radius of a given query position.

    For the neighbor search, this class is a wrapper around
    :class:`~MDAnalysis.lib.distances.capped_distance`.
    """

    def __init__(self, atom_group: AtomGroup,
                 box: Optional[npt.ArrayLike] = None) -> None:
        """

        Parameters
        ----------
        atom_list : AtomGroup
          list of atoms
        box : array-like or ``None``, optional, default ``None``
          Simulation cell dimensions in the form of
          :attr:`MDAnalysis.trajectory.timestep.Timestep.dimensions` when
          periodic boundary conditions should be taken into account for
          the calculation of contacts.
        """
        self.atom_group = atom_group
        self._u = atom_group.universe
        self._box = box

    def search(self, atoms: AtomGroup,
               radius: float,
               level: str = 'A'
               ) -> Optional[Union[AtomGroup, ResidueGroup, SegmentGroup]]:
        """
        Return all atoms/residues/segments that are within *radius* of the
        atoms in *atoms*.

        Parameters
        ----------
        atoms : AtomGroup, MDAnalysis.core.groups.AtomGroup
          AtomGroup object
        radius : float
          Radius for search in Angstrom.
        level : str
          char (A, R, S). Return atoms(A), residues(R) or segments(S) within
          *radius* of *atoms*.

        Returns
        -------
        AtomGroup : :class:`~MDAnalysis.core.groups.AtomGroup`
          When ``level='A'``, AtomGroup is being returned.
        ResidueGroup : :class:`~MDAnalysis.core.groups.ResidueGroup`
          When ``level='R'``, ResidueGroup is being returned.
        SegmentGroup : :class:`~MDAnalysis.core.groups.SegmentGroup`
          When ``level='S'``, SegmentGroup is being returned.


        .. versionchanged:: 2.0.0
           Now returns :class:`AtomGroup` (when empty this is now an empty
           :class:`AtomGroup` instead of an empty list), :class:`ResidueGroup`,
           or a :class:`SegmentGroup`
        """
        unique_idx = []
        try:
            # For atom groups, take the positions attribute
            position = atoms.positions
        except AttributeError:
            # For atom, take the position attribute
            position = atoms.position
        pairs = capped_distance(position, self.atom_group.positions,
                                radius, box=self._box, return_distances=False)

        if pairs.size > 0:
            unique_idx = unique_int_1d(np.asarray(pairs[:, 1], dtype=np.intp))
        return self._index2level(unique_idx, level)

    def _index2level(self,
                     indices: List[int],
                     level: str
                     ) -> Union[AtomGroup, ResidueGroup, SegmentGroup]:
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
        atomgroup = self.atom_group[indices]
        if level == 'A':
            return atomgroup
        elif level == 'R':
            return atomgroup.residues
        elif level == 'S':
            return atomgroup.segments
        else:
            raise NotImplementedError('{0}: level not implemented'.format(level))
