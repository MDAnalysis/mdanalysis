# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
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
#


"""
A simple atomic distance analysis --- :mod:`MDAnalysis.analysis.atomicdistances`
==========================================================

:Author: Xu Hong Chen
:Year: 2023
:Copyright: GNU Public License v3

This module provides a class to efficiently compute distances between
two groups of atoms with an equal number of atoms over a trajectory.
Specifically, for two atom groups `ag1` and `ag2`, it will return 
`|ag1[i] - ag2[i]|` for all `i` from `0` to `n_atoms - 1`, where 
`n_atoms` is the number of atoms in each atom group. These distances
are grouped by timestep in a NumPy array.

For more general functions on computing distances between atoms or
groups of atoms, please see `MDAnalysis.analysis.distances`.

See Also
--------
:mod:`MDAnalysis.analysis.distances`
:mod:`MDAnalysis.lib.distances`
"""

import numpy as np

from MDAnalysis.lib.distances import calc_bonds


import warnings
import logging
from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.atomicdistances")


class AtomicDistances(AnalysisBase):
    r"""Class to calculate atomic distances between two AtomGroups over a
    trajectory.

    Parameters
    ----------
    ag1, ag2 : AtomGroup
        :class:`~MDAnalysis.core.groups.AtomGroup` with the
        same number of atoms
    pbc : bool, optional
        If ``True``, calculates atomic distances with periodic boundary
        conditions (PBCs). Setting `pbc` to ``False``, calculates atomic
        distances without considering PBCs. Defaults to ``True``.

    Attributes
    ----------
    results : :class:`numpy.ndarray`
        The distances `|ag1[i] - ag2[i]|` for all `i` from `0` to
        `n_atoms - 1` for each frame over the trajectory.
    n_frames : int
        Number of frames included in the analysis.
    n_atoms : int
        Number of atoms in each atom group.


    .. versionadded:: 2.5.0
    """

    def __init__(self, ag1, ag2, pbc=True, **kwargs):
        # check ag1 and ag2 have the same number of atoms
        if ag1.atoms.n_atoms != ag2.atoms.n_atoms:
            raise ValueError("AtomGroups do not "
                             "have the same number of atoms")
        # check ag1 and ag2 are from the same trajectory
        elif ag1.universe.trajectory != ag2.universe.trajectory:
            raise ValueError("AtomGroups are not "
                             "from the same trajectory")

        super(AtomicDistances, self).__init__(ag1.universe.trajectory,
                                              **kwargs)

        self._ag1 = ag1
        self._ag2 = ag2
        self._pbc = pbc

    def _prepare(self):
        # initialize NumPy array of frames x distances for results
        self.results = np.zeros((self.n_frames, self._ag1.atoms.n_atoms))

    def _single_frame(self):
        # if PBCs considered, get box size
        box = self._ag1.dimensions if self._pbc else None
        self.results[self._frame_index] = calc_bonds(self._ag1.positions,
                                                     self._ag2.positions,
                                                     box)
