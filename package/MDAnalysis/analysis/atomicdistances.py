# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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


r"""
Simple atomic distance analysis --- :mod:`MDAnalysis.analysis.atomicdistances`
==============================================================================

:Author: Xu Hong Chen
:Year: 2023
:Copyright: Lesser GNU Public License v2.1+

This module provides a class to efficiently compute distances between
two groups of atoms with an equal number of atoms over a trajectory.
Specifically, for two atom groups ``ag1`` and ``ag2``, it will return
the distances

.. math::
    |ag1[i] - ag2[i]|

for all :math:`i` from :math:`0` to `n_atoms` :math:`- 1`, where
`n_atoms` is the number of atoms in each atom group. By default,
this computation is done with periodic boundary conditions, but this
can be easily turned off. These distances are grouped by time step in
a NumPy array.

For more general functions on computing distances between atoms or
groups of atoms, please see :class:`MDAnalysis.analysis.distances`.

See Also
--------
:mod:`MDAnalysis.analysis.distances`
:mod:`MDAnalysis.lib.distances`


Basic usage
-----------

This example uses files from the MDAnalysis test suite
(:data:`~MDAnalysis.tests.datafiles.GRO` and
:data:`~MDAnalysis.tests.datafiles.XTC`). To get started, execute  ::

   >>> import MDAnalysis as mda
   >>> from MDAnalysis.tests.datafiles import GRO, XTC
   >>> import MDAnalysis.analysis.atomicdistances as ad

We will calculate the distances between an atom group of atoms 101-105
and an atom group of atoms 4001-4005 with periodic boundary conditions.
To select these atoms: ::

   >>> u = mda.Universe(GRO, XTC)
   >>> ag1 = u.atoms[100:105]
   >>> ag2 = u.atoms[4000:4005]

We can run the calculations using any variable of choice such as
``my_dists`` and access our results using ``my_dists.results``: ::

   >>> my_dists = ad.AtomicDistances(ag1, ag2).run()
   >>> my_dists.results
   array([[37.80813681, 33.2594864 , 34.93676414, 34.51183299, 34.96340209],
          [27.11746625, 31.19878079, 31.69439435, 32.63446126, 33.10451345],
          [23.27210749, 30.38714688, 32.48269361, 31.91444505, 31.84583838],
          [18.40607922, 39.21993135, 39.33468192, 41.0133789 , 39.46885946],
          [26.26006981, 37.9966713 , 39.14991106, 38.13423586, 38.95451427],
          [26.83845081, 34.66255735, 35.59335027, 34.8926705 , 34.27175056],
          [37.51994763, 38.12161091, 37.56481743, 36.8488121 , 35.75278065],
          [37.27275501, 37.7831456 , 35.74359073, 34.54893794, 34.76495816],
          [38.76272761, 41.31816555, 38.81588421, 39.82491432, 38.890219  ],
          [39.20012515, 40.00563374, 40.83857688, 38.77886735, 41.45775864]])

To do the computation without periodic boundary conditions, we can enter
the keyword argument ``pbc=False`` after ``ag2``. The result is different
in this case: ::

   >>> my_dists_nopbc = ad.AtomicDistances(ag1, ag2, pbc=False).run()
   >>> my_dists_nopbc.results
   array([[37.80813681, 33.2594864 , 34.93676414, 34.51183299, 34.96340209],
          [27.11746625, 31.19878079, 31.69439435, 32.63446126, 33.10451345],
          [23.27210749, 30.38714688, 32.482695  , 31.91444505, 31.84583838],
          [18.40607922, 39.21992825, 39.33468192, 41.0133757 , 39.46885946],
          [26.26006981, 37.99666906, 39.14990985, 38.13423708, 38.95451311],
          [26.83845081, 34.66255625, 35.59335027, 34.8926705 , 34.27174827],
          [51.86981409, 48.10347964, 48.39570072, 49.14423513, 50.44804292],
          [37.27275501, 37.7831456 , 35.74359073, 34.54893794, 34.76495816],
          [56.39657447, 41.31816555, 38.81588421, 39.82491432, 38.890219  ],
          [39.20012515, 40.00563374, 40.83857688, 38.77886735, 41.45775864]])

"""

import numpy as np

from MDAnalysis.lib.distances import calc_bonds

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
        The distances :math:`|ag1[i] - ag2[i]|` for all :math:`i`
        from :math:`0` to `n_atoms` :math:`- 1` for each frame over
        the trajectory.
    n_frames : int
        Number of frames included in the analysis.
    n_atoms : int
        Number of atoms in each atom group.


    .. versionadded:: 2.5.0
    """

    def __init__(self, ag1, ag2, pbc=True, **kwargs):
        # check ag1 and ag2 have the same number of atoms
        if ag1.atoms.n_atoms != ag2.atoms.n_atoms:
            raise ValueError(
                "AtomGroups do not " "have the same number of atoms"
            )
        # check ag1 and ag2 are from the same trajectory
        elif ag1.universe.trajectory != ag2.universe.trajectory:
            raise ValueError("AtomGroups are not " "from the same trajectory")

        super(AtomicDistances, self).__init__(
            ag1.universe.trajectory, **kwargs
        )

        self._ag1 = ag1
        self._ag2 = ag2
        self._pbc = pbc

    def _prepare(self):
        # initialize NumPy array of frames x distances for results
        self.results = np.zeros((self.n_frames, self._ag1.atoms.n_atoms))

    def _single_frame(self):
        # if PBCs considered, get box size
        box = self._ag1.dimensions if self._pbc else None
        self.results[self._frame_index] = calc_bonds(
            self._ag1.positions, self._ag2.positions, box
        )
