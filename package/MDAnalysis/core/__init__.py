# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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


"""Core functions of MDAnalysis
============================

The basic class is an :class:`~MDAnalysis.core.groups.AtomGroup`; the whole
simulation is called the
:class:`~MDAnalysis.core.universe.Universe`. Selections are computed on an
:class:`~MDAnalysis.core.groups.AtomGroup` and return another
:class:`~MDAnalysis.core.groups.AtomGroup`.

To get started, load the Universe::

  u = Universe(topology_file, trajectory_file)

A simple selection of all water oxygens within 4 A of the protein::

  water_shell = u.select_atoms('name OH2 and around 4.0 protein')
  water_shell.n_atoms           # how many waters were selected
  water_shell.total_mass()       # their total mass

:class:`~MDAnalysis.core.groups.AtomGroup` instances have various methods that
allow calculation of simple properties. For more complicated analysis, obtain
the coordinates as a numpy array ::

  coords = water_shell.positions

and write your own Python code.
"""

__all__ = ["AtomGroup", "Selection"]

from .groups import AtomGroup
from .selection import Selection
