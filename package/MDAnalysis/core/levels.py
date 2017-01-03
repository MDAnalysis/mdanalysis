# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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

"""\
Hierarchy levels --- :mod:`MDAnalysis.core.levels`
==================================================

:class:`Level` instances are used to define the relationships between
different containers (see :mod:`MDAnalysis.core.groups`) in
MDAnalysis.

Levels can be navigated between using the :attr:`.child` and
:attr:`.parent` attributes of the containers. Classes for containers
at a given level are accessed using the :attr:`.singular` and
:attr:`.plural` attributes::

                     Level:
           .singular =====================  .plural
   Atom    <-        Atomlevel    ^ .child  ->  AtomGroup
   Residue <-        Residuelevel !         ->  ResidueGroup
   Segment <-        Segmentlevel v .parent ->  SegmentGroup

.. SeeAlso::
   :mod:`MDAnalysis.core.groups`

.. autoclass:: Level

"""

_LEVEL_VALUES = {'atom': 1, 'residue': 2, 'segment': 3}


class Level(object):
    """Describes the level of hierarchy within MDA objects

    Can do comparisons with either

    - :attr:`.singular` gives the Class for Components of this Level
    - :attr:`.plural` gives the Class for Groups of this Level

    or

    - :attr:`.parent` gives the :class:`Level` object above this
    - :attr:`.child` gives the :class:`Level` object below this

    """
    def __init__(self, name, singular, plural):
        self.name = name
        self.value = _LEVEL_VALUES[name]
        self.singular = singular
        self.plural = plural

    def __eq__(self, other):
        return self.value == other.value

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return self.name
