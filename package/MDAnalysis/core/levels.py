# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
Hierarchy levels --- :mod:`MDAnalysis.core.levels`
==================================================

Levels are used to define the relationships between different containers in
MDAnalysis

Levels can be navigated between using the .child and .parent attributes
Classes for containers at a given level are accessed using the .singular
and .plural attributes

                  Level:
        .singular =====================  .plural
Atom    <-        Atomlevel    ^ .child  ->  AtomGroup
Residue <-        Residuelevel !         ->  ResidueGroup
Segment <-        Segmentlevel v .parent ->  SegmentGroup
"""
import functools


_LEVEL_VALUES = {'atom': 1, 'residue': 2, 'segment': 3}


@functools.total_ordering
class Level(object):
    """Describes the level of hierarchy within MDA objects

    Can do comparisons with either

    .singular gives the Class for Components of this Level
    .plural gives the Class for Groups of this Level

    .parent gives the Level object above this
    .child gives the Level object below this
    """
    def __init__(self, name, singular, plural):
        self.name = name
        self.value = _LEVEL_VALUES[name]
        self.singular = singular
        self.plural = plural

    def __eq__(self, other):
        if isinstance(other, basestring):
            value = _LEVEL_VALUES[other]
        else:
            value = other.value
        return self.value == value

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if isinstance(other, basestring):
            value = _LEVEL_VALUES[other]
        else:
            value = other.value
        return self.value < value

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return self.name
