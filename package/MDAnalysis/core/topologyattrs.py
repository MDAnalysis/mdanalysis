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
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Topology attribute objects --- :mod:`MDAnalysis.core.topologyattrs'
===================================================================

"""

from MDAnalysis.exceptions import NoDataError


class Level(object):
    """Base class for levels represented in TopologyAttr.

    """
    def get(self, idx):
        """Get level attributes for given level indices.

        """

    def set(self, idx, values):
        """Set level attributes for given level indices.

        """

    def add(self, values):
        """Add level attributes; needed when new elements are added to the
        level.

        """

    def remove(self, idx):
        """Remove level attributes for given level indices; needed when
        elements are removed from the level.

        """


class TopologyAttr(object):
    """Base class for Topology attributes.

    :note: This class is intended to be subclassed, and mostly amounts to a
    skeleton.

    """
    attrname = 'topologyattr'

    def __init__(self, values):
        """Generate a new :class:`TopologyAttr` object with the given values.

        """
        self.atoms = Level()
        self.residues = Level()
        self.segments = Level()
