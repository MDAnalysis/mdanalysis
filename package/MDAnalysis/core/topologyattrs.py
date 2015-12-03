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


class TopologyAttr(object):
    """Base class for Topology attributes.

    .. note::   This class is intended to be subclassed, and mostly amounts to a
                skeleton. The methods here should be present in all
                :class:`TopologyAttr` child classes, but by default they raise
                appropriate exceptions.

    """
    attrname = 'topologyattr'
    topology = None

    def __init__(self, values):
        """Generate a new :class:`TopologyAttr` object with the given values.

        """
        pass

    def get_atoms(self, idx):
        """Get atom attributes for given atom indices.

        """
        raise NoDataError

    def set_atoms(self, idx, values):
        """Set atom attributes for given atom indices.

        """
        raise NotImplementedError

    def get_residues(self, idx):
        """Get residue attributes for given residue indices.

        """
        raise NoDataError

    def set_residues(self, idx, values):
        """Set residue attributes for given residue indices.

        """
        raise NotImplementedError

    def get_segments(self, idx):
        """Get segment attributes for given segment indices.

        """
        raise NoDataError

    def set_segments(self, idx, values):
        """Set segmentattributes for given segment indices.

        """
        raise NotImplementedError


class Masses(TopologyAttr):
    """
    Interface to masses for atoms, residues, and segments.
    
    Parameters
    ----------
    masses : array
        mass for each atom in the system

    """
    def __init__(self, values):
        super(self, Masses).__init__(values)

        self._values = values

    def get_atoms(self, aix):
        self._values[aix]

    def set_atoms(self, aix, values):
        self.values[aix] = values

    def get_residues(self, rix):
        resatoms = self.topology.r2a(rix)



