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


class Resids(TopologyAttr):
    """Interface to resids.
    
    Parameters
    ----------
    resids : array
        resids for residue in the system

    """
    def __init__(self, resids):
        super(self, Resids).__init__(resids)

        self.resids = resids

    def get_atoms(self, aix):
        rix = self.topology.tt.a2r(aix)
        return self.resids[rix]

    def set_atoms(self, aix, resids):
        rix = np.zeros(len(aix), dtype=np.int64)

        # get resindexes for each resid
        for i, item in enumerate(resids):
            try:
                rix[i] = np.where(self.resids == item)[0][0]
            except IndexError:
                raise NoDataError("Cannot assign atom to a residue that doesn't already exist.")

        self.topology.tt.move_atom(aix, rix)

    def get_residues(self, rix):
        self.topology.tt.

    def get_segments(self, six):
        masses = np.empty(len(six))

        segatoms = self.topology.tt.s2sa(six)

        for i, row in enumerate(segatoms):
            masses[i] = self.values[row].sum()

        return masses


#TODO: need to add cacheing
class Masses(TopologyAttr):
    """Interface to masses for atoms, residues, and segments.
    
    Parameters
    ----------
    masses : array
        mass for each atom in the system

    """
    def __init__(self, masses):
        super(self, Masses).__init__(masses)

        self.masses = masses 

    def get_atoms(self, aix):
        self.masses[aix]

    def set_atoms(self, aix, masses):
        self.masses[aix] = masses 

    def get_residues(self, rix):
        masses = np.empty(len(rix))

        resatoms = self.topology.tt.r2ra(rix)

        for i, row in enumerate(resatoms):
            masses[i] = self.masses[row].sum()

        return masses

    def get_segments(self, six):
        masses = np.empty(len(six))

        segatoms = self.topology.tt.s2sa(six)

        for i, row in enumerate(segatoms):
            masses[i] = self.masses[row].sum()

        return masses
