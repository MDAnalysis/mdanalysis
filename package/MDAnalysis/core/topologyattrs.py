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
    level = None

    def __init__(self, values):
        self.values = values

    def __len__(self):
        """Length of the TopologyAttr at its intrinsic level.

        If the level is 'atoms', then this should return a value of size
        atoms.

        """
        return len(self.values)

    def get_atoms(self, aix):
        """Get atom attributes for given atom indices.

        """
        raise NoDataError

    def set_atoms(self, aix, values):
        """Set atom attributes for given atom indices.

        """
        raise NotImplementedError

    def get_residues(self, rix):
        """Get residue attributes for given residue indices.

        """
        raise NoDataError

    def set_residues(self, rix, values):
        """Set residue attributes for given residue indices.

        """
        raise NotImplementedError

    def get_segments(self, six):
        """Get segment attributes for given segment indices.

        """
        raise NoDataError

    def set_segments(self, six, values):
        """Set segmentattributes for given segment indices.

        """
        raise NotImplementedError


## atom attributes

class AtomAttr(TopologyAttr):
    """Base class for atom attributes.

    """
    attrname = 'atomattr'
    topology = None
    level = 0 

    def get_atoms(self, aix):
        return self.values[aix]

    def set_atoms(self, aix, values):
        self.values[aix] = values

    def get_residues(self, rix):
        """By default, the values for each atom present in the set of residues
        are returned in a single array. This behavior can be overriden in child
        attributes.

        """
        aix = self.top.tt.r2a_1d(rix)
        return self.values[aix]

    def get_segments(self, six):
        """By default, the values for each atom present in the set of residues
        are returned in a single array. This behavior can be overriden in child
        attributes.

        """
        aix = self.top.tt.s2a_1d(six)
        return self.values[aix]


class Atomids(AtomAttr):
    """Interface to atomids.
    
    Parameters
    ----------
    atomids : array
        atomids for atoms in the system

    """
    attrname = 'atomids'


class Atomnames(AtomAttr):
    """Interface to atomnames.
    
    Parameters
    ----------
    atomnames : array
        atomnames for atoms in the system

    """
    attrname = 'atomnames'


#TODO: need to add cacheing
class Masses(AtomAttr):
    """Interface to masses for atoms, residues, and segments.
    
    Parameters
    ----------
    masses : array
        mass for each atom in the system

    """
    attrname = 'masses'

    def get_residues(self, rix):
        masses = np.empty(len(rix))

        resatoms = self.top.tt.r2ra(rix)

        for i, row in enumerate(resatoms):
            masses[i] = self.values[row].sum()

        return masses

    def get_segments(self, six):
        masses = np.empty(len(six))

        segatoms = self.top.tt.s2sa(six)

        for i, row in enumerate(segatoms):
            masses[i] = self.values[row].sum()

        return masses


#TODO: need to add cacheing
class Charges(AtomAttr):
    """Interface to charges for atoms, residues, and segments.
    
    Parameters
    ----------
    charges : array
        charge for each atom in the system

    """
    attrname = 'charges'

    def get_residues(self, rix):
        charges = np.empty(len(rix))

        resatoms = self.top.tt.r2a_2d(rix)

        for i, row in enumerate(resatoms):
            charges[i] = self.values[row].sum()

        return charges

    def get_segments(self, six):
        charges = np.empty(len(six))

        segatoms = self.top.tt.s2a_2d(six)

        for i, row in enumerate(segatoms):
            charges[i] = self.values[row].sum()

        return charges


## residue attributes

class ResidueAttr(TopologyAttr):
    """Base class for Topology attributes.

    .. note::   This class is intended to be subclassed, and mostly amounts to a
                skeleton. The methods here should be present in all
                :class:`TopologyAttr` child classes, but by default they raise
                appropriate exceptions.

    """
    attrname = 'residueattr'
    topology = None
    level = 1

    def get_atoms(self, aix):
        rix = self.top.tt.a2r(aix)
        return self.values[rix]

    def get_residues(self, rix):
        return self.values[rix]

    def set_residues(self, rix, values):
        self.values[rix] = values

    def get_segments(self, six):
        """By default, the values for each residue present in the set of
        segments are returned in a single array. This behavior can be overriden
        in child attributes.

        """
        rix = self.top.tt.s2r_1d(six)
        return self.values[rix]


class Resids(ResidueAttr):
    """Interface to resids.
    
    Parameters
    ----------
    resids : array
        resids for residue in the system

    """
    attrname = 'resids'

    def set_atoms(self, aix, resids):
        """Set resid for each atom given. Effectively moves each atom to
        another residue.

        """
        
        rix = np.zeros(len(aix), dtype=np.int64)

        # get resindexes for each resid
        for i, item in enumerate(resids):
            try:
                rix[i] = np.where(self.values == item)[0][0]
            except IndexError:
                raise NoDataError("Cannot assign atom to a residue that doesn't already exist.")

        self.top.tt.move_atom(aix, rix)


class Resnames(ResidueAttr):
    """Interface to resnames.
    
    Parameters
    ----------
    resnames : array
        resnames for residues in the system

    """
    attrname = 'resnames'


## segment attributes

class SegmentAttr(TopologyAttr):
    """Base class for segment attributes.

    """
    attrname = 'segmentattr'
    topology = None
    level = 2

    def get_atoms(self, aix):
        six = self.top.tt.a2s(aix)
        return self.values[six]

    def get_residues(self, rix):
        six = self.top.tt.r2s(rix)
        return self.values[six]

    def get_segments(self, six):
        return self.values[six]

    def set_segments(self, six, values):
        self.values[six] = values
