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
        pass

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


class Atomids(TopologyAttr):
    """Interface to atomids.
    
    Parameters
    ----------
    atomids : array
        atomids for atoms in the system

    """
    attrname = 'atomids'
    level = 'atoms'

    def __init__(self, atomids):
        super(Atomids, self).__init__(atomids)
        self.values = atomids 

    def get_atoms(self, aix):
        return self.values[aix]

    def set_atoms(self, aix, atomids):
        self.values[aix] = atomids 

    def get_residues(self, rix):
        aix = self.top.tt.r2a_1d(rix)
        return self.values[aix]

    def get_segments(self, six):
        aix = self.top.tt.s2a_1d(six)
        return self.values[aix]


class Atomnames(TopologyAttr):
    """Interface to atomnames.
    
    Parameters
    ----------
    atomnames : array
        atomnames for atoms in the system

    """
    attrname = 'atomnames'
    level = 'atoms'

    def __init__(self, atomnames):
        super(Atomnames, self).__init__(atomnames)

        self.values = atomnames

    def get_atoms(self, aix):
        return self.values[aix]

    def set_atoms(self, aix, atomnames):
        self.values[aix] = atomnames

    def get_residues(self, rix):
        aix = self.top.tt.r2a_1d(rix)
        return self.values[aix]

    def get_segments(self, six):
        aix = self.top.tt.s2a_1d(six)
        return self.values[aix]


class Resids(TopologyAttr):
    """Interface to resids.
    
    Parameters
    ----------
    resids : array
        resids for residue in the system

    """
    attrname = 'resids'
    level = 'residues'

    def __init__(self, resids):
        super(Resids, self).__init__(resids)

        self.values = resids

    def get_atoms(self, aix):
        rix = self.top.tt.a2r(aix)
        return self.values[rix]

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

    def get_residues(self, rix):
        return self.values[rix]

    def set_residues(self, rix, resids):
        """Change resids of the given residues.

        """
        self.values[rix] = resids

    def get_segments(self, six):
        rix = self.top.tt.s2r_1d(six)
        return self.values[rix]


class Resnames(TopologyAttr):
    """Interface to resnames.
    
    Parameters
    ----------
    resnames : array
        resnames for residues in the system

    """
    attrname = 'resnames'
    level = 'residues'

    def __init__(self, resnames):
        super(Resids, self).__init__(resnames)

        self.values = resnames

    def get_atoms(self, aix):
        rix = self.top.tt.a2r(aix)
        return self.values[rix]

    def get_residues(self, rix):
        return self.values[rix]

    def set_residues(self, rix, resnames):
        """Change resnames for the given residues.

        """
        self.values[rix] = resnames


#TODO: need to add cacheing
class Masses(TopologyAttr):
    """Interface to masses for atoms, residues, and segments.
    
    Parameters
    ----------
    masses : array
        mass for each atom in the system

    """
    attrname = 'masses'
    level = 'atoms'

    def __init__(self, masses):
        super(Masses, self).__init__(masses)

        self.values = masses 

    def get_atoms(self, aix):
        self.values[aix]

    def set_atoms(self, aix, masses):
        self.values[aix] = masses 

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
