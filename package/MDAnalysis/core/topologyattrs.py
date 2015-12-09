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

Common TopologyAttrs used by most topology parsers.

"""
from collections import defaultdict
import itertools
import numpy as np

from ..lib.util import cached
from ..exceptions import NoDataError
from .topologyobjects import TopologyGroup


class TopologyAttr(object):
    """Base class for Topology attributes.

    .. note::   This class is intended to be subclassed, and mostly amounts to a
                skeleton. The methods here should be present in all
                :class:`TopologyAttr` child classes, but by default they raise
                appropriate exceptions.

    Attributes
    ----------
    attrname : str
        the name used for the attribute when attached to a ``Topology`` object
    top : Topology
        handle for the Topology object TopologyAttr is associated with
        
    """
    attrname = 'topologyattrs'
    singular = 'topologyattr'
    top = None

    def __init__(self, values):
        self.values = values

    def __len__(self):
        """Length of the TopologyAttr at its intrinsic level."""
        return len(self.values)

    def __getitem__(self, group):
        """Accepts an AtomGroup, ResidueGroup or SegmentGroup"""
        if group.level == 'atom':
            return self.get_atoms(group)
        elif group.level == 'residue':
            return self.get_residues(group)
        elif group.level == 'segment':
            return self.get_segments(group)

    def __setitem__(self, group, values):
        if group.level == 'atom':
            return self.set_atoms(group, values)
        elif group.level == 'residue':
            return self.set_residues(group, values)
        elif group.level == 'segment':
            return self.set_segments(group, values)

    def get_atoms(self, ag):
        """Get atom attributes for a given AtomGroup"""
        # aix = ag.indices
        raise NoDataError

    def set_atoms(self, ag, values):
        """Set atom attributes for a given AtomGroup"""
        raise NotImplementedError

    def get_residues(self, rg):
        """Get residue attributes for a given ResidueGroup"""
        raise NoDataError

    def set_residues(self, rg, values):
        """Set residue attributes for a given ResidueGroup"""
        raise NotImplementedError

    def get_segments(self, sg):
        """Get segment attributes for a given SegmentGroup"""
        raise NoDataError

    def set_segments(self, sg, values):
        """Set segmentattributes for a given SegmentGroup"""
        raise NotImplementedError


## core attributes

class Atomindices(TopologyAttr):
    """Globally unique indices for each atom in the group.

    If the group is an AtomGroup, then this gives the index for each atom in
    the group. This is the unambiguous identifier for each atom in the
    topology, and it is not alterable.

    If the group is a ResidueGroup or SegmentGroup, then this gives the indices
    of each atom represented in the group in a 1-D array, in the order of the
    elements in that group. 

    """
    attrname = 'indices'
    singular = 'index'

    def __init__(self):
        pass

    def set_atoms(self, ag, values):
        raise AttributeError("Atom indices are fixed; they cannot be reset")

    def get_atoms(self, ag):
        return ag._ix

    def get_residues(self, rg):
        return self.top.tt.r2a_1d(rg._ix)

    def get_segments(self, sg):
        return self.top.tt.s2a_1d(sg._ix)


class Resindices(TopologyAttr):
    """Globally unique resindices for each residue in the group.

    If the group is an AtomGroup, then this gives the resindex for each atom in
    the group. This unambiguously determines each atom's residue membership.
    Resetting these values changes the residue membership of the atoms.

    If the group is a ResidueGroup or SegmentGroup, then this gives the
    resindices of each residue represented in the group in a 1-D array, in the
    order of the elements in that group. 

    """
    attrname = 'resindices'
    singular = 'resindex'

    def __init__(self):
        pass

    def get_atoms(self, ag):
        return self.top.tt.a2r(ag._ix)

    def set_atoms(self, ag, values):
        """Set resindex for each atom given. Effectively moves each atom to
        another residue.

        """
        self.top.tt.move_atom(ag._ix, values)

    def get_residues(self, rg):
        return rg._ix

    def set_residues(self, rg, values):
        raise AttributeError("Residue indices are fixed; they cannot be reset")

    def get_segments(self, sg):
        return rix = self.top.tt.s2r_1d(sg._ix)


class Segindices(TopologyAttr):
    """Globally unique segindices for each segment in the group.

    If the group is an AtomGroup, then this gives the segindex for each atom in
    the group. This unambiguously determines each atom's segment membership.
    It is not possible to set these, since membership in a segment is an
    attribute of each atom's residue.
    
    If the group is a ResidueGroup or SegmentGroup, then this gives the
    segindices of each segment represented in the group in a 1-D array, in the
    order of the elements in that group. 

    """
    attrname = 'segindices'
    singular = 'segindex'

    def __init__(self):
        pass

    def get_atoms(self, ag):
        return self.top.tt.a2s(ag._ix)

    def get_residues(self, rg):
        return self.top.tt.r2s(rg._ix)

    def get_segments(self, sg):
        return sg._ix

    def set_segments(self, sg, values):
        raise AttributeError("Segment indices are fixed; they cannot be reset")


## atom attributes

class AtomAttr(TopologyAttr):
    """Base class for atom attributes.

    """
    attrname = 'atomattrs'
    singular = 'atomattr'

    def get_atoms(self, ag):
        return self.values[ag._ix]

    def set_atoms(self, ag, values):
        self.values[ag._ix] = values

    def get_residues(self, rg):
        """By default, the values for each atom present in the set of residues
        are returned in a single array. This behavior can be overriden in child
        attributes.

        """
        aix = self.top.tt.r2a_1d(rg._ix)
        return self.values[aix]

    def get_segments(self, sg):
        """By default, the values for each atom present in the set of residues
        are returned in a single array. This behavior can be overriden in child
        attributes.

        """
        aix = self.top.tt.s2a_1d(sg._ix)
        return self.values[aix]


#TODO: update docs to property doc
class Atomids(AtomAttr):
    """ID for each atom.
    """
    attrname = 'ids'
    singular = 'id'


#TODO: update docs to property doc
class Atomnames(AtomAttr):
    """Name for each atom.
    """
    attrname = 'names'
    singular = 'name'


#TODO: update docs to property doc
class Atomtypes(AtomAttr):
    """Type for each atom"""
    attrname = 'types'
    singular = 'type'


#TODO: update docs to property doc
class Radii(AtomAttr):
    """Radii for each atom"""
    attrname = 'radii'
    singular = 'radius'


#TODO: need to add cacheing
#TODO: update docs to property doc
class Masses(AtomAttr):
    attrname = 'masses'
    singular = 'mass'

    def get_residues(self, rg):
        masses = np.empty(len(rg))

        resatoms = self.top.tt.r2a_2d(rg._ix)

        for i, row in enumerate(resatoms):
            masses[i] = self.values[row].sum()

        return masses

    def get_segments(self, sg):
        masses = np.empty(len(sg))

        segatoms = self.top.tt.s2a_2d(sg._ix)

        for i, row in enumerate(segatoms):
            masses[i] = self.values[row].sum()

        return masses


#TODO: need to add cacheing
#TODO: update docs to property doc
class Charges(AtomAttr):
    attrname = 'charges'
    singular = 'charge'

    def get_residues(self, rg):
        charges = np.empty(len(rg))

        resatoms = self.top.tt.r2a_2d(rg._ix)

        for i, row in enumerate(resatoms):
            charges[i] = self.values[row].sum()

        return charges

    def get_segments(self, sg):
        charges = np.empty(len(sg))

        segatoms = self.top.tt.s2a_2d(sg._ix)

        for i, row in enumerate(segatoms):
            charges[i] = self.values[row].sum()

        return charges


#TODO: update docs to property doc
class Bfactors(AtomAttr):
    """Crystallographic B-factors in A**2 for each atom"""
    attrname = 'bfactors'
    singular = 'bfactors'


#TODO: update docs to property doc
class Occupancy(AtomAttr):
    attrname = 'occupancies'
    singular = 'occupancy'


#TODO: update docs to property doc
class AltLocs(AtomAttr):
    """AltLocs for each atom"""
    attrname = 'altLocs'
    singular = 'altLoc'


## residue attributes

class ResidueAttr(TopologyAttr):
    """Base class for Topology attributes.

    .. note::   This class is intended to be subclassed, and mostly amounts to a
                skeleton. The methods here should be present in all
                :class:`TopologyAttr` child classes, but by default they raise
                appropriate exceptions.

    """
    attrname = 'residueattrs'
    singular = 'residueattr'

    def get_atoms(self, ag):
        rix = self.top.tt.a2r(ag._ix)
        return self.values[rix]

    def get_residues(self, rg):
        return self.values[rg._ix]

    def set_residues(self, rg, values):
        self.values[rg._ix] = values

    def get_segments(self, sg):
        """By default, the values for each residue present in the set of
        segments are returned in a single array. This behavior can be overriden
        in child attributes.

        """
        rix = self.top.tt.s2r_1d(sg._ix)
        return self.values[rix]


#TODO: update docs to property doc
class Resids(ResidueAttr):
    attrname = 'resids'
    singular = 'resid'


#TODO: update docs to property doc
class Resnames(ResidueAttr):
    attrname = 'resnames'
    singular = 'resname'


#TODO: update docs to property doc
class Resnums(ResidueAttr):
    attrname = 'resnums'
    singular = 'resnum'


## segment attributes

class SegmentAttr(TopologyAttr):
    """Base class for segment attributes.

    """
    attrname = 'segmentattrs'
    singular = 'segmentattr'

    def get_atoms(self, ag):
        six = self.top.tt.a2s(ag._ix)
        return self.values[six]

    def get_residues(self, rg):
        six = self.top.tt.r2s(rg._ix)
        return self.values[six]

    def get_segments(self, sg):
        return self.values[sg._ix]

    def set_segments(self, sg, values):
        self.values[sg._ix] = values


#TODO: update docs to property doc
class Segids(SegmentAttr):
    attrname = 'segids'
    singular = 'segid'


#TODO: update docs to property doc
class Bonds(AtomAttr):
    """Bonds for atoms"""
    attrname = 'bonds'
    # Singular is the same because one Atom might have
    # many bonds, so still asks for "bonds" in the plural
    singular = 'bonds'

    def __init__(self, values):
        """
        Arguments
        ---------
        values - list of tuples of indices.  Should be zero based
        to match the atom indices

        Eg:  [(0, 1), (1, 2), (2, 3)]
        """
        self.values = values
        self._cache = dict()

    def __len__(self):
        return len(self._bondDict)

    @property
    @cached('bd')
    def _bondDict(self):
        """Lazily built mapping of atoms:bonds"""
        bd = defaultdict(list)

        for b in self.values:
            # We always want the first index
            # to be less than the last
            # eg (0, 1) not (1, 0)
            # and (4, 10, 8) not (8, 10, 4)
            if b[0] < b[-1]:
                b = b[::-1]
            for a in b:
                bd[a].append(b)
        return bd

    def get_atoms(self, ag):
        unique_bonds =  set(itertools.chain(
            *[self._bondDict[a] for a in ag._ix]))
        bond_idx = np.array(sorted(unique_bonds))
        return TopologyGroup(bond_idx, ag._u, self.singular[:-1])


#TODO: update docs to property doc
class Angles(Bonds):
    """Angles for atoms"""
    attrname = 'angles'
    singular = 'angles'


#TODO: update docs to property doc
class Dihedrals(Bonds):
    """Dihedrals for atoms"""
    attrname = 'dihedrals'
    singular = 'dihedrals'


#TODO: update docs to property doc
class Impropers(Bonds):
    attrname = 'impropers'
    singular = 'impropers'
