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
Topology object --- :mod:`MDAnalysis.core.topology'
===================================================================

"""
import numpy as np

from ..lib.mdamath import one_to_many_pointers
from .topologyattrs import Atomindices, Resindices, Segindices


class TransTable(object):
    """Membership tables with methods to translate indices across levels.

    Parameters
    ----------
    n_atoms, n_residues, n_segments : int
        number of atoms, residues, segments in topology
    atom_resindex : 1-D array
        resindex for each atom in the topology; the number of unique values in this
        array must be <= `n_residues`, and the array must be length `n_atoms`;
        giving None defaults to placing all atoms in residue 0
    residue_segindex : 1-D array
        segindex for each residue in the topology; the number of unique values in this
        array must be <= `n_segments`, and the array must be length `n_residues`;
        giving None defaults to placing all residues in segment 0
 

    Attributes
    ----------
    n_atoms, n_residues, n_segments : int
        number of atoms, residues, segments in topology
    AR : 1-D array
        resindex for each atom in the topology; allows fast upward translation
        from atoms to residues 
    RA : sparse matrix
        row ``i`` corresponds to the residue with resindex ``i``, with each
        column giving 1 if the atom with atomindex ``j`` is a member or 0 if it
        is not; this matrix has dimension (nres x natoms); allows fast downward
        translation from residues to atoms
    RS : 1-D array
        segindex for each residue in the topology; allows fast upward
        translation from residues to segments 
    SR : sparse matrix
        row ``i`` corresponds to the segment with segindex ``i``, with each
        column giving 1 if the residue with resindex ``j`` is a member or 0 if
        it is not; this matrix has dimension (nseg x nres); allows fast
        downward translation from segments to residues 

    """
    def __init__(self,
                 n_atoms, n_residues, n_segments,  # Size of tables
                 atom_resindex=None, residue_segindex=None,  # Contents of tables
                ):
        self.n_atoms = n_atoms
        self.n_residues = n_residues
        self.n_segments = n_segments

        # built atom-to-residue mapping, and vice-versa
        if atom_resindex is None:
            self.AR = np.zeros(n_atoms, dtype=np.int64)
        else:
            self.AR = atom_resindex
        self._atom_order, self._res_ptrs = one_to_many_pointers(
                n_atoms, n_residues, self.AR)

        # built residue-to-segment mapping, and vice-versa
        if residue_segindex is None:
            self.RS = np.zeros(n_residues, dtype=np.int64)
        else:
            self.RS = residue_segindex
        self._res_order, self._seg_ptrs = one_to_many_pointers(
                n_residues, n_segments, self.RS)

    @property
    def size(self):
        """The shape of the table, (n_atoms, n_residues, n_segments)"""
        return (self.n_atoms, self.n_residues, self.n_segments)

    def atoms2residues(self, aix):
        """Get residue indices for each atom.

        Parameters
        ----------
        aix : array
            atom indices

        Returns
        -------
        rix : array
            residue index for each atom 

        """
        return self.AR[aix]

    def residues2atoms_1d(self, rix):
        """Get atom indices collectively represented by given residue indices.

        Parameters
        ----------
        rix : array
            residue indices

        Returns
        -------
        aix : array
            indices of atoms present in residues, collectively

        """
        try:
            return np.concatenate([self._atom_order[x:y]
                                   for x, y in self._res_ptrs[rix]])
        except TypeError:
            x, y = self._res_ptrs[rix]
            return self._atom_order[x:y]

    def residues2atoms_2d(self, rix):
        """Get atom indices represented by each residue index.

        Parameters
        ----------
        rix : array
            residue indices

        Returns
        -------
        raix : tuple
            each element corresponds to a residue index, in order given in
            `rix`, with each element being an array of the atom indices present
            in that residue

        """
        try:
            return [self._atom_order[x:y] for x, y in self._res_ptrs[rix]]
        except TypeError:
            x, y = self._res_ptrs[rix]
            return self._atom_order[x:y]

    def residues2segments(self, rix):
        """Get segment indices for each residue.

        Parameters
        ----------
        rix : array
            residue indices 

        Returns
        -------
        six : array
            segment index for each residue

        """
        return self.RS[rix]

    def segments2residues_1d(self, six):
        """Get residue indices collectively represented by given segment indices

        Parameters
        ----------
        six : array
            segment indices

        Returns
        -------
        rix : array
            sorted indices of residues present in segments, collectively

        """
        try:
            return np.concatenate([self._res_order[x:y]
                                   for x, y in self._seg_ptrs[six]])
        except TypeError:
            x, y = self._seg_ptrs[six]
            return self._res_order[x:y]

    def segments2residues_2d(self, six):
        """Get residue indices represented by each segment index.

        Parameters
        ----------
        six : array
            residue indices

        Returns
        -------
        srix : sparse matrix 
            each element corresponds to a segment index, in order given in
            `six`, with each element being an array of the residue indices present
            in that segment

        """
        try:
            return [self._res_order[x:y] for x, y in self._seg_ptrs[six]]
        except TypeError:
            x, y = self._seg_ptrs[six]
            return self._res_order[x:y]

    # Compound moves, does 2 translations
    def atoms2segments(self, aix):
        """Get segment indices for each atom.

        Parameters
        ----------
        aix : array
            atom indices

        Returns
        -------
        rix : array
            segment index for each atom

        """
        rix = self.atoms2residues(aix)
        return self.residues2segments(rix)

    def segments2atoms_1d(self, six):
        """Get atom indices collectively represented by given segment indices.

        Parameters
        ----------
        six : array
            segment indices

        Returns
        -------
        aix : array
            sorted indices of atoms present in segments, collectively

        """

        rixs = self.segments2residues_2d(six)

        try:
            return np.concatenate([self.residues2atoms_1d(rix)
                                   for rix in rixs])
        except TypeError:
            return self.residues2atoms_1d(rixs)

    def segments2atoms_2d(self, six):
        """Get atom indices represented by each segment index.

        Parameters
        ----------
        six : array
            residue indices

        Returns
        -------
        saix : sparse matrix 
            each element corresponds to a segment index, in order given in
            `six`, with each element being an array of the atom indices present
            in that segment

        """
        rixs = self.segments2residues_2d(six)

        if isinstance(rixs, np.ndarray):
            return self.residues2atoms_1d(rixs)
        else:
            return (self.residues2atoms_1d(rix) for rix in rixs)


#TODO: movers and resizers

    # Move between different groups.
    # In general, delete old address, add new address
    def move_atom(self, aix, rix):
        """Move aix to be in rix"""
        pass

    def move_residue(self, rix, six):
        """Move rix to be in six"""
        pass

    def resize(self):
        pass


class Topology(object):
    """In-memory, array-based topology database.

    The topology model of MDanalysis features atoms, which must each be a
    member of one residue. Each residue, in turn, must be a member of one
    segment. The details of maintaining this heirarchy, and mappings of atoms
    to residues, residues to segments, and vice-versa, are handled internally
    by this object.

    Parameters
    ----------
    n_atoms, n_residues, n_segments : int
        number of atoms, residues, segments in topology; there must be at least
        1 element of each level in the system
    attrs : TopologyAttr objects
        components of the topology to be included
    atom_resindex : array
        1-D array giving the resindex of each atom in the system
    residue_segindex : array
        1-D array giving the segindex of each residue in the system
 
    """

    def __init__(self, n_atoms=1, n_res=1, n_seg=1,
                 attrs=None,
                 atom_resindex=None,
                 residue_segindex=None):
        self.n_atoms = n_atoms
        self.n_residues = n_res
        self.n_segments = n_seg
        if attrs is None:
            attrs = []
        self.tt = TransTable(n_atoms, n_res, n_seg,
                             atom_resindex=atom_resindex,
                             residue_segindex=residue_segindex)

        if attrs is None:
            attrs = []
        # add core TopologyAttrs that give access to indices
        attrs.extend((Atomindices(), Resindices(), Segindices()))

        # attach the TopologyAttrs
        self.attrs = []
        for topologyattr in attrs:
            self.add_TopologyAttr(topologyattr)

    def add_TopologyAttr(self, topologyattr):
        """Add a new TopologyAttr to the Topology.

        Parameters
        ----------
        topologyattr : TopologyAttr

        """
        self.attrs.append(topologyattr)
        topologyattr.top = self
        self.__setattr__(topologyattr.attrname, topologyattr)

    @property
    def guessed_attributes(self):
        """A list of the guessed attributes in this topology"""
        return filter(lambda x: x.is_guessed, self.attrs)

    @property
    def read_attributes(self):
        """A list of the attributes read from the topology"""
        return filter(lambda x: not x.is_guessed, self.attrs)
