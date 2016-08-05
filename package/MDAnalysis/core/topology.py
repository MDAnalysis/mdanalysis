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

"""
# second docstring == multiline comment

TODO Notes:
  Could make downshift tables lazily built! This would
    a) Make these not get built when not used
    b) Optimise moving multiple atoms between residues as only built once afterwards

  Could optimise moves by only updating the two parent tables rather than rebuilding everything!


"""
from six.moves import zip
import numpy as np

from ..lib.mdamath import one_to_many_pointers
from .topologyattrs import Atomindices, Resindices, Segindices


def make_downshift_arrays(upshift):
    """From an upwards translation table, create the opposite direction

    Turns a many to one mapping to a one to many mapping

    Eg :: 
      AR = np.array([0, 1, 0, 2, 2, 0, 2])
    
      make_downshift_arrays(AR)
      array([array([0, 2, 5]), array([1]), array([3, 4, 6])], dtype=object)

    So entry 0 informs that items 0, 2 & 5 all belong to group 0.

    Returns
    -------
    Numpy array of numpy arrays (dtype object)
    Length : 
    """
    order = np.argsort(upshift)

    upshift_sorted = upshift[order]
    borders = [None] + list(np.nonzero(np.diff(upshift_sorted))[0] + 1) + [None]

    # returns an array of arrays
    return np.array([order[x:y]
                     for x, y in zip(borders[:-1], borders[1:])],
                    dtype=object)

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
    size
        tuple describing the shape of the TransTable

    Methods
    -------
    atoms2residues(aix)
        Returns the residue index for many atom indices
    residues2atoms_1d(rix)
        All atoms in the residues represented by *rix*
    residues2atoms_2d(rix)
        List of atom indices for each residue in *rix*
    residues2segments(rix)
        Segment indices for each residue in *rix*
    segments2residues_1d(six)
        Similar to `residues2atoms_1d`
    segments2residues_2d(six)
        Similar to `residues2atoms_2d`
    segments2atoms_1d(six)
        Similar to `residues2atoms_1d`
    segments2atoms_2d(six)
        Similar to `residues2atoms_2d`

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
            self._AR = np.zeros(n_atoms, dtype=np.int64)
        else:
            self._AR = atom_resindex.copy()
            if not len(self._AR) == n_atoms:
                raise ValueError("atom_resindex must be len n_atoms")
        self._RA = make_downshift_arrays(self._AR)

        # built residue-to-segment mapping, and vice-versa
        if residue_segindex is None:
            self._RS = np.zeros(n_residues, dtype=np.int64)
        else:
            self._RS = residue_segindex.copy()
            if not len(self._RS) == n_residues:
                raise ValueError("residue_segindex must be len n_residues")
        self._SR = make_downshift_arrays(self._RS)

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
        return self._AR[aix]

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
            return np.concatenate([self._RA[r] for r in rix])
        except TypeError:  # integers aren't iterable, raises TypeError
            return self._RA[rix].copy()  # don't accidentally send a view!

    def residues2atoms_2d(self, rix):
        """Get atom indices represented by each residue index.

        Parameters
        ----------
        rix : array
            residue indices

        Returns
        -------
        raix : list
            each element corresponds to a residue index, in order given in
            `rix`, with each element being an array of the atom indices present
            in that residue

        """
        try:
            return [self._RA[r].copy() for r in rix]
        except TypeError:
            return [self._RA[rix].copy()]  # why would this be singular for 2d?

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
        return self._RS[rix]

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
            return np.concatenate([self._SR[s] for s in six])
        except TypeError:
            return self._SR[six].copy()

    def segments2residues_2d(self, six):
        """Get residue indices represented by each segment index.

        Parameters
        ----------
        six : array
            residue indices

        Returns
        -------
        srix : list
            each element corresponds to a segment index, in order given in
            `six`, with each element being an array of the residue indices present
            in that segment

        """
        try:
            return [self._SR[s].copy() for s in six]
        except TypeError:
            return [self._SR[six].copy()]

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
        saix : list
            each element corresponds to a segment index, in order given in
            `six`, with each element being an array of the atom indices present
            in that segment

        """
        # residues in EACH 
        rixs = self.segments2residues_2d(six)

        if isinstance(rixs, np.ndarray):
            return self.residues2atoms_1d(rixs)
        else:
            return (self.residues2atoms_1d(rix) for rix in rixs)


    #TODO: movers and resizers

    # Move between different groups.
    def move_atom(self, aix, rix):
        """Move aix to be in rix"""
        self._AR[aix] = rix
        self._RA = make_downshift_arrays(self._AR)

    def move_residue(self, rix, six):
        """Move rix to be in six"""
        self._RS[rix] = six
        self._SR = make_downshift_arrays(self._RS)

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
