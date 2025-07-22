# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""\
Core Topology object --- :mod:`MDAnalysis.core.topology`
========================================================

.. versionadded:: 0.16.0

:class:`Topology` is the core object that holds all topology information.

TODO: Add in-depth discussion.

Notes
-----
For developers: In MDAnalysis 0.16.0 this new topology system was
introduced and discussed as issue `#363`_; this issue contains key
information and discussions on the new system. The issue number *363*
is also being used as a short-hand in discussions to refer to the new
topology system.


.. _`#363`: https://github.com/MDAnalysis/mdanalysis/issues/363

Classes
-------

.. autoclass:: Topology
   :members:
.. autoclass:: TransTable
   :members:

Helper functions
----------------

.. autofunction:: make_downshift_arrays

"""
import contextlib

import numpy as np
import typing

from .topologyattrs import Atomindices, Resindices, Segindices
from ..exceptions import NoDataError


def make_downshift_arrays(upshift, nparents):
    """From an upwards translation table, create the opposite direction

    Turns a many to one mapping (eg atoms to residues) to a one to many mapping
    (residues to atoms)

    Parameters
    ----------
    upshift : array_like
        Array of integers describing which parent each item belongs to
    nparents : integer
        Total number of parents that exist.

    Returns
    -------
    downshift : array_like (dtype object)
        An array of arrays, each containing the indices of the children
        of each parent.  Length `nparents` + 1

    Examples
    --------

    To find the residue to atom mappings for a given atom to residue mapping:

    >>> import numpy as np
    >>> import MDAnalysis as mda
    >>> from MDAnalysis.core.topology import make_downshift_arrays
    >>> atom2res = np.array([0, 1, 0, 2, 2, 0, 2])
    >>> make_downshift_arrays(atom2res, 3)
    array([array([0, 2, 5]), array([1]), array([3, 4, 6]), None], dtype=object)

    Entry 0 corresponds to residue 0 and says that this contains atoms 0, 2 & 5

    Notes
    -----
    The final entry in the return array will be ``None`` to ensure that the
    dtype of the array is :class:`object`.

    .. warning:: This means negative indexing should **never**
                 be used with these arrays.
    """
    if not len(upshift):
        return np.array([], dtype=object)

    # mergesort for a stable ordered array for the same value.
    order = np.argsort(upshift, kind="mergesort")

    upshift_sorted = upshift[order]
    u_values, indices = np.unique(upshift_sorted, return_index=True)

    # reset nparents to the larger one between input and heuristic from data
    # This is useful for creating empty Universe where default value is 1.
    nparents = np.max([nparents, u_values.max() + 1])
    residue_indices = np.zeros(nparents, dtype=int)
    missing_resids = np.sort(np.setdiff1d(np.arange(nparents), u_values))
    indices = np.append(indices, upshift_sorted.shape[0])

    residue_indices[u_values] = indices[1:]

    for missing_resid in missing_resids:
        if missing_resid == 0:
            residue_indices[missing_resid] = 0
        else:
            residue_indices[missing_resid] = residue_indices[missing_resid - 1]

    downshift = np.split(order, residue_indices[:-1])
    # Add None to end of array to force it to be of type Object
    # Without this, a rectangular array gets squashed into a single array
    downshift.append(None)
    return np.array(downshift, dtype=object)


class TransTable(object):
    """Membership tables with methods to translate indices across levels.

    There are three levels; Atom, Residue and Segment.  Each Atom **must**
    belong in a Residue, each Residue **must** belong to a Segment.

    When translating upwards, eg finding which Segment a Residue belongs in,
    a single numpy array is returned.  When translating downwards, two options
    are available; a concatenated result (suffix `_1`) or a list for each parent
    object (suffix `_2d`).

    Parameters
    ----------
    n_atoms : int
        number of atoms in topology
    n_residues : int
        number of residues in topology
    n_segments : int
        number of segments in topology
    atom_resindex : 1-D array
        resindex for each atom in the topology; the number of unique values in
        this array must be <= `n_residues`, and the array must be length
        `n_atoms`; giving None defaults to placing all atoms in residue 0
    residue_segindex : 1-D array
        segindex for each residue in the topology; the number of unique values
        in this array must be <= `n_segments`, and the array must be length
        `n_residues`; giving None defaults to placing all residues in segment 0


    Attributes
    ----------
    n_atoms : int
        number of atoms in topology
    n_residues : int
        number of residues in topology
    n_segments : int
        number of segments in topology
    size : tuple
        tuple ``(n_atoms, n_residues, n_segments)`` describing the shape of
        the TransTable

    .. versionchanged:: 2.3.0
        Lazy building RA and SR.
    """

    def __init__(
        self,
        n_atoms,
        n_residues,
        n_segments,  # Size of tables
        atom_resindex=None,
        residue_segindex=None,  # Contents of tables
    ):
        self.n_atoms = n_atoms
        self.n_residues = n_residues
        self.n_segments = n_segments

        # built atom-to-residue mapping, and vice-versa
        if atom_resindex is None:
            self._AR = np.zeros(n_atoms, dtype=np.intp)
        else:
            self._AR = np.asarray(atom_resindex, dtype=np.intp).copy()
            if not len(self._AR) == n_atoms:
                raise ValueError("atom_resindex must be len n_atoms")
        self._RA = None

        # built residue-to-segment mapping, and vice-versa
        if residue_segindex is None:
            self._RS = np.zeros(n_residues, dtype=np.intp)
        else:
            self._RS = np.asarray(residue_segindex, dtype=np.intp).copy()
            if not len(self._RS) == n_residues:
                raise ValueError("residue_segindex must be len n_residues")
        self._SR = None

    def copy(self):
        """Return a deepcopy of this Transtable"""
        return self.__class__(
            self.n_atoms,
            self.n_residues,
            self.n_segments,
            atom_resindex=self._AR,
            residue_segindex=self._RS,
        )

    @property
    def RA(self):
        if self._RA is None:
            self._RA = make_downshift_arrays(self._AR, self.n_residues)
        return self._RA

    @property
    def SR(self):
        if self._SR is None:
            self._SR = make_downshift_arrays(self._RS, self.n_segments)
        return self._SR

    @property
    def size(self):
        """The shape of the table, ``(n_atoms, n_residues, n_segments)``.

        :meta private:
        """
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
        RA = self.RA
        try:
            return np.concatenate(RA[rix])
        except ValueError:  # rix is not iterable or empty
            # don't accidentally return a view!
            return RA[rix].astype(np.intp, copy=True)

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
        RA = self.RA
        try:
            return [RA[r].copy() for r in rix]
        except TypeError:
            return [RA[rix].copy()]  # why would this be singular for 2d?

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
        SR = self.SR
        try:
            return np.concatenate(SR[six])
        except ValueError:  # six is not iterable or empty
            # don't accidentally return a view!
            return SR[six].astype(np.intp, copy=True)

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
            `six`, with each element being an array of the residue indices
            present in that segment

        """
        SR = self.SR
        try:
            return [SR[s].copy() for s in six]
        except TypeError:
            return [SR[six].copy()]

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
        rix = self.segments2residues_1d(six)
        return self.residues2atoms_1d(rix)

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
        return [self.residues2atoms_1d(rix) for rix in rixs]

    # Move between different groups.
    def move_atom(self, aix, rix):
        """Move aix to be in rix"""
        self._AR[aix] = rix
        self._RA = None

    def move_residue(self, rix, six):
        """Move rix to be in six"""
        self._RS[rix] = six
        self._SR = None

    def add_Residue(self, segidx):
        # segidx - index of parent
        self.n_residues += 1
        self._RA = None
        self._RS = np.concatenate([self._RS, np.array([segidx])])
        self._SR = None

        return self.n_residues - 1

    def add_Segment(self):
        self.n_segments += 1
        self._SR = None
        return self.n_segments - 1

    def __getstate__(self):
        # don't serialize _RA and _SR for performance.
        attrs = self.__dict__
        attrs["_RA"] = None
        attrs["_SR"] = None
        return attrs


class Topology(object):
    """In-memory, array-based topology database.

    The topology model of MDanalysis features atoms, which must each be a
    member of one residue. Each residue, in turn, must be a member of one
    segment. The details of maintaining this heirarchy, and mappings of atoms
    to residues, residues to segments, and vice-versa, are handled internally
    by this object.

    """

    def __init__(
        self,
        n_atoms=1,
        n_res=1,
        n_seg=1,
        attrs=None,
        atom_resindex=None,
        residue_segindex=None,
    ):
        """
        Parameters
        ----------
        n_atoms : int
            number of atoms in topology. Must be larger then 1 at each level
        n_residues : int
            number of residues in topology. Must be larger then 1 at each level
        n_segments : int
            number of segments in topology. Must be larger then 1 at each level
        attrs : TopologyAttr objects
            components of the topology to be included
        atom_resindex : array
            1-D array giving the resindex of each atom in the system
        residue_segindex : array
            1-D array giving the segindex of each residue in the system

        """
        self.tt = TransTable(
            n_atoms,
            n_res,
            n_seg,
            atom_resindex=atom_resindex,
            residue_segindex=residue_segindex,
        )

        if attrs is None:
            attrs = []
        # add core TopologyAttrs that give access to indices
        attrs.extend((Atomindices(), Resindices(), Segindices()))

        # attach the TopologyAttrs
        self.attrs = []
        for topologyattr in attrs:
            self.add_TopologyAttr(topologyattr)

    def copy(self):
        """Return a deepcopy of this Topology"""
        new = self.__class__(1, 1, 1)
        # copy the tt
        new.tt = self.tt.copy()
        # remove indices
        for attr in self.attrs:
            if isinstance(attr, (Atomindices, Resindices, Segindices)):
                continue
            new.add_TopologyAttr(attr.copy())
        return new

    @property
    def n_atoms(self):
        return self.tt.n_atoms

    @property
    def n_residues(self):
        return self.tt.n_residues

    @property
    def n_segments(self):
        return self.tt.n_segments

    def add_TopologyAttr(self, topologyattr):
        """Add a new TopologyAttr to the Topology.

        Parameters
        ----------
        topologyattr : TopologyAttr

        """
        self.attrs.append(topologyattr)
        topologyattr.top = self
        self.__setattr__(topologyattr.attrname, topologyattr)

    def del_TopologyAttr(self, topologyattr):
        """Remove a TopologyAttr from the Topology.

        If it is not present, nothing happens.

        Parameters
        ----------
        topologyattr : TopologyAttr


        .. versionadded:: 2.0.0
        """
        self.__delattr__(topologyattr.attrname)
        self.attrs.remove(topologyattr)

    @property
    def guessed_attributes(self):
        """A list of the guessed attributes in this topology"""
        return filter(
            lambda x: (
                x.is_guessed
                if (not isinstance(x.is_guessed, typing.Container))
                else True in x.is_guessed
            ),
            self.attrs,
        )

    @property
    def read_attributes(self):
        """A list of the attributes read from the topology"""
        return filter(
            lambda x: (
                not x.is_guessed
                if (not isinstance(x.is_guessed, typing.Container))
                else False in x.is_guessed
            ),
            self.attrs,
        )

    def add_Residue(self, segment, **new_attrs):
        """
        Returns
        -------
        residx of the new Residue

        Raises
        ------
        NoDataError
          If not all data was provided.  This error is raised before any changes


        .. versionchanged:: 2.1.0
           Added use of _add_new to TopologyAttr resize
        """
        # Check that all data is here before making any changes
        for attr in self.attrs:
            if not attr.per_object == "residue":
                continue
            if attr.singular not in new_attrs:
                missing = (
                    attr.singular
                    for attr in self.attrs
                    if (
                        attr.per_object == "residue"
                        and attr.singular not in new_attrs
                    )
                )
                raise NoDataError(
                    "Missing the following attributes for the new"
                    " Residue: {}".format(", ".join(missing))
                )

        # Resize topology table
        residx = self.tt.add_Residue(segment.segindex)

        # Add new value to each attribute
        for attr in self.attrs:
            if not attr.per_object == "residue":
                continue
            newval = new_attrs[attr.singular]
            attr._add_new(newval)

        return residx

    def add_Segment(self, **new_attrs):
        """Adds a new Segment to the Topology

        Parameters
        ----------
        new_attrs : dict
          the new attributes for the new segment, eg {'segid': 'B'}

        Raises
        -------
        NoDataError
          if an attribute wasn't specified.

        Returns
        -------
        ix : int
          the idx of the new segment


        .. versionchanged:: 2.1.0
           Added use of _add_new to resize topology attrs
        """
        for attr in self.attrs:
            if attr.per_object == "segment":
                if attr.singular not in new_attrs:
                    missing = (
                        attr.singular
                        for attr in self.attrs
                        if (
                            attr.per_object == "segment"
                            and attr.singular not in new_attrs
                        )
                    )
                    raise NoDataError(
                        "Missing the following attributes for the"
                        " new Segment: {}"
                        "".format(", ".join(missing))
                    )

        segidx = self.tt.add_Segment()

        for attr in self.attrs:
            if not attr.per_object == "segment":
                continue
            newval = new_attrs[attr.singular]
            attr._add_new(newval)

        return segidx
