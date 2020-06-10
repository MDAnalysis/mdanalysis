# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
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
Topology attribute objects --- :mod:`MDAnalysis.core.topologyattrs`
===================================================================

Common :class:`TopologyAttr` instances that are used by most topology
parsers.

TopologyAttrs are used to contain attributes such as atom names or resids.
These are usually read by the TopologyParser.
"""
from __future__ import division, absolute_import
import six
from six.moves import zip, range

import Bio.Seq
import Bio.SeqRecord
from collections import defaultdict
import copy
import functools
import itertools
import numbers
import numpy as np
import warnings


from ..lib.util import (cached, convert_aa_code, iterable, warn_if_not_unique,
                        unique_int_1d)
from ..lib import transformations, mdamath
from ..exceptions import NoDataError, SelectionError
from .topologyobjects import TopologyGroup
from . import selection
from .groups import (ComponentBase, GroupBase,
                     Atom, Residue, Segment,
                     AtomGroup, ResidueGroup, SegmentGroup,
                     check_pbc_and_unwrap)
from .. import _TOPOLOGY_ATTRS, _TOPOLOGY_TRANSPLANTS, _TOPOLOGY_ATTRNAMES


def _check_length(func):
    """Wrapper which checks the length of inputs to set_X

    Eg:

    @_check_length
    def set_X(self, group, values):

    Will check the length of *values* compared to *group* before proceeding with
    anything in the *set_X* method.

    Pseudo code for the check:

    if group in (Atom, Residue, Segment):
        values must be single values, ie int, float or string
    else:
        values must be single value OR same length as group

    """
    _SINGLE_VALUE_ERROR = ("Setting {cls} {attrname} with wrong sized input. "
                           "Must use single value, length of supplied values: {lenvalues}.")
    # Eg "Setting Residue resid with wrong sized input. Must use single value, length of supplied
    # values: 2."

    _GROUP_VALUE_ERROR = ("Setting {group} {attrname} with wrong sized array. "
                          "Length {group}: {lengroup}, length of supplied values: {lenvalues}.")
    # Eg "Setting AtomGroup masses with wrong sized array. Length AtomGroup: 100, length of
    # supplied values: 50."

    def _attr_len(values):
        # quasi len measurement
        # strings, floats, ints are len 0, ie not iterable
        # other iterables are just len'd
        if iterable(values):
            return len(values)
        else:
            return 0  # special case

    @functools.wraps(func)
    def wrapper(attr, group, values):
        val_len = _attr_len(values)

        if isinstance(group, ComponentBase):
            if not val_len == 0:
                raise ValueError(_SINGLE_VALUE_ERROR.format(
                    cls=group.__class__.__name__, attrname=attr.singular,
                    lenvalues=val_len))
        else:
            if not (val_len == 0 or val_len == len(group)):
                raise ValueError(_GROUP_VALUE_ERROR.format(
                    group=group.__class__.__name__, attrname=attr.attrname,
                    lengroup=len(group), lenvalues=val_len))
        # if everything went OK, continue with the function
        return func(attr, group, values)

    return wrapper


def _wronglevel_error(attr, group):
    """Generate an error for setting attr at wrong level

    attr : TopologyAttr that was accessed
    group : Offending Component/Group

    Eg:
    setting mass of residue, gets called with attr=Masses, group=residue

    raises a NotImplementedError with:
    'Cannot set masses from Residue.  Use 'Residue.atoms.masses'

    Mainly used to ensure consistent and helpful error messages
    """
    if isinstance(group, (Atom, AtomGroup)):
        group_level = 1
    elif isinstance(group, (Residue, ResidueGroup)):
        group_level = 2
    elif isinstance(group, (Segment, SegmentGroup)):
        group_level = 3

    # What level to go to before trying to set this attr
    if isinstance(attr, AtomAttr):
        corr_classes = ('atoms', 'atom')
        attr_level = 1
    elif isinstance(attr, ResidueAttr):
        corr_classes = ('residues', 'residue')
        attr_level = 2
    elif isinstance(attr, SegmentAttr):
        corr_classes = ('segments', 'segment')
        attr_level = 3

    if isinstance(group, ComponentBase) and (attr_level > group_level):
        # ie going downards use plurals, going upwards use singulars
        # Residue.atom!s!.mass!es! but Atom.segment!!.segid!!
        correct = corr_classes[1]
        attrname = attr.singular
    else:
        correct = corr_classes[0]
        attrname = attr.attrname

    err_msg = "Cannot set {attr} from {cls}. Use '{cls}.{correct}.{attr} = '"
    # eg "Cannot set masses from Residue.  'Use Residue.atoms.masses = '"

    return NotImplementedError(err_msg.format(
        attr=attrname, cls=group.__class__.__name__, correct=correct,
    ))


class _TopologyAttrMeta(type):
    # register TopologyAttrs
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        attrname = classdict.get('attrname')
        singular = classdict.get('singular', attrname)

        if attrname is None:
            attrname = singular

        if singular:
            _TOPOLOGY_ATTRS[singular] = _TOPOLOGY_ATTRS[attrname] = cls
            _singular = singular.lower().replace('_', '')
            _attrname = attrname.lower().replace('_', '')
            _TOPOLOGY_ATTRNAMES[_singular] = singular
            _TOPOLOGY_ATTRNAMES[_attrname] = attrname

            for clstype, transplants in cls.transplants.items():
                for name, method in transplants:
                    _TOPOLOGY_TRANSPLANTS[name] = [attrname, method, clstype]
                    clean = name.lower().replace('_', '')
                    _TOPOLOGY_ATTRNAMES[clean] = name


class TopologyAttr(six.with_metaclass(_TopologyAttrMeta, object)):
    """Base class for Topology attributes.

    Note
    ----
    This class is intended to be subclassed, and mostly amounts to
    a skeleton. The methods here should be present in all
    :class:`TopologyAttr` child classes, but by default they raise
    appropriate exceptions.


    Attributes
    ----------
    attrname : str
        the name used for the attribute when attached to a ``Topology`` object
    singular : str
        name for the attribute on a singular object (Atom/Residue/Segment)
    per_object : str
        If there is a strict mapping between Component and Attribute
    dtype : int/float/object
        Type to coerce this attribute to be.  For string use 'object'
    top : Topology
        handle for the Topology object TopologyAttr is associated with

    """
    attrname = 'topologyattrs'
    singular = 'topologyattr'
    per_object = None  # ie Resids per_object = 'residue'
    top = None  # pointer to Topology object
    transplants = defaultdict(list)
    target_classes = []

    groupdoc = None
    singledoc = None
    dtype = None

    def __init__(self, values, guessed=False):
        if self.dtype is None:
            self.values = values
        else:
            self.values = np.asarray(values, dtype=self.dtype)
        self._guessed = guessed

    @staticmethod
    def _gen_initial_values(n_atoms, n_residues, n_segments):
        """Populate an initial empty data structure for this Attribute

        The only provided parameters are the "shape" of the Universe

        Eg for charges, provide np.zeros(n_atoms)
        """
        raise NotImplementedError("No default values")

    @classmethod
    def from_blank(cls, n_atoms=None, n_residues=None, n_segments=None,
                   values=None):
        """Create a blank version of this TopologyAttribute

        Parameters
        ----------
        n_atoms : int, optional
          Size of the TopologyAttribute atoms
        n_residues: int, optional
          Size of the TopologyAttribute residues
        n_segments : int, optional
          Size of the TopologyAttribute segments
        values : optional
          Initial values for the TopologyAttribute
        """
        if values is None:
            values = cls._gen_initial_values(n_atoms, n_residues, n_segments)
        elif cls.dtype is not None:
            # if supplied starting values and statically typed
            values = np.asarray(values, dtype=cls.dtype)
        return cls(values)

    def copy(self):
        """Return a deepcopy of this attribute"""
        return self.__class__(self.values.copy(), guessed=self._guessed)

    def __len__(self):
        """Length of the TopologyAttr at its intrinsic level."""
        return len(self.values)

    def __getitem__(self, group):
        """Accepts an AtomGroup, ResidueGroup or SegmentGroup"""
        if isinstance(group, (Atom, AtomGroup)):
            return self.get_atoms(group)
        elif isinstance(group, (Residue, ResidueGroup)):
            return self.get_residues(group)
        elif isinstance(group, (Segment, SegmentGroup)):
            return self.get_segments(group)

    def __setitem__(self, group, values):
        if isinstance(group, (Atom, AtomGroup)):
            return self.set_atoms(group, values)
        elif isinstance(group, (Residue, ResidueGroup)):
            return self.set_residues(group, values)
        elif isinstance(group, (Segment, SegmentGroup)):
            return self.set_segments(group, values)

    @property
    def is_guessed(self):
        """Bool of if the source of this information is a guess"""
        return self._guessed

    def get_atoms(self, ag):
        """Get atom attributes for a given AtomGroup"""
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


# core attributes

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
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup, Atom]

    def __init__(self):
        self._guessed = False

    def set_atoms(self, ag, values):
        raise AttributeError("Atom indices are fixed; they cannot be reset")

    def get_atoms(self, ag):
        return ag.ix

    def get_residues(self, rg):
        return list(self.top.tt.residues2atoms_2d(rg.ix))

    def get_segments(self, sg):
        return list(self.top.tt.segments2atoms_2d(sg.ix))


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
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup, Atom, Residue]

    def __init__(self):
        self._guessed = False

    def get_atoms(self, ag):
        return self.top.tt.atoms2residues(ag.ix)

    def get_residues(self, rg):
        return rg.ix

    def set_residues(self, rg, values):
        raise AttributeError("Residue indices are fixed; they cannot be reset")

    def get_segments(self, sg):
        return list(self.top.tt.segments2residues_2d(sg.ix))


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
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup,
                      Atom, Residue, Segment]

    def __init__(self):
        self._guessed = False

    def get_atoms(self, ag):
        return self.top.tt.atoms2segments(ag.ix)

    def get_residues(self, rg):
        return self.top.tt.residues2segments(rg.ix)

    def get_segments(self, sg):
        return sg.ix

    def set_segments(self, sg, values):
        raise AttributeError("Segment indices are fixed; they cannot be reset")


# atom attributes

class AtomAttr(TopologyAttr):
    """Base class for atom attributes.

    """
    attrname = 'atomattrs'
    singular = 'atomattr'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup, Atom]

    def get_atoms(self, ag):
        return self.values[ag.ix]

    @_check_length
    def set_atoms(self, ag, values):
        self.values[ag.ix] = values

    def get_residues(self, rg):
        """By default, the values for each atom present in the set of residues
        are returned in a single array. This behavior can be overriden in child
        attributes.

        """
        aixs = self.top.tt.residues2atoms_2d(rg.ix)
        return [self.values[aix] for aix in aixs]

    def set_residues(self, rg, values):
        raise _wronglevel_error(self, rg)

    def get_segments(self, sg):
        """By default, the values for each atom present in the set of residues
        are returned in a single array. This behavior can be overriden in child
        attributes.

        """
        aixs = self.top.tt.segments2atoms_2d(sg.ix)
        return [self.values[aix] for aix in aixs]

    def set_segments(self, sg, values):
        raise _wronglevel_error(self, sg)


# TODO: update docs to property doc
class Atomids(AtomAttr):
    """ID for each atom.
    """
    attrname = 'ids'
    singular = 'id'
    per_object = 'atom'
    dtype = int

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.arange(1, na + 1)


# TODO: update docs to property doc
class Atomnames(AtomAttr):
    """Name for each atom.
    """
    attrname = 'names'
    singular = 'name'
    per_object = 'atom'
    dtype = object
    transplants = defaultdict(list)

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(na)], dtype=object)

    def phi_selection(residue, c_name='C', n_name='N', ca_name='CA'):
        """Select AtomGroup corresponding to the phi protein backbone dihedral
        C'-N-CA-C.

        Parameters
        ----------
        c_name: str (optional)
            name for the backbone C atom
        n_name: str (optional)
            name for the backbone N atom
        ca_name: str (optional)
            name for the alpha-carbon atom

        Returns
        -------
        AtomGroup
            4-atom selection in the correct order. If no C' found in the
            previous residue (by resid) then this method returns ``None``.

        .. versionchanged:: 1.0.0
            Added arguments for flexible atom names and refactored code for
            faster atom matching with boolean arrays.
        """
        # fnmatch is expensive. try the obv candidate first
        prev = residue.universe.residues[residue.ix-1]
        sid = residue.segment.segid
        rid = residue.resid-1
        if not (prev.segment.segid == sid and prev.resid == rid):
            sel = 'segid {} and resid {}'.format(sid, rid)
            try:
                prev = residue.universe.select_atoms(sel).residues[0]
            except IndexError:
                return None
        c_ = prev.atoms[prev.atoms.names == c_name]
        if not len(c_) == 1:
            return None

        atnames = residue.atoms.names
        ncac_names = [n_name, ca_name, c_name]
        ncac = [residue.atoms[atnames == n] for n in ncac_names]
        if not all(len(ag) == 1 for ag in ncac):
            return None

        sel = c_+sum(ncac)
        return sel

    transplants[Residue].append(('phi_selection', phi_selection))

    def phi_selections(residues, c_name='C', n_name='N', ca_name='CA'):
        """Select list of AtomGroups corresponding to the phi protein 
        backbone dihedral C'-N-CA-C.

        Parameters
        ----------
        c_name: str (optional)
            name for the backbone C atom
        n_name: str (optional)
            name for the backbone N atom
        ca_name: str (optional)
            name for the alpha-carbon atom

        Returns
        -------
        list of AtomGroups
            4-atom selections in the correct order. If no C' found in the
            previous residue (by resid) then corresponding item in the list
            is ``None``.

        .. versionadded:: 1.0.0
        """

        u = residues[0].universe
        prev = u.residues[residues.ix-1]  # obv candidates first
        rsid = residues.segids
        prid = residues.resids-1
        ncac_names = [n_name, ca_name, c_name]
        sel = 'segid {} and resid {}'

        # replace wrong residues
        wix = np.where((prev.segids != rsid) | (prev.resids != prid))[0]
        invalid = []
        if len(wix):
            prevls = list(prev)
            for s, r, i in zip(rsid[wix], prid[wix], wix):
                try:
                    prevls[i] = u.select_atoms(sel.format(s, r)).residues[0]
                except IndexError:
                    invalid.append(i)
            prev = sum(prevls)
        
        keep_prev = [sum(r.atoms.names == c_name) == 1 for r in prev]
        keep_res = [all(sum(r.atoms.names == n) == 1 for n in ncac_names) 
                    for r in residues]
        keep = np.array(keep_prev) & np.array(keep_res)
        keep[invalid] = False
        results = np.zeros_like(residues, dtype=object)
        results[~keep] = None
        prev = prev[keep]
        residues = residues[keep]
        keepix = np.where(keep)[0]
        
        c_ = prev.atoms[prev.atoms.names == c_name]
        n = residues.atoms[residues.atoms.names == n_name]
        ca = residues.atoms[residues.atoms.names == ca_name]
        c = residues.atoms[residues.atoms.names == c_name]
        results[keepix] = [sum(atoms) for atoms in zip(c_, n, ca, c)]
        return list(results)

    transplants[ResidueGroup].append(('phi_selections', phi_selections))


    def psi_selection(residue, c_name='C', n_name='N', ca_name='CA'):
        """Select AtomGroup corresponding to the psi protein backbone dihedral
        N-CA-C-N'.

        Parameters
        ----------
        c_name: str (optional)
            name for the backbone C atom
        n_name: str (optional)
            name for the backbone N atom
        ca_name: str (optional)
            name for the alpha-carbon atom

        Returns
        -------
        AtomGroup
            4-atom selection in the correct order. If no N' found in the
            following residue (by resid) then this method returns ``None``.

        .. versionchanged:: 1.0.0
            Added arguments for flexible atom names and refactored code for
            faster atom matching with boolean arrays.
        """

        # fnmatch is expensive. try the obv candidate first
        _manual_sel = False
        sid = residue.segment.segid
        rid = residue.resid+1
        try:
            nxt = residue.universe.residues[residue.ix+1]
        except IndexError:
            _manual_sel = True
        else:
            if not (nxt.segment.segid == sid and nxt.resid == rid):
                _manual_sel = True

        if _manual_sel:
            sel = 'segid {} and resid {}'.format(sid, rid)
            try:
                nxt = residue.universe.select_atoms(sel).residues[0]
            except IndexError:
                return None
        n_ = nxt.atoms[nxt.atoms.names == n_name]
        if not len(n_) == 1:
            return None

        atnames = residue.atoms.names
        ncac_names = [n_name, ca_name, c_name]
        ncac = [residue.atoms[atnames == n] for n in ncac_names]
        if not all(len(ag) == 1 for ag in ncac):
            return None

        sel = sum(ncac) + n_
        return sel

    transplants[Residue].append(('psi_selection', psi_selection))

    def _get_next_residues_by_resid(residues):
        """Select list of Residues corresponding to the next resid for each 
        residue in `residues`.

        Returns
        -------
        List of Residues
            List of the next residues in the Universe, by resid and segid. 
            If not found, the corresponding item in the list is ``None``.

        .. versionadded:: 1.0.0
        """
        u = residues[0].universe
        nxres = np.array([None]*len(residues), dtype=object)
        ix = np.arange(len(residues))
        # no guarantee residues is ordered or unique
        last = max(residues.ix)
        if last == len(u.residues)-1:
            notlast = residues.ix != last
            ix = ix[notlast]
            residues = residues[notlast]

        nxres[ix] = nxt = u.residues[residues.ix+1]
        rsid = residues.segids
        nrid = residues.resids+1
        sel = 'segid {} and resid {}'

        # replace wrong residues
        wix = np.where((nxt.segids != rsid) | (nxt.resids != nrid))[0]
        if len(wix):
            for s, r, i in zip(rsid[wix], nrid[wix], wix):
                try:
                    nxres[ix[i]] = u.select_atoms(sel.format(s, r)).residues[0]
                except IndexError:
                    nxres[ix[i]] = None
        return nxres

    transplants[ResidueGroup].append(('_get_next_residues_by_resid',
                                      _get_next_residues_by_resid))

    def _get_prev_residues_by_resid(residues):
        """Select list of Residues corresponding to the previous resid for each 
        residue in `residues`.

        Returns
        -------
        List of Residues
            List of the previous residues in the Universe, by resid and segid. 
            If not found, the corresponding item in the list is ``None``.

        .. versionadded:: 1.0.0
        """
        u = residues[0].universe
        pvres = np.array([None]*len(residues))
        pvres[:] = prev = u.residues[residues.ix-1]
        rsid = residues.segids
        prid = residues.resids-1
        sel = 'segid {} and resid {}'

        # replace wrong residues
        wix = np.where((prev.segids != rsid) | (prev.resids != prid))[0]
        if len(wix):
            for s, r, i in zip(rsid[wix], prid[wix], wix):
                try:
                    pvres[i] = u.select_atoms(sel.format(s, r)).residues[0]
                except IndexError:
                    pvres[i] = None
        return pvres

    
    transplants[ResidueGroup].append(('_get_prev_residues_by_resid',
                                      _get_prev_residues_by_resid))

    def psi_selections(residues, c_name='C', n_name='N', ca_name='CA'):
        """Select list of AtomGroups corresponding to the psi protein 
        backbone dihedral N-CA-C-N'.

        Parameters
        ----------
        c_name: str (optional)
            name for the backbone C atom
        n_name: str (optional)
            name for the backbone N atom
        ca_name: str (optional)
            name for the alpha-carbon atom

        Returns
        -------
        List of AtomGroups
            4-atom selections in the correct order. If no N' found in the
            following residue (by resid) then the corresponding item in the 
            list is ``None``.

        .. versionadded:: 1.0.0
        """
        results = np.array([None]*len(residues), dtype=object)
        nxtres = residues._get_next_residues_by_resid()
        rix = np.where(nxtres)[0]
        nxt = sum(nxtres[rix])
        residues = residues[rix]
        ncac_names = [n_name, ca_name, c_name]

        keep_nxt = [sum(r.atoms.names == n_name) == 1 for r in nxt]
        keep_res = [all(sum(r.atoms.names == n) == 1 for n in ncac_names) 
                    for r in residues]
        keep = np.array(keep_nxt) & np.array(keep_res)
        nxt = nxt[keep]
        residues = residues[keep]
        keepix = np.where(keep)[0]
        
        n = residues.atoms[residues.atoms.names == n_name]
        ca = residues.atoms[residues.atoms.names == ca_name]
        c = residues.atoms[residues.atoms.names == c_name]
        n_ = nxt.atoms[nxt.atoms.names == n_name]
        results[rix[keepix]] = [sum(atoms) for atoms in zip(n, ca, c, n_)]
        return list(results)

    transplants[ResidueGroup].append(('psi_selections', psi_selections))

    def omega_selection(residue, c_name='C', n_name='N', ca_name='CA'):
        """Select AtomGroup corresponding to the omega protein backbone dihedral
        CA-C-N'-CA'.

        omega describes the -C-N- peptide bond. Typically, it is trans (180
        degrees) although cis-bonds (0 degrees) are also occasionally observed
        (especially near Proline).

        Parameters
        ----------
        c_name: str (optional)
            name for the backbone C atom
        n_name: str (optional)
            name for the backbone N atom
        ca_name: str (optional)
            name for the alpha-carbon atom

        Returns
        -------
        AtomGroup
            4-atom selection in the correct order. If no C' found in the
            previous residue (by resid) then this method returns ``None``.

        .. versionchanged:: 1.0.0
            Added arguments for flexible atom names and refactored code for
            faster atom matching with boolean arrays.
        """
        # fnmatch is expensive. try the obv candidate first
        _manual_sel = False
        sid = residue.segment.segid
        rid = residue.resid+1
        try:
            nxt = residue.universe.residues[residue.ix+1]
        except IndexError:
            _manual_sel = True
        else:
            if not (nxt.segment.segid == sid and nxt.resid == rid):
                _manual_sel = True

        if _manual_sel:
            sel = 'segid {} and resid {}'.format(sid, rid)
            try:
                nxt = residue.universe.select_atoms(sel).residues[0]
            except IndexError:
                return None

        ca = residue.atoms[residue.atoms.names == ca_name]
        c = residue.atoms[residue.atoms.names == c_name]
        n_ = nxt.atoms[nxt.atoms.names == n_name]
        ca_ = nxt.atoms[nxt.atoms.names == ca_name]

        if not all(len(ag) == 1 for ag in [ca_, n_, ca, c]):
            return None

        return ca+c+n_+ca_

    transplants[Residue].append(('omega_selection', omega_selection))

    def omega_selections(residues, c_name='C', n_name='N', ca_name='CA'):
        """Select list of AtomGroups corresponding to the omega protein 
        backbone dihedral CA-C-N'-CA'.

        omega describes the -C-N- peptide bond. Typically, it is trans (180
        degrees) although cis-bonds (0 degrees) are also occasionally observed
        (especially near Proline).

        Parameters
        ----------
        c_name: str (optional)
            name for the backbone C atom
        n_name: str (optional)
            name for the backbone N atom
        ca_name: str (optional)
            name for the alpha-carbon atom

        Returns
        -------
        List of AtomGroups
            4-atom selections in the correct order. If no C' found in the
            previous residue (by resid) then the corresponding item in the 
            list is ``None``.

        .. versionadded:: 1.0.0
        """
        results = np.array([None]*len(residues), dtype=object)
        nxtres = residues._get_next_residues_by_resid()
        rix = np.where(nxtres)[0]
        nxt = sum(nxtres[rix])
        residues = residues[rix]
        
        nxtatoms = [ca_name, n_name]
        resatoms = [ca_name, c_name]
        keep_nxt = [all(sum(r.atoms.names == n) == 1 for n in nxtatoms) 
                    for r in nxt]
        keep_res = [all(sum(r.atoms.names == n) == 1 for n in resatoms) 
                    for r in residues]
        keep = np.array(keep_nxt) & np.array(keep_res)
        nxt = nxt[keep]
        residues = residues[keep]
        keepix = np.where(keep)[0]
        
        c = residues.atoms[residues.atoms.names == c_name]
        ca = residues.atoms[residues.atoms.names == ca_name]
        n_ = nxt.atoms[nxt.atoms.names == n_name]
        ca_ = nxt.atoms[nxt.atoms.names == ca_name]

        results[rix[keepix]] = [sum(atoms) for atoms in zip(ca, c, n_, ca_)]
        return list(results)

    transplants[ResidueGroup].append(('omega_selections', omega_selections))

    def chi1_selection(residue, n_name='N', ca_name='CA', cb_name='CB', 
                       cg_name='CG'):
        """Select AtomGroup corresponding to the chi1 sidechain dihedral N-CA-CB-CG.

        Parameters
        ----------
        c_name: str (optional)
            name for the backbone C atom
        ca_name: str (optional)
            name for the alpha-carbon atom
        cb_name: str (optional)
            name for the beta-carbon atom
        cg_name: str (optional)
            name for the gamma-carbon atom

        Returns
        -------
        AtomGroup
            4-atom selection in the correct order. If no CB and/or CG is found
            then this method returns ``None``.

        .. versionchanged:: 1.0.0
            Added arguments for flexible atom names and refactored code for
            faster atom matching with boolean arrays.

        .. versionadded:: 0.7.5
        """
        names = [n_name, ca_name, cb_name, cg_name]
        ags = [residue.atoms[residue.atoms.names == n] for n in names]
        if any(len(ag)!= 1 for ag in ags):
            return None
        return sum(ags)

    transplants[Residue].append(('chi1_selection', chi1_selection))

    def chi1_selections(residues, n_name='N', ca_name='CA', cb_name='CB', 
                        cg_name='CG'):
        """Select list of AtomGroups corresponding to the chi1 sidechain dihedral 
        N-CA-CB-CG.

        Parameters
        ----------
        c_name: str (optional)
            name for the backbone C atom
        ca_name: str (optional)
            name for the alpha-carbon atom
        cb_name: str (optional)
            name for the beta-carbon atom
        cg_name: str (optional)
            name for the gamma-carbon atom

        Returns
        -------
        List of AtomGroups
            4-atom selections in the correct order. If no CB and/or CG is found
            then the corresponding item in the list is ``None``.

        .. versionadded:: 1.0.0
        """
        results = np.array([None]*len(residues))
        names = [n_name, ca_name, cb_name, cg_name]
        keep = [all(sum(r.atoms.names == n) == 1 for n in names) 
                for r in residues]
        keepix = np.where(keep)[0]
        residues = residues[keep]

        atnames = residues.atoms.names
        ags = [residues.atoms[atnames == n] for n in names]
        results[keepix] = [sum(atoms) for atoms in zip(*ags)]
        return list(results)
    
    transplants[ResidueGroup].append(('chi1_selections', chi1_selections))


# TODO: update docs to property doc
class Atomtypes(AtomAttr):
    """Type for each atom"""
    attrname = 'types'
    singular = 'type'
    per_object = 'atom'
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(na)], dtype=object)


# TODO: update docs to property doc
class Elements(AtomAttr):
    """Element for each atom"""
    attrname = 'elements'
    singular = 'element'
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(na)], dtype=object)


# TODO: update docs to property doc
class Radii(AtomAttr):
    """Radii for each atom"""
    attrname = 'radii'
    singular = 'radius'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


class RecordTypes(AtomAttr):
    """For PDB-like formats, indicates if ATOM or HETATM

    Defaults to 'ATOM'

    .. versionchanged:: 0.20.0
       Now stores array of dtype object rather than boolean mapping
    """
    attrname = 'record_types'
    singular = 'record_type'
    per_object = 'atom'
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['ATOM'] * na, dtype=object)


class ChainIDs(AtomAttr):
    """ChainID per atom

    Note
    ----
    This is an attribute of the Atom, not Residue or Segment
    """
    attrname = 'chainIDs'
    singular = 'chainID'
    per_object = 'atom'
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(na)], dtype=object)


class Tempfactors(AtomAttr):
    """Tempfactor for atoms"""
    attrname = 'tempfactors'
    singular = 'tempfactor'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


class Masses(AtomAttr):
    attrname = 'masses'
    singular = 'mass'
    per_object = 'atom'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup,
                      Atom, Residue, Segment]
    transplants = defaultdict(list)
    dtype = np.float64

    groupdoc = """Mass of each component in the Group.

    If the Group is an AtomGroup, then the masses are for each atom. If the
    Group is a ResidueGroup or SegmentGroup, the masses are for each residue or
    segment, respectively. These are obtained by summation of the member atoms
    for each component.
    """

    singledoc = """Mass of the component."""

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)

    def get_residues(self, rg):
        resatoms = self.top.tt.residues2atoms_2d(rg.ix)

        if isinstance(rg._ix, numbers.Integral):
            # for a single residue
            masses = self.values[resatoms].sum()
        else:
            # for a residuegroup
            masses = np.empty(len(rg))
            for i, row in enumerate(resatoms):
                masses[i] = self.values[row].sum()

        return masses

    def get_segments(self, sg):
        segatoms = self.top.tt.segments2atoms_2d(sg.ix)

        if isinstance(sg._ix, numbers.Integral):
            # for a single segment
            masses = self.values[tuple(segatoms)].sum()
        else:
            # for a segmentgroup
            masses = np.array([self.values[row].sum() for row in segatoms])

        return masses

    @warn_if_not_unique
    @check_pbc_and_unwrap
    def center_of_mass(group, pbc=False, compound='group', unwrap=False):
        """Center of mass of (compounds of) the group.

        Computes the center of mass of :class:`Atoms<Atom>` in the group.
        Centers of mass per :class:`Residue`, :class:`Segment`, molecule, or
        fragment can be obtained by setting the `compound` parameter
        accordingly. If the masses of a compound sum up to zero, the
        center of mass coordinates of that compound will be ``nan`` (not a
        number).

        Parameters
        ----------
        pbc : bool, optional
            If ``True`` and `compound` is ``'group'``, move all atoms to the
            primary unit cell before calculation.
            If ``True`` and `compound` is ``'segments'`` or ``'residues'``, the
            centers of mass of each compound will be calculated without moving
            any atoms to keep the compounds intact. Instead, the resulting
            center-of-mass position vectors will be moved to the primary unit
            cell after calculation.
        compound : {'group', 'segments', 'residues', 'molecules', 'fragments'},\
                   optional
            If ``'group'``, the center of mass of all atoms in the group will
            be returned as a single position vector. Otherwise, the centers of
            mass of each :class:`Segment`, :class:`Residue`, molecule, or
            fragment will be returned as an array of position vectors, i.e. a 2d
            array.
            Note that, in any case, *only* the positions of :class:`Atoms<Atom>`
            *belonging to the group* will be taken into account.
        unwrap : bool, optional
            If ``True``, compounds will be unwrapped before computing their centers.

        Returns
        -------
        center : numpy.ndarray
            Position vector(s) of the center(s) of mass of the group.
            If `compound` was set to ``'group'``, the output will be a single
            position vector.
            If `compound` was set to ``'segments'`` or ``'residues'``, the
            output will be a 2d coordinate array of shape ``(n, 3)`` where ``n``
            is the number of compounds.

        Note
        ----
        * This method can only be accessed if the underlying topology has
          information about atomic masses.


        .. versionchanged:: 0.8 Added `pbc` parameter
        .. versionchanged:: 0.19.0 Added `compound` parameter
        .. versionchanged:: 0.20.0 Added ``'molecules'`` and ``'fragments'``
            compounds
        .. versionchanged:: 0.20.0 Added `unwrap` parameter
        """
        atoms = group.atoms
        return atoms.center(weights=atoms.masses, pbc=pbc, compound=compound, unwrap=unwrap)

    transplants[GroupBase].append(
        ('center_of_mass', center_of_mass))

    @warn_if_not_unique
    def total_mass(group, compound='group'):
        """Total mass of (compounds of) the group.

        Computes the total mass of :class:`Atoms<Atom>` in the group.
        Total masses per :class:`Residue`, :class:`Segment`, molecule, or
        fragment can be obtained by setting the `compound` parameter
        accordingly.

        Parameters
        ----------
        compound : {'group', 'segments', 'residues', 'molecules', 'fragments'},\
                   optional
            If ``'group'``, the total mass of all atoms in the group will be
            returned as a single value. Otherwise, the total masses per
            :class:`Segment`, :class:`Residue`, molecule, or fragment will be
            returned as a 1d array.
            Note that, in any case, *only* the masses of :class:`Atoms<Atom>`
            *belonging to the group* will be taken into account.

        Returns
        -------
        float or numpy.ndarray
            Total mass of (compounds of) the group.
            If `compound` was set to ``'group'``, the output will be a single
            value. Otherwise, the output will be a 1d array of shape ``(n,)``
            where ``n`` is the number of compounds.


        .. versionchanged:: 0.20.0 Added `compound` parameter
        """
        return group.accumulate("masses", compound=compound)

    transplants[GroupBase].append(
        ('total_mass', total_mass))

    @warn_if_not_unique
    @check_pbc_and_unwrap
    def moment_of_inertia(group, pbc=False, **kwargs):
        """Tensor moment of inertia relative to center of mass as 3x3 numpy
        array.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]


        .. versionchanged:: 0.8 Added *pbc* keyword
        .. versionchanged:: 0.20.0 Added `unwrap` parameter

        """
        atomgroup = group.atoms
        unwrap = kwargs.pop('unwrap', False)
        compound = kwargs.pop('compound', 'group')

        com = atomgroup.center_of_mass(
            pbc=pbc, unwrap=unwrap, compound=compound)
        if compound != 'group':
            com = (com * group.masses[:, None]
                   ).sum(axis=0) / group.masses.sum()

        if pbc:
            pos = atomgroup.pack_into_box(inplace=False) - com
        elif unwrap:
            pos = atomgroup.unwrap(compound=compound, inplace=False) - com
        else:
            pos = atomgroup.positions - com

        masses = atomgroup.masses
        # Create the inertia tensor
        # m_i = mass of atom i
        # (x_i, y_i, z_i) = pos of atom i
        # Ixx = sum(m_i*(y_i^2+z_i^2));
        # Iyy = sum(m_i*(x_i^2+z_i^2));
        # Izz = sum(m_i*(x_i^2+y_i^2))
        # Ixy = Iyx = -1*sum(m_i*x_i*y_i)
        # Ixz = Izx = -1*sum(m_i*x_i*z_i)
        # Iyz = Izy = -1*sum(m_i*y_i*z_i)
        tens = np.zeros((3, 3), dtype=np.float64)
        # xx
        tens[0][0] = (masses * (pos[:, 1] ** 2 + pos[:, 2] ** 2)).sum()
        # xy & yx
        tens[0][1] = tens[1][0] = - (masses * pos[:, 0] * pos[:, 1]).sum()
        # xz & zx
        tens[0][2] = tens[2][0] = - (masses * pos[:, 0] * pos[:, 2]).sum()
        # yy
        tens[1][1] = (masses * (pos[:, 0] ** 2 + pos[:, 2] ** 2)).sum()
        # yz + zy
        tens[1][2] = tens[2][1] = - (masses * pos[:, 1] * pos[:, 2]).sum()
        # zz
        tens[2][2] = (masses * (pos[:, 0] ** 2 + pos[:, 1] ** 2)).sum()

        return tens

    transplants[GroupBase].append(
        ('moment_of_inertia', moment_of_inertia))

    @warn_if_not_unique
    def radius_of_gyration(group, pbc=False, **kwargs):
        """Radius of gyration.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]


        .. versionchanged:: 0.8 Added *pbc* keyword

        """
        atomgroup = group.atoms
        masses = atomgroup.masses

        com = atomgroup.center_of_mass(pbc=pbc)
        if pbc:
            recenteredpos = atomgroup.pack_into_box(inplace=False) - com
        else:
            recenteredpos = atomgroup.positions - com

        rog_sq = np.sum(masses * np.sum(recenteredpos**2,
                                        axis=1)) / atomgroup.total_mass()

        return np.sqrt(rog_sq)

    transplants[GroupBase].append(
        ('radius_of_gyration', radius_of_gyration))

    @warn_if_not_unique
    def shape_parameter(group, pbc=False, **kwargs):
        """Shape parameter.

        See [Dima2004a]_ for background information.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]


        References
        ----------
        .. [Dima2004a] Dima, R. I., & Thirumalai, D. (2004). Asymmetry
           in the shapes of folded and denatured states of
           proteins. *J Phys Chem B*, 108(21),
           6564-6570. doi:`10.1021/jp037128y
           <https://doi.org/10.1021/jp037128y>`_


        .. versionadded:: 0.7.7
        .. versionchanged:: 0.8 Added *pbc* keyword

        """
        atomgroup = group.atoms
        masses = atomgroup.masses

        com = atomgroup.center_of_mass(pbc=pbc)
        if pbc:
            recenteredpos = atomgroup.pack_into_box(inplace=False) - com
        else:
            recenteredpos = atomgroup.positions - com
        tensor = np.zeros((3, 3))

        for x in range(recenteredpos.shape[0]):
            tensor += masses[x] * np.outer(recenteredpos[x, :],
                                           recenteredpos[x, :])
        tensor /= atomgroup.total_mass()
        eig_vals = np.linalg.eigvalsh(tensor)
        shape = 27.0 * np.prod(eig_vals - np.mean(eig_vals)
                               ) / np.power(np.sum(eig_vals), 3)

        return shape

    transplants[GroupBase].append(
        ('shape_parameter', shape_parameter))

    @warn_if_not_unique
    @check_pbc_and_unwrap
    def asphericity(group, pbc=False, unwrap=None, compound='group'):
        """Asphericity.

        See [Dima2004b]_ for background information.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]
        unwrap : bool, optional
            If ``True``, compounds will be unwrapped before computing their centers.
        compound : {'group', 'segments', 'residues', 'molecules', 'fragments'}, optional
            Which type of component to keep together during unwrapping.


        References
        ----------

        .. [Dima2004b] Dima, R. I., & Thirumalai, D. (2004). Asymmetry
           in the shapes of folded and denatured states of
           proteins. *J Phys Chem B*, 108(21),
           6564-6570. doi:`10.1021/jp037128y
           <https://doi.org/10.1021/jp037128y>`_


        .. versionadded:: 0.7.7
        .. versionchanged:: 0.8 Added *pbc* keyword
        .. versionchanged:: 0.20.0 Added *unwrap* and *compound* parameter

        """
        atomgroup = group.atoms
        masses = atomgroup.masses

        com = atomgroup.center_of_mass(
            pbc=pbc, unwrap=unwrap, compound=compound)
        if compound != 'group':
            com = (com * group.masses[:, None]
                   ).sum(axis=0) / group.masses.sum()

        if pbc:
            recenteredpos = (atomgroup.pack_into_box(inplace=False) - com)
        elif unwrap:
            recenteredpos = (atomgroup.unwrap(inplace=False) - com)
        else:
            recenteredpos = (atomgroup.positions - com)

        tensor = np.zeros((3, 3))
        for x in range(recenteredpos.shape[0]):
            tensor += masses[x] * np.outer(recenteredpos[x],
                                           recenteredpos[x])

        tensor /= atomgroup.total_mass()
        eig_vals = np.linalg.eigvalsh(tensor)
        shape = (3.0 / 2.0) * (np.sum((eig_vals - np.mean(eig_vals))**2) /
                               np.sum(eig_vals)**2)

        return shape

    transplants[GroupBase].append(
        ('asphericity', asphericity))

    @warn_if_not_unique
    def principal_axes(group, pbc=False):
        """Calculate the principal axes from the moment of inertia.

        e1,e2,e3 = AtomGroup.principal_axes()

        The eigenvectors are sorted by eigenvalue, i.e. the first one
        corresponds to the highest eigenvalue and is thus the first principal
        axes.

        The eigenvectors form a right-handed coordinate system.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]

        Returns
        -------
        axis_vectors : array
            3 x 3 array with ``v[0]`` as first, ``v[1]`` as second, and
            ``v[2]`` as third eigenvector.


        .. versionchanged:: 0.8 Added *pbc* keyword
        .. versionchanged:: 1.0.0 
            Always return principal axes in right-hand convention.

        """
        atomgroup = group.atoms
        e_val, e_vec = np.linalg.eig(atomgroup.moment_of_inertia(pbc=pbc))

        # Sort
        indices = np.argsort(e_val)[::-1]
        # Make transposed in more logical form. See Issue 33.
        e_vec = e_vec[:, indices].T

        # Make sure the right hand convention is followed
        if np.dot(np.cross(e_vec[0], e_vec[1]), e_vec[2]) < 0:
            e_vec *= -1

        return e_vec

    transplants[GroupBase].append(
        ('principal_axes', principal_axes))

    def align_principal_axis(group, axis, vector):
        """Align principal axis with index `axis` with `vector`.

        Parameters
        ----------
        axis : {0, 1, 2}
            Index of the principal axis (0, 1, or 2), as produced by
            :meth:`~principal_axes`.
        vector : array_like
            Vector to align principal axis with.

        Notes
        -----
        To align the long axis of a channel (the first principal axis, i.e.
        *axis* = 0) with the z-axis::

          u.atoms.align_principal_axis(0, [0,0,1])
          u.atoms.write("aligned.pdb")
        """
        p = group.principal_axes()[axis]
        angle = np.degrees(mdamath.angle(p, vector))
        ax = transformations.rotaxis(p, vector)
        # print "principal[%d] = %r" % (axis, p)
        # print "axis = %r, angle = %f deg" % (ax, angle)
        return group.rotateby(angle, ax)

    transplants[GroupBase].append(
        ('align_principal_axis', align_principal_axis))


# TODO: update docs to property doc
class Charges(AtomAttr):
    attrname = 'charges'
    singular = 'charge'
    per_object = 'atom'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup,
                      Atom, Residue, Segment]
    transplants = defaultdict(list)
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)

    def get_residues(self, rg):
        resatoms = self.top.tt.residues2atoms_2d(rg.ix)

        if isinstance(rg._ix, numbers.Integral):
            charges = self.values[resatoms].sum()
        else:
            charges = np.empty(len(rg))
            for i, row in enumerate(resatoms):
                charges[i] = self.values[row].sum()

        return charges

    def get_segments(self, sg):
        segatoms = self.top.tt.segments2atoms_2d(sg.ix)

        if isinstance(sg._ix, numbers.Integral):
            # for a single segment
            charges = self.values[tuple(segatoms)].sum()
        else:
            # for a segmentgroup
            charges = np.array([self.values[row].sum() for row in segatoms])

        return charges

    @warn_if_not_unique
    def total_charge(group, compound='group'):
        """Total charge of (compounds of) the group.

        Computes the total charge of :class:`Atoms<Atom>` in the group.
        Total charges per :class:`Residue`, :class:`Segment`, molecule, or
        fragment can be obtained by setting the `compound` parameter
        accordingly.

        Parameters
        ----------
        compound : {'group', 'segments', 'residues', 'molecules', 'fragments'},\
                   optional
            If 'group', the total charge of all atoms in the group will
            be returned as a single value. Otherwise, the total charges per
            :class:`Segment`, :class:`Residue`, molecule, or fragment
            will be returned as a 1d array.
            Note that, in any case, *only* the charges of :class:`Atoms<Atom>`
            *belonging to the group* will be taken into account.

        Returns
        -------
        float or numpy.ndarray
            Total charge of (compounds of) the group.
            If `compound` was set to ``'group'``, the output will be a single
            value. Otherwise, the output will be a 1d array of shape ``(n,)``
            where ``n`` is the number of compounds.


        .. versionchanged:: 0.20.0 Added `compound` parameter
        """
        return group.accumulate("charges", compound=compound)

    transplants[GroupBase].append(
        ('total_charge', total_charge))


# TODO: update docs to property doc
class Bfactors(AtomAttr):
    """Crystallographic B-factors in A**2 for each atom"""
    attrname = 'bfactors'
    singular = 'bfactor'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


# TODO: update docs to property doc
class Occupancies(AtomAttr):
    attrname = 'occupancies'
    singular = 'occupancy'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


# TODO: update docs to property doc
class AltLocs(AtomAttr):
    """AltLocs for each atom"""
    attrname = 'altLocs'
    singular = 'altLoc'
    per_object = 'atom'
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(na)], dtype=object)


class GBScreens(AtomAttr):
    """Generalized Born screening factor"""
    attrname = 'gbscreens'
    singular = 'gbscreen'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


class SolventRadii(AtomAttr):
    """Intrinsic solvation radius"""
    attrname = 'solventradii'
    singular = 'solventradius'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


class NonbondedIndices(AtomAttr):
    """Nonbonded index (AMBER)"""
    attrname = 'nbindices'
    singular = 'nbindex'
    per_object = 'atom'
    dtype = int

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na, dtype=np.int32)


class RMins(AtomAttr):
    """The Rmin/2 LJ parameter"""
    attrname = 'rmins'
    singular = 'rmin'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


class Epsilons(AtomAttr):
    """The epsilon LJ parameter"""
    attrname = 'epsilons'
    singular = 'epsilon'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


class RMin14s(AtomAttr):
    """The Rmin/2 LJ parameter for 1-4 interactions"""
    attrname = 'rmin14s'
    singular = 'rmin14'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


class Epsilon14s(AtomAttr):
    """The epsilon LJ parameter for 1-4 interactions"""
    attrname = 'epsilon14s'
    singular = 'epsilon14'
    per_object = 'atom'
    dtype = float

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


class ResidueAttr(TopologyAttr):
    attrname = 'residueattrs'
    singular = 'residueattr'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup, Residue]
    per_object = 'residue'

    def get_atoms(self, ag):
        rix = self.top.tt.atoms2residues(ag.ix)
        return self.values[rix]

    def set_atoms(self, ag, values):
        raise _wronglevel_error(self, ag)

    def get_residues(self, rg):
        return self.values[rg.ix]

    @_check_length
    def set_residues(self, rg, values):
        self.values[rg.ix] = values

    def get_segments(self, sg):
        """By default, the values for each residue present in the set of
        segments are returned in a single array. This behavior can be overriden
        in child attributes.

        """
        rixs = self.top.tt.segments2residues_2d(sg.ix)
        return [self.values[rix] for rix in rixs]

    def set_segments(self, sg, values):
        raise _wronglevel_error(self, sg)


# TODO: update docs to property doc
class Resids(ResidueAttr):
    """Residue ID"""
    attrname = 'resids'
    singular = 'resid'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup, Atom, Residue]
    dtype = int

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.arange(1, nr + 1)


# TODO: update docs to property doc
class Resnames(ResidueAttr):
    attrname = 'resnames'
    singular = 'resname'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup, Atom, Residue]
    transplants = defaultdict(list)
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)

    def sequence(self, **kwargs):
        """Returns the amino acid sequence.

        The format of the sequence is selected with the keyword *format*:

        ============== ============================================
        *format*       description
        ============== ============================================
        'SeqRecord'    :class:`Bio.SeqRecord.SeqRecord` (default)
        'Seq'          :class:`Bio.Seq.Seq`
        'string'       string
        ============== ============================================

        The sequence is returned by default (keyword ``format = 'SeqRecord'``)
        as a :class:`Bio.SeqRecord.SeqRecord` instance, which can then be
        further processed. In this case, all keyword arguments (such as the
        *id* string or the *name* or the *description*) are directly passed to
        :class:`Bio.SeqRecord.SeqRecord`.

        If the keyword *format* is set to ``'Seq'``, all *kwargs* are ignored
        and a :class:`Bio.Seq.Seq` instance is returned. The difference to the
        record is that the record also contains metadata and can be directly
        used as an input for other functions in :mod:`Bio`.

        If the keyword *format* is set to ``'string'``, all *kwargs* are
        ignored and a Python string is returned.

        .. rubric:: Example: Write FASTA file

        Use :func:`Bio.SeqIO.write`, which takes sequence records::

           import Bio.SeqIO

           # get the sequence record of a protein component of a Universe
           protein = u.select_atoms("protein")
           record = protein.sequence(id="myseq1", name="myprotein")

           Bio.SeqIO.write(record, "single.fasta", "fasta")

        A FASTA file with multiple entries can be written with ::

           Bio.SeqIO.write([record1, record2, ...], "multi.fasta", "fasta")

        Parameters
        ----------
        format : string, optional
           - ``"string"``: return sequence as a string of 1-letter codes
           - ``"Seq"``: return a :class:`Bio.Seq.Seq` instance
           - ``"SeqRecord"``: return a :class:`Bio.SeqRecord.SeqRecord`
             instance

            Default is ``"SeqRecord"``
        id : optional
           Sequence ID for SeqRecord (should be different for different
           sequences)
        name : optional
           Name of the protein.
        description : optional
           Short description of the sequence.
        kwargs : optional
           Any other keyword arguments that are understood by
           class:`Bio.SeqRecord.SeqRecord`.

        Raises
        ------
        :exc:`ValueError` if a residue name cannot be converted to a
        1-letter IUPAC protein amino acid code; make sure to only
        select protein residues.

        :exc:`TypeError` if an unknown *format* is selected.


        .. versionadded:: 0.9.0
        """
        formats = ('string', 'Seq', 'SeqRecord')

        format = kwargs.pop("format", "SeqRecord")
        if format not in formats:
            raise TypeError("Unknown format='{0}': must be one of: {1}".format(
                format, ", ".join(formats)))
        try:
            sequence = "".join([convert_aa_code(r)
                                for r in self.residues.resnames])
        except KeyError as err:
            six.raise_from(ValueError("AtomGroup contains a residue name '{0}' that "
                                      "does not have a IUPAC protein 1-letter "
                                      "character".format(err.message)), None)
        if format == "string":
            return sequence
        seq = Bio.Seq.Seq(sequence)
        if format == "Seq":
            return seq
        return Bio.SeqRecord.SeqRecord(seq, **kwargs)

    transplants[ResidueGroup].append(
        ('sequence', sequence))


# TODO: update docs to property doc
class Resnums(ResidueAttr):
    attrname = 'resnums'
    singular = 'resnum'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup, Atom, Residue]
    dtype = int

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.arange(1, nr + 1)


class ICodes(ResidueAttr):
    """Insertion code for Atoms"""
    attrname = 'icodes'
    singular = 'icode'
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)


class Moltypes(ResidueAttr):
    """Name of the molecule type

    Two molecules that share a molecule type share a common template topology.
    """
    attrname = 'moltypes'
    singular = 'moltype'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup, Atom, Residue]
    dtype = object


class Molnums(ResidueAttr):
    """Name of the molecule type

    Two molecules that share a molecule type share a common template topology.
    """
    attrname = 'molnums'
    singular = 'molnum'
    target_classes = [AtomGroup, ResidueGroup, Atom, Residue]
    dtype = np.intp

# segment attributes


class SegmentAttr(TopologyAttr):
    """Base class for segment attributes.

    """
    attrname = 'segmentattrs'
    singular = 'segmentattr'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup, Segment]
    per_object = 'segment'

    def get_atoms(self, ag):
        six = self.top.tt.atoms2segments(ag.ix)
        return self.values[six]

    def set_atoms(self, ag, values):
        raise _wronglevel_error(self, ag)

    def get_residues(self, rg):
        six = self.top.tt.residues2segments(rg.ix)
        return self.values[six]

    def set_residues(self, rg, values):
        raise _wronglevel_error(self, rg)

    def get_segments(self, sg):
        return self.values[sg.ix]

    @_check_length
    def set_segments(self, sg, values):
        self.values[sg.ix] = values


# TODO: update docs to property doc
class Segids(SegmentAttr):
    attrname = 'segids'
    singular = 'segid'
    target_classes = [AtomGroup, ResidueGroup, SegmentGroup,
                      Atom, Residue, Segment]
    transplants = defaultdict(list)
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(ns)], dtype=object)


def _check_connection_values(func):
    """
    Checks values passed to _Connection methods for:
     - appropriate number of atom indices
     - coerces them to tuples of ints (for hashing)
     - ensures that first value is less than last (reversibility & hashing)

    .. versionadded:: 1.0.0

    """
    @functools.wraps(func)
    def wrapper(self, values, *args, **kwargs):
        if not all(len(x) == self._n_atoms
                   and all(isinstance(y, (int, np.integer)) for y in x)
                   for x in values):
            raise ValueError(("{} must be an iterable of tuples with {}"
                              " atom indices").format(self.attrname,
                                                      self._n_atoms))
        clean = []
        for v in values:
            if v[0] > v[-1]:
                v = v[::-1]
            clean.append(tuple(v))

        return func(self, clean, *args, **kwargs)
    return wrapper


class _Connection(AtomAttr):
    """Base class for connectivity between atoms

    .. versionchanged:: 1.0.0
        Added type checking to atom index values.
    """

    @_check_connection_values
    def __init__(self, values, types=None, guessed=False, order=None):
        self.values = values
        if types is None:
            types = [None] * len(values)
        self.types = types
        if guessed in (True, False):
            # if single value passed, multiply this across
            # all bonds
            guessed = [guessed] * len(values)
        self._guessed = guessed
        if order is None:
            order = [None] * len(values)
        self.order = order
        self._cache = dict()

    def copy(self):
        """Return a deepcopy of this attribute"""
        return self.__class__(copy.copy(self.values),
                              copy.copy(self.types),
                              copy.copy(self._guessed),
                              copy.copy(self.order))

    def __len__(self):
        return len(self._bondDict)

    @property
    @cached('bd')
    def _bondDict(self):
        """Lazily built mapping of atoms:bonds"""
        bd = defaultdict(list)

        for b, t, g, o in zip(self.values, self.types,
                              self._guessed, self.order):
            for a in b:
                bd[a].append((b, t, g, o))
        return bd

    def set_atoms(self, ag):
        return NotImplementedError("Cannot set bond information")

    def get_atoms(self, ag):
        try:
            unique_bonds = set(itertools.chain(
                *[self._bondDict[a] for a in ag.ix]))
        except TypeError:
            # maybe we got passed an Atom
            unique_bonds = self._bondDict[ag.ix]
        bond_idx, types, guessed, order = np.hsplit(
            np.array(sorted(unique_bonds)), 4)
        bond_idx = np.array(bond_idx.ravel().tolist(), dtype=np.int32)
        types = types.ravel()
        guessed = guessed.ravel()
        order = order.ravel()
        return TopologyGroup(bond_idx, ag.universe,
                             self.singular[:-1],
                             types,
                             guessed,
                             order)

    @_check_connection_values
    def _add_bonds(self, values, types=None, guessed=True, order=None):
        if types is None:
            types = itertools.cycle((None,))
        if guessed in (True, False):
            guessed = itertools.cycle((guessed,))
        if order is None:
            order = itertools.cycle((None,))

        existing = set(self.values)
        for v, t, g, o in zip(values, types, guessed, order):
            if v not in existing:
                self.values.append(v)
                self.types.append(t)
                self._guessed.append(g)
                self.order.append(o)
        # kill the old cache of bond Dict
        try:
            del self._cache['bd']
        except KeyError:
            pass

    @_check_connection_values
    def _delete_bonds(self, values):
        """
        .. versionadded:: 1.0.0
        """

        to_check = set(values)
        self_values = set(self.values)
        if not to_check.issubset(self_values):
            missing = to_check-self_values
            indices = ', '.join(map(str, missing))
            raise ValueError(('Cannot delete nonexistent '
                              '{attrname} with atom indices:'
                              '{indices}').format(attrname=self.attrname,
                                                  indices=indices))
        idx = [self.values.index(v) for v in to_check]
        for i in sorted(idx, reverse=True):
            del self.values[i]

        for attr in ('types', '_guessed', 'order'):
            arr = np.array(getattr(self, attr), dtype='object')
            new = np.delete(arr, idx)
            setattr(self, attr, list(new))
        # kill the old cache of bond Dict
        try:
            del self._cache['bd']
        except KeyError:
            pass


class Bonds(_Connection):
    """Bonds between two atoms

    Must be initialised by a list of zero based tuples.
    These indices refer to the atom indices.
    E.g., ` [(0, 1), (1, 2), (2, 3)]`

    Also adds the `bonded_atoms`, `fragment` and `fragments`
    attributes.
    """
    attrname = 'bonds'
    # Singular is the same because one Atom might have
    # many bonds, so still asks for "bonds" in the plural
    singular = 'bonds'
    transplants = defaultdict(list)
    _n_atoms = 2

    def bonded_atoms(self):
        """An :class:`~MDAnalysis.core.groups.AtomGroup` of all
        :class:`Atoms<MDAnalysis.core.groups.Atom>` bonded to this
        :class:`~MDAnalysis.core.groups.Atom`."""
        idx = [b.partner(self).index for b in self.bonds]
        return self.universe.atoms[idx]

    transplants[Atom].append(
        ('bonded_atoms', property(bonded_atoms, None, None,
                                  bonded_atoms.__doc__)))

    def fragindex(self):
        """The index (ID) of the
        :class:`~MDAnalysis.core.topologyattrs.Bonds.fragment` this
        :class:`~MDAnalysis.core.groups.Atom` is part of.

        Note
        ----
        This property is only accessible if the underlying topology contains
        bond information.


        .. versionadded:: 0.20.0
        """
        return self.universe._fragdict[self.ix].ix

    def fragindices(self):
        """The
        :class:`fragment indices<MDAnalysis.core.topologyattrs.Bonds.fragindex>`
        of all :class:`Atoms<MDAnalysis.core.groups.Atom>` in this
        :class:`~MDAnalysis.core.groups.AtomGroup`.

        A :class:`numpy.ndarray` with
        :attr:`~numpy.ndarray.shape`\ ``=(``\ :attr:`~AtomGroup.n_atoms`\ ``,)``
        and :attr:`~numpy.ndarray.dtype`\ ``=numpy.int64``.

        Note
        ----
        This property is only accessible if the underlying topology contains
        bond information.


        .. versionadded:: 0.20.0
        """
        fragdict = self.universe._fragdict
        return np.array([fragdict[aix].ix for aix in self.ix], dtype=np.intp)

    def fragment(self):
        """An :class:`~MDAnalysis.core.groups.AtomGroup` representing the
        fragment this :class:`~MDAnalysis.core.groups.Atom` is part of.

        A fragment is a
        :class:`group of atoms<MDAnalysis.core.groups.AtomGroup>` which are
        interconnected by :class:`~MDAnalysis.core.topologyattrs.Bonds`, i.e.,
        there exists a path along one
        or more :class:`~MDAnalysis.core.topologyattrs.Bonds` between any pair
        of :class:`Atoms<MDAnalysis.core.groups.Atom>`
        within a fragment. Thus, a fragment typically corresponds to a molecule.

        Note
        ----
        This property is only accessible if the underlying topology contains
        bond information.


        .. versionadded:: 0.9.0
        """
        return self.universe._fragdict[self.ix].fragment

    def fragments(self):
        """Read-only :class:`tuple` of
        :class:`fragments<MDAnalysis.core.topologyattrs.Bonds.fragment>`.

        Contains all fragments that
        any :class:`~MDAnalysis.core.groups.Atom` in this
        :class:`~MDAnalysis.core.groups.AtomGroup` is part of.

        A fragment is a
        :class:`group of atoms<MDAnalysis.core.groups.AtomGroup>` which are
        interconnected by :class:`~MDAnalysis.core.topologyattrs.Bonds`, i.e.,
        there exists a path along one
        or more :class:`~MDAnalysis.core.topologyattrs.Bonds` between any pair
        of :class:`Atoms<MDAnalysis.core.groups.Atom>`
        within a fragment. Thus, a fragment typically corresponds to a molecule.

        Note
        ----
        * This property is only accessible if the underlying topology contains
          bond information.
        * The contents of the fragments may extend beyond the contents of this
          :class:`~MDAnalysis.core.groups.AtomGroup`.


        .. versionadded:: 0.9.0
        """
        fragdict = self.universe._fragdict
        return tuple(sorted(set(fragdict[aix].fragment for aix in self.ix),
                            key=lambda x: x[0].ix))

    def n_fragments(self):
        """The number of unique
        :class:`~MDAnalysis.core.topologyattrs.Bonds.fragments` the
        :class:`Atoms<MDAnalysis.core.groups.Atom>` of this
        :class:`~MDAnalysis.core.groups.AtomGroup` are part of.

        Note
        ----
        This property is only accessible if the underlying topology contains
        bond information.


        .. versionadded:: 0.20.0
        """
        return len(unique_int_1d(self.fragindices))

    transplants[Atom].append(
        ('fragment', property(fragment, None, None,
                              fragment.__doc__)))

    transplants[Atom].append(
        ('fragindex', property(fragindex, None, None,
                               fragindex.__doc__)))

    transplants[AtomGroup].append(
        ('fragments', property(fragments, None, None,
                               fragments.__doc__)))

    transplants[AtomGroup].append(
        ('fragindices', property(fragindices, None, None,
                                 fragindices.__doc__)))

    transplants[AtomGroup].append(
        ('n_fragments', property(n_fragments, None, None,
                                 n_fragments.__doc__)))


class UreyBradleys(_Connection):
    """Angles between two atoms

    Initialise with a list of 2 long tuples

    These indices refer to the atom indices.

    .. versionadded:: 1.0.0
    """
    attrname = 'ureybradleys'
    singular = 'ureybradleys'
    transplants = defaultdict(list)
    _n_atoms = 2


class Angles(_Connection):
    """Angles between three atoms

    Initialise with a list of 3 long tuples
    E.g.,  `[(0, 1, 2), (1, 2, 3), (2, 3, 4)]`

    These indices refer to the atom indices.
    """
    attrname = 'angles'
    singular = 'angles'
    transplants = defaultdict(list)
    _n_atoms = 3


class Dihedrals(_Connection):
    """A connection between four sequential atoms"""
    attrname = 'dihedrals'
    singular = 'dihedrals'
    transplants = defaultdict(list)
    _n_atoms = 4


class Impropers(_Connection):
    """An imaginary dihedral between four atoms"""
    attrname = 'impropers'
    singular = 'impropers'
    transplants = defaultdict(list)
    _n_atoms = 4


class CMaps(_Connection):
    """
    A connection between five atoms
    .. versionadded:: 1.0.0
    """
    attrname = 'cmaps'
    singular = 'cmaps'
    transplants = defaultdict(list)
    _n_atoms = 5
