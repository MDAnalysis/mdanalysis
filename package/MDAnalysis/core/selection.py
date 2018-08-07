# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""Atom selection Hierarchy --- :mod:`MDAnalysis.core.selection`
=============================================================

This module contains objects that represent selections. They are
constructed and then applied to the group.

In general, :meth:`Parser.parse` creates a :class:`Selection` object
from a selection string. This :class:`Selection` object is then passed
an :class:`~MDAnalysis.core.groups.AtomGroup` through its
:meth:`~MDAnalysis.core.groups.AtomGroup.apply` method to apply the
``Selection`` to the ``AtomGroup``.

This is all invisible to the user through the
:meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` method of an
:class:`~MDAnalysis.core.groups.AtomGroup`.

"""
from __future__ import division, absolute_import
import six
from six.moves import zip

import collections
import re
import functools
import warnings

import numpy as np
from numpy.lib.utils import deprecate

from MDAnalysis.lib.pkdtree import PeriodicKDTree
from MDAnalysis.lib.util import unique_int_1d
from MDAnalysis.core import flags
from ..lib import distances
from ..exceptions import SelectionError, NoDataError


def is_keyword(val):
    """Is val a selection keyword?

    Returns False on any of the following strings:
      - keys in SELECTIONDICT (tokens from Selection objects)
      - keys in OPERATIONS (tokens from LogicOperations)
      - (Parentheses)
      - The value `None` (used as EOF in selection strings)
    """
    return (val in _SELECTIONDICT or
            val in _OPERATIONS or
            val in ['(', ')'] or
            val is None)


def grab_not_keywords(tokens):
    """Pop tokens from the left until you hit a keyword

    Parameters
    ----------
    tokens : collections.deque
        deque of strings, some tokens some not

    Returns
    -------
    values : list of strings
        All non keywords found until a keyword was hit

    Note
    ----
    This function pops the values from the deque

    Examples
    --------
    grab_not_keywords(['H', 'and','resname', 'MET'])
    >>> ['H']

    grab_not_keywords(['H', 'Ca', 'N', 'and','resname', 'MET'])
    >>> ['H', 'Ca' ,'N']

    grab_not_keywords(['and','resname', 'MET'])
    >>> []
    """
    values = []
    while not is_keyword(tokens[0]):
        val = tokens.popleft()
        # Insert escape characters here to use keywords as names?
        values.append(val)
    return values


_SELECTIONDICT = {}
_OPERATIONS = {}
# These are named args to select_atoms that have a special meaning and must
# not be allowed as names for the 'group' keyword.
_RESERVED_KWARGS=('updating',)


# And and Or are exception and aren't strictly a Selection
# as they work on other Selections rather than doing work themselves.
# So their init is a little strange too....
class _Operationmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            _OPERATIONS[classdict['token']] = cls
        except KeyError:
            pass


class LogicOperation(six.with_metaclass(_Operationmeta, object)):
    def __init__(self, lsel, rsel):
        self.rsel = rsel
        self.lsel = lsel


class AndOperation(LogicOperation):
    token = 'and'
    precedence = 3

    def apply(self, group):
        rsel = self.rsel.apply(group)
        lsel = self.lsel.apply(group)

        # Mask which lsel indices appear in rsel
        mask = np.in1d(rsel.indices, lsel.indices)
        # and mask rsel according to that
        return rsel[mask].unique


class OrOperation(LogicOperation):
    token = 'or'
    precedence = 3

    def apply(self, group):
        lsel = self.lsel.apply(group)
        rsel = self.rsel.apply(group)

        # Find unique indices from both these AtomGroups
        # and slice master list using them
        idx = np.union1d(lsel.indices, rsel.indices).astype(np.int32)

        return group.universe.atoms[idx]


class _Selectionmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            _SELECTIONDICT[classdict['token']] = cls
        except KeyError:
            pass


class Selection(six.with_metaclass(_Selectionmeta, object)):
    pass


class AllSelection(Selection):
    token = 'all'

    def __init__(self, parser, tokens):
        pass

    def apply(self, group):
        # Check whether group is identical to the one stored
        # in the corresponding universe, in which case this
        # is returned directly. This works since the Universe.atoms
        # are unique by construction.
        if group is group.universe.atoms:
            return group
        return group[:].unique


class UnarySelection(Selection):
    def __init__(self, parser, tokens):
        sel = parser.parse_expression(self.precedence)
        self.sel = sel


class NotSelection(UnarySelection):
    token = 'not'
    precedence = 5

    def apply(self, group):
        notsel = self.sel.apply(group)
        return group[~np.in1d(group.indices, notsel.indices)].unique


class GlobalSelection(UnarySelection):
    token = 'global'
    precedence = 5

    def apply(self, group):
        return self.sel.apply(group.universe.atoms).unique


class ByResSelection(UnarySelection):
    token = 'byres'
    precedence = 1

    def apply(self, group):
        res = self.sel.apply(group)
        unique_res = unique_int_1d(res.resids.astype(np.int64))
        mask = np.in1d(group.resids, unique_res)

        return group[mask].unique


class DistanceSelection(Selection):
    """Base class for distance search based selections

    Grabs the flags for this selection
     - 'use_KDTree_routines'
     - 'use_periodic_selections'

    Populates the `apply` method with either
     - _apply_KDTree
     - _apply_distmat
    """
    def __init__(self):
        if flags['use_KDTree_routines'] in (True, 'fast', 'always'):
            self.apply = self._apply_KDTree
        else:
            self.apply = self._apply_distmat

        self.periodic = flags['use_periodic_selections']

    def validate_dimensions(self, dimensions):
        r"""Check if the system is periodic in all three-dimensions.

        Parameters
        ----------
        dimensions : numpy.ndarray
            6-item array denoting system size and angles

        Returns
        -------
        None or numpy.ndarray
            Returns argument dimensions if system is periodic in all
            three-dimensions, otherwise returns None
        """
        if self.periodic and all(dimensions[:3]):
            return dimensions
        return None

class AroundSelection(DistanceSelection):
    token = 'around'
    precedence = 1

    def __init__(self, parser, tokens):
        super(AroundSelection, self).__init__()
        self.cutoff = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)

    def _apply_KDTree(self, group):
        """KDTree based selection is about 7x faster than distmat
        for typical problems.
        """
        sel = self.sel.apply(group)
        # All atoms in group that aren't in sel
        sys = group[~np.in1d(group.indices, sel.indices)]

        if not sys:
            return sys[[]]

        box = self.validate_dimensions(group.dimensions)

        cut = self.cutoff if box is not None else None
        kdtree = PeriodicKDTree(box=box, leafsize=10)
        kdtree.set_coords(sys.positions, cutoff=cut)
        kdtree.search(sel.positions, self.cutoff)
        unique_idx = np.asarray(kdtree.get_indices())

        return sys[unique_idx.astype(np.int32)].unique

    def _apply_distmat(self, group):
        sel = self.sel.apply(group)
        sys = group[~np.in1d(group.indices, sel.indices)]

        box = self.validate_dimensions(group.dimensions)
        dist = distances.distance_array(
            sys.positions, sel.positions, box)

        mask = (dist <= self.cutoff).any(axis=1)

        return sys[mask].unique


class SphericalLayerSelection(DistanceSelection):
    token = 'sphlayer'
    precedence = 1

    def __init__(self, parser, tokens):
        super(SphericalLayerSelection, self).__init__()
        self.inRadius = float(tokens.popleft())
        self.exRadius = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)

    def _apply_KDTree(self, group):
        """Selection using KDTree and PeriodicKDTree for aperiodic and
        fully-periodic systems, respectively.
        """
        sel = self.sel.apply(group)
        box = self.validate_dimensions(group.dimensions)
        periodic = box is not None
        ref = sel.center_of_geometry(pbc=periodic)
        kdtree = PeriodicKDTree(box=box)

        cutoff = self.exRadius if box is not None else None
        kdtree.set_coords(group.positions, cutoff=cutoff)
        kdtree.search(ref, self.exRadius)
        found_ExtIndices = kdtree.get_indices()
        kdtree.search(ref, self.inRadius)
        found_IntIndices = kdtree.get_indices()
        found_indices = list(set(found_ExtIndices) - set(found_IntIndices))
        return group[found_indices].unique

    def _apply_distmat(self, group):
        sel = self.sel.apply(group)
        box = self.validate_dimensions(group.dimensions)
        periodic = box is not None
        ref = sel.center_of_geometry(pbc=periodic).reshape(1, 3).astype(
            np.float32)
        d = distances.distance_array(ref,
                                     group.positions,
                                     box=box)[0]
        mask = d < self.exRadius
        mask &= d > self.inRadius

        return group[mask].unique


class SphericalZoneSelection(DistanceSelection):
    token = 'sphzone'
    precedence = 1

    def __init__(self, parser, tokens):
        super(SphericalZoneSelection, self).__init__()
        self.cutoff = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)

    def _apply_KDTree(self, group):
        """Selection using KDTree and PeriodicKDTree for aperiodic and
        fully-periodic systems, respectively.
        """
        sel = self.sel.apply(group)
        box = self.validate_dimensions(group.dimensions)
        periodic = box is not None
        ref = sel.center_of_geometry(pbc=periodic)

        cut = self.cutoff if box is not None else None
        kdtree = PeriodicKDTree(box=box)
        kdtree.set_coords(group.positions, cutoff=cut)
        kdtree.search(ref, self.cutoff)
        found_indices = kdtree.get_indices()
        return group[found_indices].unique

    def _apply_distmat(self, group):
        sel = self.sel.apply(group)
        box = self.validate_dimensions(group.dimensions)
        periodic = box is not None
        ref = sel.center_of_geometry(pbc=periodic).reshape(1, 3).\
            astype(np.float32)
        d = distances.distance_array(ref,
                                     group.positions,
                                     box=box)[0]
        idx = d < self.cutoff
        return group[idx].unique


class CylindricalSelection(Selection):
    def __init__(self):
        self.periodic = flags['use_periodic_selections']

    def apply(self, group):
        sel = self.sel.apply(group)

        # Calculate vectors between point of interest and our group
        vecs = group.positions - sel.center_of_geometry()

        if self.periodic and not np.any(group.dimensions[:3] == 0):
            box = group.dimensions[:3]
            cyl_z_hheight = self.zmax - self.zmin

            if 2 * self.exRadius > box[0]:
                raise NotImplementedError(
                    "The diameter of the cylinder selection ({:.3f}) is larger "
                    "than the unit cell's x dimension ({:.3f}). Can only do "
                    "selections where it is smaller or equal."
                    "".format(2*self.exRadius, box[0]))
            if 2 * self.exRadius > box[1]:
                raise NotImplementedError(
                    "The diameter of the cylinder selection ({:.3f}) is larger "
                    "than the unit cell's y dimension ({:.3f}). Can only do "
                    "selections where it is smaller or equal."
                    "".format(2*self.exRadius, box[1]))
            if cyl_z_hheight > box[2]:
                raise NotImplementedError(
                    "The total length of the cylinder selection in z ({:.3f}) "
                    "is larger than the unit cell's z dimension ({:.3f}). Can "
                    "only do selections where it is smaller or equal."
                    "".format(cyl_z_hheight, box[2]))

            if np.all(group.dimensions[3:] == 90.):
                # Orthogonal version
                vecs -= box[:3] * np.rint(vecs / box[:3])
            else:
                # Triclinic version
                tribox = group.universe.trajectory.ts.triclinic_dimensions
                vecs -= tribox[2] * np.rint(vecs[:, 2] / tribox[2][2])[:, None]
                vecs -= tribox[1] * np.rint(vecs[:, 1] / tribox[1][1])[:, None]
                vecs -= tribox[0] * np.rint(vecs[:, 0] / tribox[0][0])[:, None]

        # First deal with Z dimension criteria
        mask = (vecs[:, 2] > self.zmin) & (vecs[:, 2] < self.zmax)
        # Mask out based on height to reduce number of radii comparisons
        vecs = vecs[mask]
        group = group[mask]

        # Radial vectors from sel to each in group
        radii = vecs[:, 0]**2 + vecs[:, 1]**2
        mask = radii < self.exRadius**2
        try:
            mask &= radii > self.inRadius**2
        except AttributeError:
            # Only for cylayer, cyzone doesn't have inRadius
            pass

        return group[mask].unique


class CylindricalZoneSelection(CylindricalSelection):
    token = 'cyzone'
    precedence = 1

    def __init__(self, parser, tokens):
        super(CylindricalZoneSelection, self).__init__()
        self.exRadius = float(tokens.popleft())
        self.zmax = float(tokens.popleft())
        self.zmin = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)


class CylindricalLayerSelection(CylindricalSelection):
    token = 'cylayer'
    precedence = 1

    def __init__(self, parser, tokens):
        super(CylindricalLayerSelection, self).__init__()
        self.inRadius = float(tokens.popleft())
        self.exRadius = float(tokens.popleft())
        self.zmax = float(tokens.popleft())
        self.zmin = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)


class PointSelection(DistanceSelection):
    token = 'point'

    def __init__(self, parser, tokens):
        super(PointSelection, self).__init__()
        x = float(tokens.popleft())
        y = float(tokens.popleft())
        z = float(tokens.popleft())
        self.ref = np.array([x, y, z], dtype=np.float32)
        self.cutoff = float(tokens.popleft())

    def _apply_KDTree(self, group):
        box = self.validate_dimensions(group.dimensions)
        kdtree = PeriodicKDTree(box=box)
        cut = self.cutoff if box is not None else None
        kdtree.set_coords(group.positions, cutoff=cut)
        kdtree.search(self.ref, self.cutoff)
        found_indices = kdtree.get_indices()

        return group[found_indices].unique

    def _apply_distmat(self, group):
        ref_coor = self.ref[np.newaxis, ...]

        ref_coor = np.asarray(ref_coor, dtype=np.float32)
        box = self.validate_dimensions(group.dimensions)

        dist = distances.distance_array(group.positions, ref_coor, box)
        mask = (dist <= self.cutoff).any(axis=1)
        return group[mask].unique


class AtomSelection(Selection):
    token = 'atom'

    def __init__(self, parser, tokens):
        self.segid = tokens.popleft()
        self.resid = int(tokens.popleft())
        self.name = tokens.popleft()

    def apply(self, group):
        sub = group[group.names == self.name]
        if sub:
            sub = sub[sub.resids == self.resid]
        if sub:
            sub = sub[sub.segids == self.segid]
        return sub.unique


class BondedSelection(Selection):
    token = 'bonded'
    precedence = 1

    def __init__(self, parser, tokens):
        self.sel = parser.parse_expression(self.precedence)

    def apply(self, group):
        grp = self.sel.apply(group)
        # Check if we have bonds
        if not group.bonds:
            warnings.warn("Bonded selection has 0 bonds")
            return group[[]]

        grpidx = grp.indices

        # (n, 2) array of bond indices
        bix = np.array(group.bonds.to_indices())

        idx = []
        # left side
        idx.append(bix[:, 0][np.in1d(bix[:, 1], grpidx)])
        # right side
        idx.append(bix[:, 1][np.in1d(bix[:, 0], grpidx)])

        idx = np.union1d(*idx)

        return group.universe.atoms[np.unique(idx)]


class SelgroupSelection(Selection):
    token = 'group'

    def __init__(self, parser, tokens):
        grpname = tokens.popleft()
        if grpname in _RESERVED_KWARGS:
            raise TypeError("The '{}' keyword is reserved and cannot be "
                            "used as a selection group name."
                            .format(grpname))
        try:
            self.grp = parser.selgroups[grpname]
        except KeyError:
            raise ValueError("Failed to find group: {0}".format(grpname))

    def apply(self, group):
        mask = np.in1d(group.indices, self.grp.indices)
        return group[mask]

# TODO: remove in 1.0 (should have been removed in 0.15.0)
class FullSelgroupSelection(Selection):
    token = 'fullgroup'

    def __init__(self, parser, tokens):
        grpname = tokens.popleft()
        try:
            self.grp = parser.selgroups[grpname]
        except KeyError:
            raise ValueError("Failed to find group: {0}".format(grpname))

    @deprecate(old_name='fullgroup', new_name='global group',
               message=' This will be removed in v0.15.0')
    def apply(self, group):
        return self.grp.unique


class StringSelection(Selection):
    """Selections based on text attributes

    Supports the use of wildcards at the end of strings
    """
    def __init__(self, parser, tokens):
        vals = grab_not_keywords(tokens)
        if not vals:
            raise ValueError("Unexpected token '{0}'".format(tokens[0]))

        self.values = vals

    def apply(self, group):
        mask = np.zeros(len(group), dtype=np.bool)
        for val in self.values:
            wc_pos = val.find('*')
            if wc_pos == -1:  # No wildcard found
                mask |= getattr(group, self.field) == val
            else:
                values = getattr(group, self.field).astype(np.str_)
                mask |= np.char.startswith(values, val[:wc_pos])

        return group[mask].unique


class AtomNameSelection(StringSelection):
    """Select atoms based on 'names' attribute"""
    token = 'name'
    field = 'names'


class AtomTypeSelection(StringSelection):
    """Select atoms based on 'types' attribute"""
    token = 'type'
    field = 'types'


class AtomICodeSelection(StringSelection):
    """Select atoms based on icode attribute"""
    token = 'icode'
    field = 'icodes'


class ResidueNameSelection(StringSelection):
    """Select atoms based on 'resnames' attribute"""
    token = 'resname'
    field = 'resnames'


class MoleculeTypeSelection(StringSelection):
    """Select atoms based on 'moltypes' attribute"""
    token = 'moltype'
    field = 'moltypes'


class SegmentNameSelection(StringSelection):
    """Select atoms based on 'segids' attribute"""
    token = 'segid'
    field = 'segids'


class AltlocSelection(StringSelection):
    """Select atoms based on 'altLoc' attribute"""
    token = 'altloc'
    field = 'altLocs'


class ResidSelection(Selection):
    """Select atoms based on numerical fields

    Allows the use of ':' and '-' to specify a range of values
    For example

      resid 1:10
    """
    token = 'resid'

    def __init__(self, parser, tokens):
        values = grab_not_keywords(tokens)
        if not values:
            raise ValueError("Unexpected token: '{0}'".format(tokens[0]))

        # each value in uppers and lowers is a tuple of (resid, icode)
        uppers = []
        lowers = []

        for val in values:
            m1 = re.match("(\d+)(\w?)$", val)
            if not m1 is None:
                res = m1.groups()
                lower = int(res[0]), res[1]
                upper = None, None
            else:
                # check if in appropriate format 'lower:upper' or 'lower-upper'
                # each val is one or more digits, maybe a letter
                selrange = re.match("(\d+)(\w?)[:-](\d+)(\w?)", val)
                if selrange is None:  # re.match returns None on failure
                    raise ValueError("Failed to parse value: {0}".format(val))
                res = selrange.groups()
                # resid and icode
                lower = int(res[0]), res[1]
                upper = int(res[2]), res[3]

            lowers.append(lower)
            uppers.append(upper)

        self.lowers = lowers
        self.uppers = uppers

    def apply(self, group):
        # Grab arrays here to reduce number of calls to main topology
        vals = group.resids
        try:  # optional attribute
            icodes = group.icodes
        except (AttributeError, NoDataError):
            icodes = None
            # if no icodes and icodes are part of selection, cause a fuss
            if (any(v[1] for v in self.uppers) or
                any(v[1] for v in self.lowers)):
                raise ValueError("Selection specified icodes, while the "
                                 "topology doesn't have any.")

        if not icodes is None:
            mask = self._sel_with_icodes(vals, icodes)
        else:
            mask = self._sel_without_icodes(vals)

        return group[mask].unique

    def _sel_without_icodes(self, vals):
        # Final mask that gets applied to group
        mask = np.zeros(len(vals), dtype=np.bool)

        for (u_resid, _), (l_resid, _) in zip(self.uppers, self.lowers):
            if u_resid is not None:  # range selection
                thismask = vals >= l_resid
                thismask &= vals <= u_resid
            else:  # single residue selection
                thismask = vals == l_resid

            mask |= thismask

        return mask

    def _sel_with_icodes(self, vals, icodes):
        # Final mask that gets applied to group
        mask = np.zeros(len(vals), dtype=np.bool)

        for (u_resid, u_icode), (l_resid, l_icode) in zip(self.uppers, self.lowers):
            if u_resid is not None:  # Selecting a range
                # Special case, if l_resid == u_resid, ie 163A-163C, this simplifies to:
                # all 163, and A <= icode <= C
                if l_resid == u_resid:
                    thismask = vals == l_resid
                    thismask &= icodes >= l_icode
                    thismask &= icodes <= u_icode
                # For 163A to 166B we want:
                # [START]  all 163 and icode >= 'A'
                # [MIDDLE] all of 164 and 165, any icode
                # [END]    166 and icode <= 'B'
                else:
                    # start of range
                    startmask = vals == l_resid
                    startmask &= icodes >= l_icode
                    thismask = startmask

                    # middle of range
                    mid = np.arange(l_resid + 1, u_resid)
                    if len(mid):  # if there are any resids in the middle
                        mid_beg, mid_end = mid[0], mid[-1]
                        midmask = vals >= mid_beg
                        midmask &= vals <= mid_end

                        thismask |= midmask

                    # end of range
                    endmask = vals == u_resid
                    endmask &= icodes <= u_icode

                    thismask |= endmask
            else:  # Selecting a single residue
                thismask = vals == l_resid
                thismask &= icodes == l_icode

            mask |= thismask

        return mask


class RangeSelection(Selection):
    value_offset=0

    def __init__(self, parser, tokens):
        values = grab_not_keywords(tokens)
        if not values:
            raise ValueError("Unexpected token: '{0}'".format(tokens[0]))

        uppers = []  # upper limit on any range
        lowers = []  # lower limit on any range

        for val in values:
            try:
                lower = int(val)
                upper = None
            except ValueError:
                # check if in appropriate format 'lower:upper' or 'lower-upper'
                selrange = re.match("(\d+)[:-](\d+)", val)
                if not selrange:
                    raise ValueError(
                        "Failed to parse number: {0}".format(val))
                lower, upper = np.int64(selrange.groups())

            lowers.append(lower)
            uppers.append(upper)

        self.lowers = lowers
        self.uppers = uppers

    def apply(self, group):
        mask = np.zeros(len(group), dtype=np.bool)
        vals = getattr(group, self.field) + self.value_offset

        for upper, lower in zip(self.uppers, self.lowers):
            if upper is not None:
                thismask = vals >= lower
                thismask &= vals <= upper
            else:
                thismask = vals == lower

            mask |= thismask
        return group[mask].unique


class ResnumSelection(RangeSelection):
    token = 'resnum'
    field = 'resnums'


class ByNumSelection(RangeSelection):
    token = 'bynum'
    field = 'indices'
    value_offset = 1  # queries are in 1 based indices


class MolidSelection(RangeSelection):
    token = 'molnum'
    field = 'molnums'


class ProteinSelection(Selection):
    """Consists of all residues with  recognized residue names.

    Recognized residue names in :attr:`ProteinSelection.prot_res`.

      * from the CHARMM force field::
         awk '/RESI/ {printf "'"'"%s"'"',",$2 }' top_all27_prot_lipid.rtf

      * manually added special CHARMM, OPLS/AA and Amber residue names.

    See Also
    --------
    :func:`MDAnalysis.lib.util.convert_aa_code`
    """
    token = 'protein'

    prot_res = np.array([
        # CHARMM top_all27_prot_lipid.rtf
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HSD',
        'HSE', 'HSP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR',
        'TRP', 'TYR', 'VAL', 'ALAD',
        ## 'CHO','EAM', # -- special formyl and ethanolamine termini of gramicidin
        # PDB
        'HIS', 'MSE',
        # from Gromacs 4.5.3 oplsaa.ff/aminoacids.rtp
        'ARGN', 'ASPH', 'CYS2', 'CYSH', 'QLN', 'PGLU', 'GLUH', 'HIS1', 'HISD',
        'HISE', 'HISH', 'LYSH',
        # from Gromacs 4.5.3 gromos53a6.ff/aminoacids.rtp
        'ASN1', 'CYS1', 'HISA', 'HISB', 'HIS2',
        # from Gromacs 4.5.3 amber03.ff/aminoacids.rtp
        'HID', 'HIE', 'HIP', 'ORN', 'DAB', 'LYN', 'HYP', 'CYM', 'CYX', 'ASH',
        'GLH', 'ACE', 'NME',
        # from Gromacs 2016.3 amber99sb-star-ildn.ff/aminoacids.rtp
        'NALA', 'NGLY', 'NSER', 'NTHR', 'NLEU', 'NILE', 'NVAL', 'NASN', 'NGLN',
        'NARG', 'NHID', 'NHIE', 'NHIP', 'NTRP', 'NPHE', 'NTYR', 'NGLU', 'NASP',
        'NLYS', 'NPRO', 'NCYS', 'NCYX', 'NMET', 'CALA', 'CGLY', 'CSER', 'CTHR',
        'CLEU', 'CILE', 'CVAL', 'CASF', 'CASN', 'CGLN', 'CARG', 'CHID', 'CHIE',
        'CHIP', 'CTRP', 'CPHE', 'CTYR', 'CGLU', 'CASP', 'CLYS', 'CPRO', 'CCYS',
        'CCYX', 'CMET', 'CME', 'ASF',
    ])

    def __init__(self, parser, tokens):
        pass

    def apply(self, group):
        mask = np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class NucleicSelection(Selection):
    """All atoms in nucleic acid residues with recognized residue names.

    Recognized residue names:

    * from the CHARMM force field ::
        awk '/RESI/ {printf "'"'"%s"'"',",$2 }' top_all27_prot_na.rtf
    * recognized: 'ADE', 'URA', 'CYT', 'GUA', 'THY'
    * recognized (CHARMM in Gromacs): 'DA', 'DU', 'DC', 'DG', 'DT'

    .. versionchanged:: 0.8
       additional Gromacs selections
    """
    token = 'nucleic'

    nucl_res = np.array([
        'ADE', 'URA', 'CYT', 'GUA', 'THY', 'DA', 'DC', 'DG', 'DT', 'RA',
        'RU', 'RG', 'RC', 'A', 'T', 'U', 'C', 'G',
        'DA5', 'DC5', 'DG5', 'DT5',
        'DA3', 'DC3', 'DG3', 'DT3',
        'RA5', 'RU5', 'RG5', 'RC5',
        'RA3', 'RU3', 'RG3', 'RC3'
    ])

    def __init__(self, parser, tokens):
        pass

    def apply(self, group):
        mask = np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class BackboneSelection(ProteinSelection):
    """A BackboneSelection contains all atoms with name 'N', 'CA', 'C', 'O'.

    This excludes OT* on C-termini
    (which are included by, eg VMD's backbone selection).
    """
    token = 'backbone'
    bb_atoms = np.array(['N', 'CA', 'C', 'O'])

    def apply(self, group):
        mask = np.in1d(group.names, self.bb_atoms)
        mask &= np.in1d(group.resnames, self.prot_res)
        return group[mask].unique


class NucleicBackboneSelection(NucleicSelection):
    """Contains all atoms with name "P", "C5'", C3'", "O3'", "O5'".

    These atoms are only recognized if they are in a residue matched
    by the :class:`NucleicSelection`.
    """
    token = 'nucleicbackbone'
    bb_atoms = np.array(["P", "C5'", "C3'", "O3'", "O5'"])

    def apply(self, group):
        mask = np.in1d(group.names, self.bb_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class BaseSelection(NucleicSelection):
    """Selection of atoms in nucleobases.

    Recognized atom names (from CHARMM):

     'N9', 'N7', 'C8', 'C5', 'C4', 'N3', 'C2', 'N1', 'C6',
     'O6','N2','N6', 'O2','N4','O4','C5M'
    """
    token = 'nucleicbase'
    base_atoms = np.array([
        'N9', 'N7', 'C8', 'C5', 'C4', 'N3', 'C2', 'N1', 'C6',
        'O6', 'N2', 'N6',
        'O2', 'N4', 'O4', 'C5M'])

    def apply(self, group):
        mask = np.in1d(group.names, self.base_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class NucleicSugarSelection(NucleicSelection):
    """Contains all atoms with name C1', C2', C3', C4', O2', O4', O3'.
    """
    token = 'nucleicsugar'
    sug_atoms = np.array(["C1'", "C2'", "C3'", "C4'", "O4'"])

    def apply(self, group):
        mask = np.in1d(group.names, self.sug_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return group[mask].unique


class PropertySelection(Selection):
    """Some of the possible properties:
    x, y, z, radius, mass,
    """
    token = 'prop'
    ops = dict([
        ('>', np.greater),
        ('<', np.less),
        ('>=', np.greater_equal),
        ('<=', np.less_equal),
        ('==', np.equal),
        ('!=', np.not_equal),
    ])
    # order here is important, need to check <= before < so the
    # larger (in terms of string length) symbol is considered first
    _op_symbols = ('<=', '>=', '==', '!=', '<', '>')

    # symbols to replace with when flipping
    # eg 6 > x -> x <= 6, 5 == x -> x == 5
    opposite_ops = {
        '==': '==', '!=': '!=',
        '<': '>=', '>=': '<',
        '>': '<=', '<=': '>',
    }

    props = {'mass', 'charge', 'x', 'y', 'z'}

    def __init__(self, parser, tokens):
        """
        Possible splitting around operator:

        prop x < 5
        prop x< 5
        prop x <5
        prop x<5
        """
        prop = tokens.popleft()
        oper = None
        value = None
        if prop == "abs":
            self.absolute = True
            prop = tokens.popleft()
        else:
            self.absolute = False

        # check if prop has any extra information atm
        for possible in self._op_symbols:
            try:
                x, y = prop.split(possible)
            except ValueError:
                # won't unpack into 2 args unless *possible* is present
                pass
            else:
                prop = x
                oper = possible + y  # add back after splitting
                break

        if oper is None:
            oper = tokens.popleft()
        # check if oper has the value appended
        for possible in self._op_symbols:
            if possible in oper:
                x, y = oper.split(possible)
                if y:  # '<='.split('<=') == ['', ''], therefore y won't exist
                    oper = possible
                    value = y
                break

        if value is None:
            value = tokens.popleft()

        # check if we flip prop and value
        # eg 5 > x -> x <= 5
        if value in self.props:
            prop, value = value, prop
            oper = self.opposite_ops[oper]

        self.prop = prop
        try:
            self.operator = self.ops[oper]
        except KeyError:
            raise ValueError(
                "Invalid operator : '{0}' Use one of : '{1}'"
                "".format(oper, self.ops.keys()))
        self.value = float(value)

    def apply(self, group):
        try:
            col = {'x': 0, 'y': 1, 'z': 2}[self.prop]
        except KeyError:
            if self.prop == 'mass':
                values = group.masses
            elif self.prop == 'charge':
                values = group.charges
            else:
                raise SelectionError(
                    "Expected one of : {0}"
                    "".format(['x', 'y', 'z', 'mass', 'charge']))
        else:
            values = group.positions[:, col]

        if self.absolute:
            values = np.abs(values)
        mask = self.operator(values, self.value)

        return group[mask].unique


class SameSelection(Selection):
    token = 'same'
    precedence = 1

    prop_trans = {
        'fragment': None,
        'x': None,
        'y': None,
        'z': None,
        'residue': 'resids',
        'segment': 'segids',
        'name': 'names',
        'type': 'types',
        'resname': 'resnames',
        'resid': 'resids',
        'segid': 'segids',
        'mass': 'masses',
        'charge': 'charges',
        'radius': 'radii',
        'bfactor': 'bfactors',
        'resnum': 'resnums',
    }

    def __init__(self, parser, tokens):
        prop = tokens.popleft()
        if prop not in self.prop_trans:
            raise ValueError("Unknown same property : {0}"
                             "Choose one of : {1}"
                             "".format(prop, self.prop_trans.keys()))
        self.prop = prop
        parser.expect("as")
        self.sel = parser.parse_expression(self.precedence)
        self.prop = prop

    def apply(self, group):
        res = self.sel.apply(group)
        if not res:
            return group[[]]  # empty selection

        # Fragment must come before self.prop_trans lookups!
        if self.prop == 'fragment':
            # Combine all fragments together, then check where group
            # indices are same as fragment(s) indices
            allfrags = functools.reduce(lambda x, y: x + y, res.fragments)

            mask = np.in1d(group.indices, allfrags.indices)
            return group[mask].unique
        # [xyz] must come before self.prop_trans lookups too!
        try:
            pos_idx = {'x': 0, 'y': 1, 'z': 2}[self.prop]
        except KeyError:
            # The self.prop string was already checked,
            # so don't need error checking here.
            # KeyError at this point is impossible!
            attrname = self.prop_trans[self.prop]
            vals = getattr(res, attrname)
            mask = np.in1d(getattr(group, attrname), vals)

            return group[mask].unique
        else:
            vals = res.positions[:, pos_idx]
            pos = group.positions[:, pos_idx]

            # isclose only does one value at a time
            mask = np.vstack([np.isclose(pos, v)
                              for v in vals]).any(axis=0)
            return group[mask].unique


class SelectionParser(object):
    """A small parser for selection expressions.  Demonstration of
    recursive descent parsing using Precedence climbing (see
    http://www.engr.mun.ca/~theo/Misc/exp_parsing.htm).  Transforms
    expressions into nested Selection tree.

    For reference, the grammar that we parse is ::

       E(xpression)--> Exp(0)
       Exp(p) -->      P {B Exp(q)}
       P -->           U Exp(q) | "(" E ")" | v
       B(inary) -->    "and" | "or"
       U(nary) -->     "not"
       T(erms) -->     segid [value]
                       | resname [value]
                       | resid [value]
                       | name [value]
                       | type [value]
   """
    # Borg pattern: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66531
    _shared_state = {}

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def expect(self, token):
        """Anticipate and remove a given token"""
        if self.tokens[0] == token:
            self.tokens.popleft()
        else:
            raise SelectionError(
                "Unexpected token: '{0}' Expected: '{1}'"
                "".format(self.tokens[0], token))

    def parse(self, selectstr, selgroups):
        """Create a Selection object from a string.

        Parameters
        ----------
        selectstr : str
            The string that describes the selection
        selgroups : AtomGroups
            AtomGroups to be used in `group` selections

        Returns
        -------
        The appropriate Selection object.  Use the .apply method on
        this to perform the selection.

        Raises
        ------
        SelectionError
            If anything goes wrong in creating the Selection object.
        """
        self.selectstr = selectstr
        self.selgroups = selgroups
        tokens = selectstr.replace('(', ' ( ').replace(')', ' ) ')
        self.tokens = collections.deque(tokens.split() + [None])
        parsetree = self.parse_expression(0)
        if self.tokens[0] is not None:
            raise SelectionError(
                "Unexpected token at end of selection string: '{0}'"
                "".format(self.tokens[0]))
        return parsetree

    def parse_expression(self, p):
        exp1 = self._parse_subexp()
        while (self.tokens[0] in _OPERATIONS and
               _OPERATIONS[self.tokens[0]].precedence >= p):
            op = _OPERATIONS[self.tokens.popleft()]
            q = 1 + op.precedence
            exp2 = self.parse_expression(q)
            exp1 = op(exp1, exp2)
        return exp1

    def _parse_subexp(self):
        op = self.tokens.popleft()

        if op == '(':
            exp = self.parse_expression(0)
            self.expect(')')
            return exp

        try:
            return _SELECTIONDICT[op](self, self.tokens)
        except KeyError:
            raise SelectionError("Unknown selection token: '{0}'".format(op))
        except ValueError as e:
            raise SelectionError("Selection failed: '{0}'".format(e))


# The module level instance
Parser = SelectionParser()
