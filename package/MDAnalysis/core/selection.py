# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
import collections
import re
import fnmatch
import functools
import sys
import types
import warnings

import numpy as np


from ..lib.util import unique_int_1d
from ..lib import distances
from ..exceptions import SelectionError, NoDataError, SelectionWarning

#: Regular expression for recognizing a floating point number in a selection.
#: Numbers such as 1.2, 1.2e-01, -1.2 are all parsed as Python floats.
FLOAT_PATTERN = r"-?\d*\.?\d*(?:e[-+]?\d+)?"

#: Regular expression for recognizing un/signed integers in a selection.
INT_PATTERN = r"-?\d+"

#: Regular expression for recognising a range separator.
#: Delimiters include ":", "-", "to" and can have arbitrary whitespace.
RANGE_PATTERN = r"\s*(?:[:-]| to )\s*"


def is_keyword(val):
    """Is val a selection keyword?

    Returns False on any of the following strings:
      - keys in SELECTIONDICT (tokens from Selection objects)
      - keys in OPERATIONS (tokens from LogicOperations)
      - (Parentheses)
      - The value `None` (used as EOF in selection strings)
    """
    return (
        val in _SELECTIONDICT
        or val in _OPERATIONS
        or val in ["(", ")"]
        or val is None
    )


def grab_not_keywords(tokens):
    """Pop tokens from the left until you hit a keyword

    Keywords can be escaped with a backslash to allow their use as names,
    e.g. '\\protein'

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
        # Remove escape sequence to allow use of keywords as names
        values.append(val.removeprefix("\\"))
    return values


def join_separated_values(values):
    """Join range values that are separated by whitespace

    Parameters
    ----------
    values: list
        list of value strings

    Returns
    -------
    values: list of strings

    Examples
    --------
    join_separated_values(['37', 'to', '22'])
    >>> ['37 to 22']

    .. versionadded:: 2.0.0
    """
    _values = []
    DELIMITERS = ("to", ":", "-")
    while values:
        v = values.pop(0)

        if v in DELIMITERS:
            try:
                _values[-1] = f"{_values[-1]} {v} {values.pop(0)}"
            except IndexError:
                given = f"{' '.join(_values)} {v} {' '.join(values)}"
                raise SelectionError(f"Invalid expression given: {given}")
        elif _values and (
            v[:2] in ("--", "to")
            or v[0] == ":"
            or any(_values[-1].endswith(x) for x in DELIMITERS)
        ):
            _values[-1] = f"{_values[-1]} {v}"
        else:
            _values.append(v)
    return _values


_SELECTIONDICT = {}
_OPERATIONS = {}
# These are named args to select_atoms that have a special meaning and must
# not be allowed as names for the 'group' keyword.
_RESERVED_KWARGS = ("updating",)


# And and Or are exception and aren't strictly a Selection
# as they work on other Selections rather than doing work themselves.
# So their init is a little strange too....
class _Operationmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            _OPERATIONS[classdict["token"]] = cls
        except KeyError:
            pass


class LogicOperation(object, metaclass=_Operationmeta):
    def __init__(self, lsel, rsel, parser):
        self.rsel = rsel
        self.lsel = lsel
        self.parser = parser

    def apply(self, *args, **kwargs):
        return self._apply(*args, **kwargs).asunique(sorted=self.parser.sorted)


class AndOperation(LogicOperation):
    token = "and"
    precedence = 3

    def _apply(self, group):
        rsel = self.rsel.apply(group)
        lsel = self.lsel.apply(group)

        # Mask which lsel indices appear in rsel
        mask = np.isin(rsel.indices, lsel.indices)
        # and mask rsel according to that
        return rsel[mask]


class OrOperation(LogicOperation):
    token = "or"
    precedence = 3

    def _apply(self, group):
        lsel = self.lsel.apply(group)
        rsel = self.rsel.apply(group)

        # Find unique indices from both these AtomGroups
        # and slice main list using them
        idx = np.union1d(lsel.indices, rsel.indices).astype(np.int32)

        return group.universe.atoms[idx]


def return_empty_on_apply(func):
    """
    Decorator to return empty AtomGroups from the apply() function
    without evaluating it
    """

    @functools.wraps(func)
    def _apply(self, group):
        if len(group) == 0:
            return group
        return func(self, group)

    return _apply


class _Selectionmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            _SELECTIONDICT[classdict["token"]] = cls
            _SELECTIONDICT[classdict["token"].lower()] = cls
        except KeyError:
            pass


class Selection(object, metaclass=_Selectionmeta):

    def __init__(self, parser, tokens):
        self.parser = parser

    def apply(self, *args, **kwargs):
        return self._apply(*args, **kwargs).asunique(sorted=self.parser.sorted)


class AllSelection(Selection):
    token = "all"

    def _apply(self, group):
        # Check whether group is identical to the one stored
        # in the corresponding universe, in which case this
        # is returned directly. This works since the Universe.atoms
        # are unique by construction.
        if group is group.universe.atoms:
            return group
        return group[:]


class UnarySelection(Selection):
    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        sel = parser.parse_expression(self.precedence)
        self.sel = sel


class NotSelection(UnarySelection):
    token = "not"
    precedence = 5

    def _apply(self, group):
        notsel = self.sel.apply(group)
        return group[~np.isin(group.indices, notsel.indices)]


class GlobalSelection(UnarySelection):
    token = "global"
    precedence = 5

    def _apply(self, group):
        return self.sel.apply(group.universe.atoms).unique


class ByResSelection(UnarySelection):
    """
    Selects all atoms that are in the same segment and residue as selection

    .. versionchanged:: 1.0.0
       Use :code:`"resindices"` instead of :code:`"resids"` (see #2669 and #2672)
    """

    token = "byres"
    precedence = 1

    def _apply(self, group):
        res = self.sel.apply(group)
        unique_res = unique_int_1d(res.resindices)
        mask = np.isin(group.resindices, unique_res)

        return group[mask]


class AroundSelection(Selection):
    token = "around"
    precedence = 1

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.periodic = parser.periodic
        self.cutoff = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)

    @return_empty_on_apply
    def _apply(self, group):
        indices = []
        sel = self.sel.apply(group)
        # All atoms in group that aren't in sel
        sys = group[~np.isin(group.indices, sel.indices)]

        if not sys or not sel:
            return sys[[]]

        box = group.dimensions if self.periodic else None
        pairs = distances.capped_distance(
            sel.positions,
            sys.positions,
            self.cutoff,
            box=box,
            return_distances=False,
        )
        if pairs.size > 0:
            indices = np.sort(pairs[:, 1])

        return sys[np.asarray(indices, dtype=np.int64)]


class SphericalLayerSelection(Selection):
    token = "sphlayer"
    precedence = 1

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.periodic = parser.periodic
        self.inRadius = float(tokens.popleft())
        self.exRadius = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)

    @return_empty_on_apply
    def _apply(self, group):
        indices = []
        sel = self.sel.apply(group)
        if len(sel) == 0:
            return group[[]]

        box = group.dimensions if self.periodic else None
        ref = sel.center_of_geometry().reshape(1, 3).astype(np.float32)
        pairs = distances.capped_distance(
            ref,
            group.positions,
            self.exRadius,
            min_cutoff=self.inRadius,
            box=box,
            return_distances=False,
        )
        if pairs.size > 0:
            indices = np.sort(pairs[:, 1])

        return group[np.asarray(indices, dtype=np.int64)]


class IsoLayerSelection(Selection):
    token = "isolayer"
    precedence = 1

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.periodic = parser.periodic
        self.inRadius = float(tokens.popleft())
        self.exRadius = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)

    @return_empty_on_apply
    def _apply(self, group):
        indices = []
        sel = self.sel.apply(group)
        # All atoms in group that aren't in sel
        sys = group[~np.isin(group.indices, sel.indices)]

        if not sys or not sel:
            return sys[[]]

        box = group.dimensions if self.periodic else None
        pairs_outer = distances.capped_distance(
            sel.positions,
            sys.positions,
            self.exRadius,
            box=box,
            return_distances=False,
        )
        pairs_inner = distances.capped_distance(
            sel.positions,
            sys.positions,
            self.inRadius,
            box=box,
            return_distances=False,
        )

        if pairs_outer.size > 0:
            sys_ind_outer = np.sort(np.unique(pairs_outer[:, 1]))
            if pairs_inner.size > 0:
                sys_ind_inner = np.sort(np.unique(pairs_inner[:, 1]))
                indices = sys_ind_outer[~np.isin(sys_ind_outer, sys_ind_inner)]
            else:
                indices = sys_ind_outer

        return sys[np.asarray(indices, dtype=np.int64)]


class SphericalZoneSelection(Selection):
    token = "sphzone"
    precedence = 1

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.periodic = parser.periodic
        self.cutoff = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)

    @return_empty_on_apply
    def _apply(self, group):
        indices = []
        sel = self.sel.apply(group)
        if len(sel) == 0:
            return group[[]]

        box = group.dimensions if self.periodic else None
        ref = sel.center_of_geometry().reshape(1, 3).astype(np.float32)
        pairs = distances.capped_distance(
            ref, group.positions, self.cutoff, box=box, return_distances=False
        )
        if pairs.size > 0:
            indices = np.sort(pairs[:, 1])

        return group[np.asarray(indices, dtype=np.int64)]


class CylindricalSelection(Selection):
    @return_empty_on_apply
    def _apply(self, group):
        sel = self.sel.apply(group)
        if len(sel) == 0:
            return group[[]]
        # Calculate vectors between point of interest and our group
        vecs = group.positions - sel.center_of_geometry()

        if self.periodic and not group.dimensions is None:
            box = group.dimensions[:3]
            cyl_z_hheight = self.zmax - self.zmin

            if 2 * self.exRadius > box[0]:
                raise NotImplementedError(
                    "The diameter of the cylinder selection ({:.3f}) is larger "
                    "than the unit cell's x dimension ({:.3f}). Can only do "
                    "selections where it is smaller or equal."
                    "".format(2 * self.exRadius, box[0])
                )
            if 2 * self.exRadius > box[1]:
                raise NotImplementedError(
                    "The diameter of the cylinder selection ({:.3f}) is larger "
                    "than the unit cell's y dimension ({:.3f}). Can only do "
                    "selections where it is smaller or equal."
                    "".format(2 * self.exRadius, box[1])
                )
            if cyl_z_hheight > box[2]:
                raise NotImplementedError(
                    "The total length of the cylinder selection in z ({:.3f}) "
                    "is larger than the unit cell's z dimension ({:.3f}). Can "
                    "only do selections where it is smaller or equal."
                    "".format(cyl_z_hheight, box[2])
                )

            if np.all(group.dimensions[3:] == 90.0):
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
        radii = vecs[:, 0] ** 2 + vecs[:, 1] ** 2
        mask = radii < self.exRadius**2
        try:
            mask &= radii > self.inRadius**2
        except AttributeError:
            # Only for cylayer, cyzone doesn't have inRadius
            pass

        return group[mask]


class CylindricalZoneSelection(CylindricalSelection):
    token = "cyzone"
    precedence = 1

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.periodic = parser.periodic
        self.exRadius = float(tokens.popleft())
        self.zmax = float(tokens.popleft())
        self.zmin = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)


class CylindricalLayerSelection(CylindricalSelection):
    token = "cylayer"
    precedence = 1

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.periodic = parser.periodic
        self.inRadius = float(tokens.popleft())
        self.exRadius = float(tokens.popleft())
        self.zmax = float(tokens.popleft())
        self.zmin = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)


class PointSelection(Selection):
    token = "point"

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.periodic = parser.periodic
        x = float(tokens.popleft())
        y = float(tokens.popleft())
        z = float(tokens.popleft())
        self.ref = np.array([x, y, z], dtype=np.float32)
        self.cutoff = float(tokens.popleft())

    @return_empty_on_apply
    def _apply(self, group):
        indices = []

        box = group.dimensions if self.periodic else None
        pairs = distances.capped_distance(
            self.ref[None, :],
            group.positions,
            self.cutoff,
            box=box,
            return_distances=False,
        )
        if pairs.size > 0:
            indices = np.sort(pairs[:, 1])

        return group[np.asarray(indices, dtype=np.int64)]


class AtomSelection(Selection):
    token = "atom"

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.segid = tokens.popleft()
        self.resid = int(tokens.popleft())
        self.name = tokens.popleft()

    def _apply(self, group):
        sub = group[group.names == self.name]
        if sub:
            sub = sub[sub.resids == self.resid]
        if sub:
            sub = sub[sub.segids == self.segid]
        return sub.unique


class BondedSelection(Selection):
    token = "bonded"
    precedence = 1

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.sel = parser.parse_expression(self.precedence)

    def _apply(self, group):
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
        idx.append(bix[:, 0][np.isin(bix[:, 1], grpidx)])
        # right side
        idx.append(bix[:, 1][np.isin(bix[:, 0], grpidx)])

        idx = np.union1d(*idx)

        return group.universe.atoms[np.unique(idx)]


class SelgroupSelection(Selection):
    token = "group"

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        grpname = tokens.popleft()
        if grpname in _RESERVED_KWARGS:
            raise TypeError(
                "The '{}' keyword is reserved and cannot be "
                "used as a selection group name.".format(grpname)
            )
        try:
            self.grp = parser.selgroups[grpname]
        except KeyError:
            errmsg = f"Failed to find group: {grpname}"
            raise ValueError(errmsg) from None

    def _apply(self, group):
        mask = np.isin(group.indices, self.grp.indices)
        return group[mask]


class SingleCharSelection(Selection):
    """for when an attribute is just a single character, eg RSChirality

    .. versionadded:: 2.1.0
    """

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        vals = grab_not_keywords(tokens)
        if not vals:
            raise ValueError("Unexpected token '{0}'".format(tokens[0]))

        self.values = vals

    @return_empty_on_apply
    def _apply(self, group):
        attr = getattr(group, self.field)

        mask = np.isin(attr, self.values)

        return group[mask]


class _ProtoStringSelection(Selection):
    """Selections based on text attributes

    .. versionchanged:: 1.0.0
        Supports multiple wildcards, based on fnmatch
    """

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        vals = grab_not_keywords(tokens)
        if not vals:
            raise ValueError("Unexpected token '{0}'".format(tokens[0]))

        self.values = vals

    @return_empty_on_apply
    def _apply(self, group):
        # rather than work on group.names, cheat and look at the lookup table
        nmattr = getattr(group.universe._topology, self.field)

        matches = []  # list of passing indices
        # iterate through set of known atom names, check which pass
        for nm, ix in nmattr.namedict.items():
            if any(fnmatch.fnmatchcase(nm, val) for val in self.values):
                matches.append(ix)

        # atomname indices for members of this group
        nmidx = nmattr.nmidx[getattr(group, self.level)]

        return group[np.isin(nmidx, matches)]


class AromaticSelection(Selection):
    """Select aromatic atoms.

    Aromaticity is available in the `aromaticities` attribute and is made
    available through RDKit"""

    token = "aromatic"
    field = "aromaticities"

    def _apply(self, group):
        return group[group.aromaticities]


class SmartsSelection(Selection):
    """Select atoms based on SMARTS queries.

    Uses RDKit to run the query and converts the result to MDAnalysis.
    Supports chirality.

    .. versionchanged:: 2.2.0
       ``rdkit_wargs`` and ``smarts_kwargs`` can now be passed to control
       the behaviour of the RDKit converter and RDKit's ``GetSubstructMatches``
       respectively.
       The default ``maxMatches`` value passed to ``GetSubstructMatches`` has
       been changed from ``1000`` to ``max(1000, n_atoms * 10)`` in order to
       limit cases where too few matches were generated. A warning is now also
       thrown if ``maxMatches`` has been reached.
    """

    token = "smarts"

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        # The parser will add spaces around parentheses and then split the
        # selection based on spaces to create the tokens
        # If the input SMARTS query contained parentheses, the query will be
        # split because of that and we need to reconstruct it
        # We also need to keep the parentheses that are not part of the smarts
        # query intact
        pattern = []
        counter = {"(": 0, ")": 0}
        # loop until keyword but ignore parentheses as a keyword
        while tokens[0] in ["(", ")"] or not is_keyword(tokens[0]):
            # keep track of the number of open and closed parentheses
            if tokens[0] in ["(", ")"]:
                counter[tokens[0]] += 1
                # if the char is a closing ")" but there's no corresponding
                # open "(" then we've reached then end of the smarts query and
                # the current token ")" is part of a grouping parenthesis
                if tokens[0] == ")" and counter["("] < (counter[")"]):
                    break
            # add the token to the pattern and remove it from the tokens
            val = tokens.popleft()
            pattern.append(val)
        self.pattern = "".join(pattern)
        self.rdkit_kwargs = parser.rdkit_kwargs
        self.smarts_kwargs = parser.smarts_kwargs

    def _apply(self, group):
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError(
                "RDKit is required for SMARTS-based atom "
                "selection but it's not installed. Try "
                "installing it with \n"
                "conda install -c conda-forge rdkit"
            )
        pattern = Chem.MolFromSmarts(self.pattern)
        if not pattern:
            raise ValueError(f"{self.pattern!r} is not a valid SMARTS query")
        mol = group.convert_to("RDKIT", **self.rdkit_kwargs)
        self.smarts_kwargs.setdefault("useChirality", True)
        self.smarts_kwargs.setdefault("maxMatches", max(1000, len(group) * 10))
        matches = mol.GetSubstructMatches(pattern, **self.smarts_kwargs)
        if len(matches) == self.smarts_kwargs["maxMatches"]:
            warnings.warn(
                "Your smarts-based atom selection returned the max"
                "number of matches. This indicates that not all"
                "matching atoms were selected. When calling"
                "atom_group.select_atoms(), the default value"
                "of maxMatches is max(100, len(atom_group * 10)). "
                "To fix this, add the following argument to "
                "select_atoms: \n"
                "smarts_kwargs={maxMatches: <higher_value>}"
            )
        # flatten all matches and remove duplicated indices
        indices = np.unique([idx for match in matches for idx in match])
        # create boolean mask for atoms based on index
        mask = np.isin(range(group.n_atoms), indices)
        return group[mask]


class ResidSelection(Selection):
    """Select atoms based on numerical fields

    Allows the use of ':', '-' and 'to' to specify a range of values
    For example

      resid 1:10
    """

    token = "resid"

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        values = grab_not_keywords(tokens)
        if not values:
            raise ValueError("Unexpected token: '{0}'".format(tokens[0]))

        values = join_separated_values(values)

        # each value in uppers and lowers is a tuple of (resid, icode)
        uppers = []
        lowers = []

        for val in values:
            m1 = re.match(f"({INT_PATTERN})" + r"(\w?)$", val)
            if not m1 is None:
                res = m1.groups()
                lower = int(res[0]), res[1]
                upper = None, None
            else:
                # check if in appropriate format 'lower:upper' or 'lower-upper'
                # each val is one or more digits, maybe a letter
                pattern = f"({INT_PATTERN})(\\w?){RANGE_PATTERN}"
                pattern += f"({INT_PATTERN})(\\w?)"
                selrange = re.match(pattern, val)
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

    def _apply(self, group):
        # Grab arrays here to reduce number of calls to main topology
        vals = group.resids
        try:  # optional attribute
            icodes = group.icodes
        except (AttributeError, NoDataError):
            icodes = None
            # if no icodes and icodes are part of selection, cause a fuss
            if any(v[1] for v in self.uppers) or any(
                v[1] for v in self.lowers
            ):
                errmsg = (
                    "Selection specified icodes, while the topology "
                    "doesn't have any."
                )
                raise ValueError(errmsg) from None

        if not icodes is None:
            mask = self._sel_with_icodes(vals, icodes)
        else:
            mask = self._sel_without_icodes(vals)

        return group[mask]

    def _sel_without_icodes(self, vals):
        # Final mask that gets applied to group
        mask = np.zeros(len(vals), dtype=bool)

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
        mask = np.zeros(len(vals), dtype=bool)

        for (u_resid, u_icode), (l_resid, l_icode) in zip(
            self.uppers, self.lowers
        ):
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


class BoolSelection(Selection):
    """Selection for boolean values"""

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        values = grab_not_keywords(tokens)
        if not values:
            values = ["true"]

        self.values = []
        for val in values:
            lower = val.lower()
            if lower == "false":
                bval = False
            elif lower == "true":
                bval = True
            else:
                raise ValueError(
                    f"'{val}' is an invalid value "
                    "for boolean selection. "
                    "Use 'True' or 'False'"
                )
            self.values.append(bval)

    def _apply(self, group):
        vals = getattr(group, self.field)
        mask = np.zeros(len(vals), dtype=bool)
        for val in self.values:
            mask |= vals == val
        return group[mask]


class RangeSelection(Selection):
    """Range selection for int values"""

    value_offset = 0
    pattern = f"({INT_PATTERN}){RANGE_PATTERN}({INT_PATTERN})"
    dtype = int

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)
        self.rtol = parser.rtol
        self.atol = parser.atol

        values = grab_not_keywords(tokens)
        if not values:
            raise ValueError("Unexpected token: '{0}'".format(tokens[0]))

        values = join_separated_values(values)

        uppers = []  # upper limit on any range
        lowers = []  # lower limit on any range

        for val in values:
            try:
                lower = self.dtype(val)
                upper = None
            except ValueError:
                # check if in appropriate format 'lower:upper' or 'lower-upper'
                selrange = re.match(self.pattern, val)
                if not selrange:
                    errmsg = f"Failed to parse number: {val}"
                    raise ValueError(errmsg) from None
                lower, upper = map(self.dtype, selrange.groups())

            lowers.append(lower)
            uppers.append(upper)

        self.lowers = lowers
        self.uppers = uppers

    def _apply(self, group):
        mask = np.zeros(len(group), dtype=bool)
        vals = getattr(group, self.field) + self.value_offset

        for upper, lower in zip(self.uppers, self.lowers):
            if upper is not None:
                thismask = vals >= lower
                thismask &= vals <= upper
            else:
                thismask = vals == lower

            mask |= thismask
        return group[mask]


class FloatRangeSelection(RangeSelection):
    """Range selection for float values"""

    pattern = f"({FLOAT_PATTERN}){RANGE_PATTERN}({FLOAT_PATTERN})"
    dtype = float

    def _apply(self, group):
        mask = np.zeros(len(group), dtype=bool)
        vals = getattr(group, self.field) + self.value_offset

        for upper, lower in zip(self.uppers, self.lowers):
            if upper is not None:
                thismask = vals >= lower
                thismask &= vals <= upper
            else:
                low, high = lower - 1, lower + 1
                fp = "https://docs.python.org/3.8/tutorial/floatingpoint.html"
                msg = (
                    "Using float equality to select atoms is "
                    "not recommended because of inherent "
                    "limitations in representing numbers on "
                    f"computers (see {fp} for more). "
                    "Instead, we recommend using a range "
                    f"to select, e.g. '{self.token} {low} to {high}'. "
                    "If you still want to compare floats, use the "
                    "`atol` and `rtol` keywords to modify the tolerance "
                    "for `np.isclose`."
                )
                warnings.warn(msg, category=SelectionWarning)
                thismask = np.isclose(
                    vals, lower, atol=self.atol, rtol=self.rtol
                )

            mask |= thismask
        return group[mask]


class ByNumSelection(RangeSelection):
    token = "bynum"
    field = "indices"
    value_offset = 1  # queries are in 1 based indices


class ProteinSelection(Selection):
    """Consists of all residues with  recognized residue names.

    Recognized residue names in :attr:`ProteinSelection.prot_res`.

      * from the CHARMM force field::
         awk '/RESI/ {printf "'"'"%s"'"',",$2 }' top_all27_prot_lipid.rtf

      * manually added special CHARMM, OPLS/AA and Amber residue names.

    See Also
    --------
    :func:`MDAnalysis.lib.util.convert_aa_code`


    .. versionchanged:: 1.0.1
       prot_res changed to set (from numpy array)
       performance improved by ~100x on larger systems
    """

    token = "protein"

    prot_res = {
        # CHARMM top_all27_prot_lipid.rtf
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HSD",
        "HSE",
        "HSP",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "ALAD",
        ## 'CHO','EAM', # -- special formyl and ethanolamine termini of gramicidin
        # PDB
        "HIS",
        "MSE",
        # from Gromacs 4.5.3 oplsaa.ff/aminoacids.rtp
        "ARGN",
        "ASPH",
        "CYS2",
        "CYSH",
        "QLN",
        "PGLU",
        "GLUH",
        "HIS1",
        "HISD",
        "HISE",
        "HISH",
        "LYSH",
        # from Gromacs 4.5.3 gromos53a6.ff/aminoacids.rtp
        "ASN1",
        "CYS1",
        "HISA",
        "HISB",
        "HIS2",
        # from Gromacs 4.5.3 amber03.ff/aminoacids.rtp
        "HID",
        "HIE",
        "HIP",
        "ORN",
        "DAB",
        "LYN",
        "HYP",
        "CYM",
        "CYX",
        "ASH",
        "GLH",
        "ACE",
        "NME",
        # from Gromacs 2016.3 amber99sb-star-ildn.ff/aminoacids.rtp
        "NALA",
        "NGLY",
        "NSER",
        "NTHR",
        "NLEU",
        "NILE",
        "NVAL",
        "NASN",
        "NGLN",
        "NARG",
        "NHID",
        "NHIE",
        "NHIP",
        "NTRP",
        "NPHE",
        "NTYR",
        "NGLU",
        "NASP",
        "NLYS",
        "NPRO",
        "NCYS",
        "NCYX",
        "NMET",
        "CALA",
        "CGLY",
        "CSER",
        "CTHR",
        "CLEU",
        "CILE",
        "CVAL",
        "CASF",
        "CASN",
        "CGLN",
        "CARG",
        "CHID",
        "CHIE",
        "CHIP",
        "CTRP",
        "CPHE",
        "CTYR",
        "CGLU",
        "CASP",
        "CLYS",
        "CPRO",
        "CCYS",
        "CCYX",
        "CMET",
        "CME",
        "ASF",
    }

    def _apply(self, group):
        resname_attr = group.universe._topology.resnames
        # which values in resname attr are in prot_res?
        matches = [
            ix
            for (nm, ix) in resname_attr.namedict.items()
            if nm in self.prot_res
        ]
        # index of each atom's resname
        nmidx = resname_attr.nmidx[group.resindices]
        # intersect atom's resname index and matches to prot_res
        return group[np.isin(nmidx, matches)]


class NucleicSelection(Selection):
    """All atoms in nucleic acid residues with recognized residue names.

    Recognized residue names:

    * from the CHARMM force field ::
        awk '/RESI/ {printf "'"'"%s"'"',",$2 }' top_all27_prot_na.rtf
    * recognized: 'ADE', 'URA', 'CYT', 'GUA', 'THY'
    * recognized (CHARMM in Gromacs): 'DA', 'DU', 'DC', 'DG', 'DT'

    .. versionchanged:: 0.8
       additional Gromacs selections
    .. versionchanged:: 1.0.1
       nucl_res changed to set (from numpy array)
       performance improved by ~100x on larger systems
    """

    token = "nucleic"

    nucl_res = {
        "ADE",
        "URA",
        "CYT",
        "GUA",
        "THY",
        "DA",
        "DC",
        "DG",
        "DT",
        "RA",
        "RU",
        "RG",
        "RC",
        "A",
        "T",
        "U",
        "C",
        "G",
        "DA5",
        "DC5",
        "DG5",
        "DT5",
        "DA3",
        "DC3",
        "DG3",
        "DT3",
        "RA5",
        "RU5",
        "RG5",
        "RC5",
        "RA3",
        "RU3",
        "RG3",
        "RC3",
    }

    def _apply(self, group):
        resnames = group.universe._topology.resnames
        nmidx = resnames.nmidx[group.resindices]

        matches = [
            ix for (nm, ix) in resnames.namedict.items() if nm in self.nucl_res
        ]
        mask = np.isin(nmidx, matches)

        return group[mask]


class WaterSelection(Selection):
    """All atoms in water residues with recognized water residue names.

    Recognized residue names:

    * recognized 3 Letter resnames: 'H2O', 'HOH', 'OH2', 'HHO', 'OHH'
        'TIP', 'T3P', 'T4P', 'T5P', 'SOL', 'WAT'

    * recognized 4 Letter resnames: 'TIP2', 'TIP3', 'TIP4'

    .. versionadded:: 2.9.0
    """

    token = "water"

    # Recognized water resnames
    water_res = {
        "H2O",
        "HOH",
        "OH2",
        "HHO",
        "OHH",
        "T3P",
        "T4P",
        "T5P",
        "SOL",
        "WAT",
        "TIP",
        "TIP2",
        "TIP3",
        "TIP4",
    }

    def _apply(self, group):
        resnames = group.universe._topology.resnames
        nmidx = resnames.nmidx[group.resindices]

        matches = [
            ix
            for (nm, ix) in resnames.namedict.items()
            if nm in self.water_res
        ]
        mask = np.isin(nmidx, matches)

        return group[mask]


class BackboneSelection(ProteinSelection):
    """A BackboneSelection contains all atoms with name 'N', 'CA', 'C', 'O'.

    This excludes OT* on C-termini
    (which are included by, eg VMD's backbone selection).


    .. versionchanged:: 1.0.1
       bb_atoms changed to set (from numpy array)
       performance improved by ~100x on larger systems
    """

    token = "backbone"
    bb_atoms = {"N", "CA", "C", "O"}

    def _apply(self, group):
        atomnames = group.universe._topology.names
        resnames = group.universe._topology.resnames

        # filter by atom names
        name_matches = [
            ix
            for (nm, ix) in atomnames.namedict.items()
            if nm in self.bb_atoms
        ]
        nmidx = atomnames.nmidx[group.ix]
        group = group[np.isin(nmidx, name_matches)]

        # filter by resnames
        resname_matches = [
            ix for (nm, ix) in resnames.namedict.items() if nm in self.prot_res
        ]
        nmidx = resnames.nmidx[group.resindices]
        group = group[np.isin(nmidx, resname_matches)]

        return group.unique


class NucleicBackboneSelection(NucleicSelection):
    """Contains all atoms with name "P", "C5'", C3'", "O3'", "O5'".

    These atoms are only recognized if they are in a residue matched
    by the :class:`NucleicSelection`.


    .. versionchanged:: 1.0.1
       bb_atoms changed to set (from numpy array)
       performance improved by ~100x on larger systems
    """

    token = "nucleicbackbone"
    bb_atoms = {"P", "C5'", "C3'", "O3'", "O5'"}

    def _apply(self, group):
        atomnames = group.universe._topology.names
        resnames = group.universe._topology.resnames

        # filter by atom names
        name_matches = [
            ix
            for (nm, ix) in atomnames.namedict.items()
            if nm in self.bb_atoms
        ]
        nmidx = atomnames.nmidx[group.ix]
        group = group[np.isin(nmidx, name_matches)]

        # filter by resnames
        resname_matches = [
            ix for (nm, ix) in resnames.namedict.items() if nm in self.nucl_res
        ]
        nmidx = resnames.nmidx[group.resindices]
        group = group[np.isin(nmidx, resname_matches)]

        return group.unique


class BaseSelection(NucleicSelection):
    """Selection of atoms in nucleobases.

    Recognized atom names (from CHARMM):

     'N9', 'N7', 'C8', 'C5', 'C4', 'N3', 'C2', 'N1', 'C6',
     'O6','N2','N6', 'O2','N4','O4','C5M'


    .. versionchanged:: 1.0.1
       base_atoms changed to set (from numpy array)
       performance improved by ~100x on larger systems
    """

    token = "nucleicbase"
    base_atoms = {
        "N9",
        "N7",
        "C8",
        "C5",
        "C4",
        "N3",
        "C2",
        "N1",
        "C6",
        "O6",
        "N2",
        "N6",
        "O2",
        "N4",
        "O4",
        "C5M",
    }

    def _apply(self, group):
        atomnames = group.universe._topology.names
        resnames = group.universe._topology.resnames

        # filter by atom names
        name_matches = [
            ix
            for (nm, ix) in atomnames.namedict.items()
            if nm in self.base_atoms
        ]
        nmidx = atomnames.nmidx[group.ix]
        group = group[np.isin(nmidx, name_matches)]

        # filter by resnames
        resname_matches = [
            ix for (nm, ix) in resnames.namedict.items() if nm in self.nucl_res
        ]
        nmidx = resnames.nmidx[group.resindices]
        group = group[np.isin(nmidx, resname_matches)]

        return group.unique


class NucleicSugarSelection(NucleicSelection):
    """Contains all atoms with name C1', C2', C3', C4', O2', O4', O3'.


    .. versionchanged:: 1.0.1
       sug_atoms changed to set (from numpy array)
       performance improved by ~100x on larger systems
    """

    token = "nucleicsugar"
    sug_atoms = {"C1'", "C2'", "C3'", "C4'", "O4'"}

    def _apply(self, group):
        atomnames = group.universe._topology.names
        resnames = group.universe._topology.resnames

        # filter by atom names
        name_matches = [
            ix
            for (nm, ix) in atomnames.namedict.items()
            if nm in self.sug_atoms
        ]
        nmidx = atomnames.nmidx[group.ix]
        group = group[np.isin(nmidx, name_matches)]

        # filter by resnames
        resname_matches = [
            ix for (nm, ix) in resnames.namedict.items() if nm in self.nucl_res
        ]
        nmidx = resnames.nmidx[group.resindices]
        group = group[np.isin(nmidx, resname_matches)]

        return group.unique


class PropertySelection(Selection):
    """Some of the possible properties:
    x, y, z, radius, mass,

    .. versionchanged:: 2.0.0
        changed == operator to use np.isclose instead of np.equals.
        Added ``atol`` and ``rtol`` keywords to control ``np.isclose``
        tolerance.
    """

    token = "prop"
    ops = dict(
        [
            (">", np.greater),
            ("<", np.less),
            (">=", np.greater_equal),
            ("<=", np.less_equal),
            ("==", np.isclose),
            ("!=", np.not_equal),
        ]
    )
    # order here is important, need to check <= before < so the
    # larger (in terms of string length) symbol is considered first
    _op_symbols = ("<=", ">=", "==", "!=", "<", ">")

    # symbols to replace with when flipping
    # eg 6 > x -> x <= 6, 5 == x -> x == 5
    opposite_ops = {
        "==": "==",
        "!=": "!=",
        "<": ">=",
        ">=": "<",
        ">": "<=",
        "<=": ">",
    }

    props = {"x": "positions", "y": "positions", "z": "positions"}

    def __init__(self, parser, tokens):
        """
        Possible splitting around operator:

        prop x < 5
        prop x< 5
        prop x <5
        prop x<5
        """
        super().__init__(parser, tokens)

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
            errmsg = (
                f"Invalid operator : '{oper}' Use one of : "
                f"'{self.ops.keys()}'"
            )
            raise ValueError(errmsg) from None
        else:
            if oper == "==":
                self.operator = functools.partial(
                    self.operator, atol=parser.atol, rtol=parser.rtol
                )
        self.value = float(value)

    def _apply(self, group):
        try:
            values = getattr(group, self.props[self.prop])
        except KeyError:
            errmsg = f"Expected one of {list(self.props.keys())}"
            raise SelectionError(errmsg) from None
        except NoDataError:
            attr = self.props[self.prop]
            errmsg = f"This Universe does not contain {attr} information"
            raise SelectionError(errmsg) from None

        try:
            col = {"x": 0, "y": 1, "z": 2}[self.prop]
        except KeyError:
            pass
        else:
            values = values[:, col]

        if self.absolute:
            values = np.abs(values)
        mask = self.operator(values, self.value)

        return group[mask]


class SameSelection(Selection):
    """
    Selects all atoms that have the same subkeyword value as any atom in selection

    .. versionchanged:: 1.0.0
       Map :code:`"residue"` to :code:`"resindices"` and :code:`"segment"` to
       :code:`"segindices"` (see #2669 and #2672)
    """

    token = "same"
    precedence = 1

    prop_trans = {
        "fragment": None,
        "x": None,
        "y": None,
        "z": None,
        "residue": "resindices",
        "segment": "segindices",
        "name": "names",
        "type": "types",
        "resname": "resnames",
        "resid": "resids",
        "segid": "segids",
        "mass": "masses",
        "charge": "charges",
        "radius": "radii",
        "bfactor": "bfactors",
        "resnum": "resnums",
    }

    def __init__(self, parser, tokens):
        super().__init__(parser, tokens)

        prop = tokens.popleft()
        if prop not in self.prop_trans:
            raise ValueError(
                "Unknown same property : {0}"
                "Choose one of : {1}"
                "".format(prop, self.prop_trans.keys())
            )
        self.prop = prop
        parser.expect("as")
        self.sel = parser.parse_expression(self.precedence)
        self.prop = prop

    def _apply(self, group):
        res = self.sel.apply(group)
        if not res:
            return group[[]]  # empty selection

        # Fragment must come before self.prop_trans lookups!
        if self.prop == "fragment":
            # Combine all fragments together, then check where group
            # indices are same as fragment(s) indices
            allfrags = functools.reduce(lambda x, y: x + y, res.fragments)

            mask = np.isin(group.indices, allfrags.indices)
            return group[mask]
        # [xyz] must come before self.prop_trans lookups too!
        try:
            pos_idx = {"x": 0, "y": 1, "z": 2}[self.prop]
        except KeyError:
            # The self.prop string was already checked,
            # so don't need error checking here.
            # KeyError at this point is impossible!
            attrname = self.prop_trans[self.prop]
            vals = getattr(res, attrname)
            mask = np.isin(getattr(group, attrname), vals)

            return group[mask]
        else:
            vals = res.positions[:, pos_idx]
            pos = group.positions[:, pos_idx]

            # isclose only does one value at a time
            mask = np.vstack([np.isclose(pos, v) for v in vals]).any(axis=0)
            return group[mask]


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
                "".format(self.tokens[0], token)
            )

    def parse(
        self,
        selectstr,
        selgroups,
        periodic=None,
        atol=1e-08,
        rtol=1e-05,
        sorted=True,
        rdkit_kwargs=None,
        smarts_kwargs=None,
    ):
        """Create a Selection object from a string.

        Parameters
        ----------
        selectstr : str
            The string that describes the selection
        selgroups : AtomGroups
            AtomGroups to be used in `group` selections
        periodic : bool, optional
            for distance based selections, whether to consider
            periodic boundary conditions
        atol : float, optional
            The absolute tolerance parameter for float comparisons.
            Passed to :func:`numpy.isclose`.
        rtol : float, optional
            The relative tolerance parameter for float comparisons.
            Passed to :func:`numpy.isclose`.
        sorted : bool, optional
            Whether to sorted the output AtomGroup.
        rdkit_kwargs : dict, optional
            Arguments passed to the RDKitConverter when using selection based
            on SMARTS queries
        smarts_kwargs : dict, optional
          Arguments passed internally to RDKit's `GetSubstructMatches
          <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol.GetSubstructMatches>`_.

        Returns
        -------
        The appropriate Selection object.  Use the .apply method on
        this to perform the selection.

        Raises
        ------
        SelectionError
            If anything goes wrong in creating the Selection object.


        .. versionchanged:: 2.0.0
            Added `atol` and `rtol` keywords to select float values. Added
            `rdkit_kwargs` to pass arguments to the RDKitConverter
        .. versionchanged:: 2.2.0
            Added ``smarts_kwargs`` argument, allowing users to pass a
            a dictionary of arguments to RDKit's ``GetSubstructMatches``.
        """
        self.periodic = periodic
        self.atol = atol
        self.rtol = rtol
        self.sorted = sorted
        self.rdkit_kwargs = rdkit_kwargs or {}
        self.smarts_kwargs = smarts_kwargs or {}

        self.selectstr = selectstr
        self.selgroups = selgroups
        tokens = selectstr.replace("(", " ( ").replace(")", " ) ")
        self.tokens = collections.deque(tokens.split() + [None])
        parsetree = self.parse_expression(0)
        if self.tokens[0] is not None:
            raise SelectionError(
                "Unexpected token at end of selection string: '{0}'"
                "".format(self.tokens[0])
            )
        return parsetree

    def parse_expression(self, p):
        exp1 = self._parse_subexp()
        while (
            self.tokens[0] in _OPERATIONS
            and _OPERATIONS[self.tokens[0]].precedence >= p
        ):
            op = _OPERATIONS[self.tokens.popleft()]
            q = 1 + op.precedence
            exp2 = self.parse_expression(q)
            exp1 = op(exp1, exp2, self)
        return exp1

    def _parse_subexp(self):
        op = self.tokens.popleft()

        if op == "(":
            exp = self.parse_expression(0)
            self.expect(")")
            return exp

        try:
            return _SELECTIONDICT[op](self, self.tokens)
        except KeyError:
            errmsg = f"Unknown selection token: '{op}'"
            raise SelectionError(errmsg) from None
        except ValueError as e:
            errmsg = f"Selection failed: '{e}'"
            raise SelectionError(errmsg) from None


# The module level instance
Parser = SelectionParser()

# create a module container to avoid name clashes of autogenerated classes
_selectors = types.ModuleType(
    f"{__name__}._selectors", doc="Automatically generated selectors"
)
# stick it in sys.modules so pickle can find it
sys.modules[_selectors.__name__] = _selectors


def gen_selection_class(singular, attrname, dtype, per_object):
    """Selection class factory for arbitrary TopologyAttrs.

    Normally this should not be used except within the codebase
    or by developers; it is called by the metaclass
    :class:`MDAnalysis.core.topologyattrs._TopologyAttrMeta` to
    auto-generate suitable selection classes by creating a token
    with the topology attribute (singular) name. The function
    uses the provided ``dtype`` to choose which Selection class
    to subclass:

    * :class:`BoolSelection` for booleans
    * :class:`RangeSelection` for integers
    * :class:`FloatRangeSelection` for floats
    * :class:`_ProtoStringSelection` for strings

    Other value types are not yet supported and will raise a
    ValueError. The classes are created in the :mod:`_selectors`
    module to avoid namespace clashes.

    Parameters
    ----------
    singular: str
        singular name of TopologyAttr
    attrname: str
        attribute name of TopologyAttr
    dtype: type
        type of TopologyAttr
    per_object: str
        level of TopologyAttr

    Returns
    -------
    selection: subclass of Selection

    Raises
    ------
    ValueError
        If ``dtype`` is not one of the supported types


    Example
    -------

    The function creates a class inside ``_selectors`` and returns it.
    Normally it should not need to be manually called, as it is created
    for each TopologyAttr::

        >>> gen_selection_class("resname", "resnames", object, "residue")
        <class 'MDAnalysis.core.selection._selectors.ResnameSelection'>

    Simply generating this selector is sufficient for the keyword to be
    accessible by :meth:`~MDAnalysis.core.universe.Universe.select_atoms`,
    as that is automatically handled by
    :class:`~MDAnalysis.core.selections._Selectionmeta`.

    See also
    --------
    :class:`MDAnalysis.core.topologyattrs._TopologyAttrMeta`

    .. versionadded:: 2.0.0
    """
    basedct = {
        "token": singular,
        "field": attrname,
        # manually make modules the _selectors wrapper
        "__module__": _selectors.__name__,
    }
    name = f"{singular.capitalize()}Selection"

    if dtype == "U1":  # order is important here, U1 will trip up issubclass
        basecls = SingleCharSelection
    elif issubclass(dtype, bool):
        basecls = BoolSelection
    elif np.issubdtype(dtype, np.integer):
        basecls = RangeSelection
    elif np.issubdtype(dtype, np.floating):
        basecls = FloatRangeSelection
    elif issubclass(dtype, str) or dtype == object:
        basecls = _ProtoStringSelection
        if per_object == "segment":
            basedct["level"] = "segindices"
        elif per_object == "residue":
            basedct["level"] = "resindices"
        else:
            basedct["level"] = "ix"
    else:
        raise ValueError(
            f"No base class defined for dtype {dtype}. "
            "Define a Selection class manually by "
            "subclassing core.selection.Selection"
        )

    cls = type(name, (basecls,), basedct)
    setattr(_selectors, name, cls)  # stick it in _selectors
    return cls
