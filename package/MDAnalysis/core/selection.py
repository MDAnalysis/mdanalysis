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


class SugarSelection(Selection):
    """Consists of sugar residues with recognized abbreviations.

    Recognized residue names in :attr:`SugarSelection.sugar_res`.

      * from glycam.org server::
         Abbreviations for PDB, CHARMM and GLYCAM
         https://glycam.org/docs/othertoolsservice/2016/06/09/3d-snfg-list-of-residue-names/index.html

      * manually added aglycans 'ROH', 'OME', 'TBT'
        from GLYCAM-Web generated files.

    .. versionadded:: 2.11.0
    """

    token = "sugar"

    sugar_res = {
        # https://glycam.org/docs/othertoolsservice/2016/06/09/3d-snfg-list-of-residue-names/index.html
        # Full PDB Abbreviations Nomenclature
        "GLC",
        "MAL",
        "BGC",
        "NAG",
        "4YS",
        "SGN",
        "BGLN",
        "NDG",
        "GCS",
        "GCU",
        "QUI",
        "OLI",
        "MAN",
        "BMA",
        "MAV",
        "BEM",
        "RAM",
        "TYV",
        "ARA",
        "AHR",
        "GAL",
        "GLA",
        "NGA",
        "ADA",
        "GUL",
        "GUP",
        "GL0",
        "LGU",
        "ALT",
        "ALL",
        "WOO",
        "TAL",
        "IDO",
        "IDS",
        "FUC",
        "FUL",
        "LYX",
        "ABE",
        "XYL",
        "XYS",
        "LXC",
        "XYP",
        "PAR",
        "RIB",
        "DIG",
        "COL",
        "BAC",
        "API",
        "FRU",
        "TAG",
        "SOR",
        "PSI",
        "DHA",
        "KDN",
        "KDO",
        "NEU",
        "SIA",
        "MUR",
        "GMH",
        # Full CHARMM Abbreviations Nomenclature
        "AGLC",
        "BGLC",
        "AGLCNA",
        "BGLCNA",
        "BGLCN0",
        "AGLCA",
        "BGLCA",
        "BGLCA0",
        "AMAN",
        "BMAN",
        "ARHM",
        "BRHM",
        "AARB",
        "BARB",
        "AGAL",
        "BGAL",
        "AGALNA",
        "BGALNA",
        "AGUL",
        "BGUL",
        "AALT",
        "BALT",
        "AALL",
        "ATAL",
        "BTAL",
        "AIDO",
        "BIDO",
        "AIDOA",
        "BIDOA",
        "AFUC",
        "BFUC",
        "ALYF",
        "BLYF",
        "AXYL",
        "BXYL",
        "AXYF",
        "BXYF",
        "ARIB",
        "BRIB",
        "AFRU",
        "BFRU",
        "ANE5AC",
        "BNE5AC",
        # GLYCAM Abbreviations
        # Glucose Nomenculature
        "0GA",
        "0GB",
        "1GA",
        "1GB",
        "2GA",
        "2GB",
        "3GA",
        "3GB",
        "4GA",
        "4GB",
        "6GA",
        "6GB",
        "ZGA",
        "ZGB",
        "YGA",
        "YGB",
        "XGA",
        "XGB",
        "WGA",
        "WGB",
        "VGA",
        "VGB",
        "UGA",
        "UGB",
        "TGA",
        "TGB",
        "SGA",
        "SGB",
        "RGA",
        "RGB",
        "QGA",
        "QGB",
        "PGA",
        "PGB",
        "0gA",
        "0gB",
        "1gA",
        "1gB",
        "2gA",
        "2gB",
        "3gA",
        "3gB",
        "4gA",
        "4gB",
        "6gA",
        "6gB",
        "ZgA",
        "ZgB",
        "YgA",
        "YgB",
        "XgA",
        "XgB",
        "WgA",
        "WgB",
        "VgA",
        "VgB",
        "UgA",
        "UgB",
        "TgA",
        "TgB",
        "SgA",
        "SgB",
        "RgA",
        "RgB",
        "QgA",
        "QgB",
        "PgA",
        "PgB",
        # N-Acetyl Glucosamine Nomenculature
        "0YA",
        "0YB",
        "1YA",
        "1YB",
        "3YA",
        "3YB",
        "4YA",
        "4YB",
        "6YA",
        "6YB",
        "WYA",
        "WYB",
        "VYA",
        "VYB",
        "UYA",
        "UYB",
        "QYA",
        "QYB",
        "0yA",
        "0yB",
        "1yA",
        "1yB",
        "3yA",
        "3yB",
        "4yA",
        "4yB",
        "6yA",
        "6yB",
        "WyA",
        "WyB",
        "VyA",
        "VyB",
        "UyA",
        "UyB",
        "QyA",
        "QyB",
        # Glucosamine Nomenculature
        "0YN",
        "0Yn",
        "0YNP",
        "0YnP",
        "0YS",
        "0Ys",
        "3YS",
        "3Ys",
        "4YS",
        "4Ys",
        "6YS",
        "6Ys",
        "QYS",
        "QYs",
        "UYS",
        "UYs",
        "VYS",
        "VYs",
        "WYS",
        "WYs",
        "0yS",
        "0ys",
        "3yS",
        "3ys",
        "4yS",
        "4ys",
        # Glucuronic Acid Nomenculature
        "0ZA",
        "0ZB",
        "1ZA",
        "1ZB",
        "2ZA",
        "2ZB",
        "3ZA",
        "3ZB",
        "4ZA",
        "4ZB",
        "ZZA",
        "ZZB",
        "YZA",
        "YZB",
        "WZA",
        "WZB",
        "TZA",
        "TZB",
        "0zA",
        "0zB",
        "1zA",
        "1zB",
        "2zA",
        "2zB",
        "3zA",
        "3zB",
        "4zA",
        "4zB",
        "ZzA",
        "ZzB",
        "YzA",
        "YzB",
        "WzA",
        "WzB",
        "TzA",
        "TzB",
        "0ZBP",
        # Quinovose Nomenculature
        "0QA",
        "0QB",
        "1QA",
        "1QB",
        "2QA",
        "2QB",
        "3QA",
        "3QB",
        "4QA",
        "4QB",
        "ZQA",
        "ZQB",
        "YQA",
        "YQB",
        "WQA",
        "WQB",
        "TQA",
        "TQB",
        "0qA",
        "0qB",
        "1qA",
        "1qB",
        "2qA",
        "2qB",
        "3qA",
        "3qB",
        "4qA",
        "4qB",
        "ZqA",
        "ZqB",
        "YqA",
        "YqB",
        "WqA",
        "WqB",
        "TqA",
        "TqB",
        # Mannose Nomenculature
        "0MA",
        "0MB",
        "1MA",
        "1MB",
        "2MA",
        "2MB",
        "3MA",
        "3MB",
        "4MA",
        "4MB",
        "6MA",
        "6MB",
        "ZMA",
        "ZMB",
        "YMA",
        "YMB",
        "XMA",
        "XMB",
        "WMA",
        "WMB",
        "VMA",
        "VMB",
        "UMA",
        "UMB",
        "TMA",
        "TMB",
        "SMA",
        "SMB",
        "RMA",
        "RMB",
        "QMA",
        "QMB",
        "PMA",
        "PMB",
        "0mA",
        "0mB",
        "1mA",
        "1mB",
        "2mA",
        "2mB",
        "3mA",
        "3mB",
        "4mA",
        "4mB",
        "6mA",
        "6mB",
        "ZmA",
        "ZmB",
        "YmA",
        "YmB",
        "XmA",
        "XmB",
        "WmA",
        "WmB",
        "VmA",
        "VmB",
        "UmA",
        "UmB",
        "TmA",
        "TmB",
        "SmA",
        "SmB",
        "RmA",
        "RmB",
        "QmA",
        "QmB",
        "PmA",
        "PmB",
        # N-Acetyl Mannosamine Nomenculature
        "0WA",
        "0WB",
        "1WA",
        "1WB",
        "3WA",
        "3WB",
        "4WA",
        "4WB",
        "6WA",
        "6WB",
        "WWA",
        "WWB",
        "VWA",
        "VWB",
        "UWA",
        "UWB",
        "QWA",
        "QWB",
        "0wA",
        "0wB",
        "1wA",
        "1wB",
        "3wA",
        "3wB",
        "4wA",
        "4wB",
        "6wA",
        "6wB",
        "WwA",
        "WwB",
        "VwA",
        "VwB",
        "UwA",
        "UwB",
        "QwA",
        "QwB",
        # Rhamnose Nomenculature
        "0HA",
        "0HB",
        "1HA",
        "1HB",
        "2HA",
        "2HB",
        "3HA",
        "3HB",
        "4HA",
        "4HB",
        "ZHA",
        "ZHB",
        "YHA",
        "YHB",
        "WHA",
        "WHB",
        "THA",
        "THB",
        "0hA",
        "0hB",
        "1hA",
        "1hB",
        "2hA",
        "2hB",
        "3hA",
        "3hB",
        "4hA",
        "4hB",
        "ZhA",
        "ZhB",
        "YhA",
        "YhB",
        "WhA",
        "WhB",
        "ThA",
        "ThB",
        # Tyvelose Nomenculature
        "0TV",
        "0Tv",
        "1TV",
        "1Tv",
        "2TV",
        "2Tv",
        "4TV",
        "4Tv",
        "YTV",
        "YTv",
        "0tV",
        "0tv",
        "1tV",
        "1tv",
        "2tV",
        "2tv",
        "4tV",
        "4tv",
        "YtV",
        "Ytv",
        # Arabinose Nomenculature
        "0AA",
        "0AB",
        "1AA",
        "1AB",
        "2AA",
        "2AB",
        "3AA",
        "3AB",
        "4AA",
        "4AB",
        "ZAA",
        "ZAB",
        "YAA",
        "YAB",
        "WAA",
        "WAB",
        "TAA",
        "TAB",
        "0AD",
        "0AU",
        "1AD",
        "1AU",
        "2AD",
        "2AU",
        "3AD",
        "3AU",
        "5AD",
        "5AU",
        "ZAD",
        "ZAU",
        "0aA",
        "0aB",
        "1aA",
        "1aB",
        "2aA",
        "2aB",
        "3aA",
        "3aB",
        "4aA",
        "4aB",
        "ZaA",
        "ZaB",
        "YaA",
        "YaB",
        "WaA",
        "WaB",
        "TaA",
        "TaB",
        "0aD",
        "0aU",
        "1aD",
        "1aU",
        "2aD",
        "2aU",
        "3aD",
        "3aU",
        "5aD",
        "5aU",
        "ZaD",
        "ZaU",
        # Galactose Nomenculature
        "0LA",
        "0LB",
        "1LA",
        "1LB",
        "2LA",
        "2LB",
        "3LA",
        "3LB",
        "4LA",
        "4LB",
        "6LA",
        "6LB",
        "ZLA",
        "ZLB",
        "YLA",
        "YLB",
        "XLA",
        "XLB",
        "WLA",
        "WLB",
        "VLA",
        "VLB",
        "ULA",
        "ULB",
        "TLA",
        "TLB",
        "SLA",
        "SLB",
        "RLA",
        "RLB",
        "QLA",
        "QLB",
        "PLA",
        "PLB",
        "0lA",
        "0lB",
        "1lA",
        "1lB",
        "2lA",
        "2lB",
        "3lA",
        "3lB",
        "4lA",
        "4lB",
        "6lA",
        "6lB",
        "ZlA",
        "ZlB",
        "YlA",
        "YlB",
        "XlA",
        "XlB",
        "WlA",
        "WlB",
        "VlA",
        "VlB",
        "UlA",
        "UlB",
        "TlA",
        "TlB",
        "SlA",
        "SlB",
        "RlA",
        "RlB",
        "QlA",
        "QlB",
        "PlA",
        "PlB",
        # N-Acetyl Galactosamine Nomenculature
        "0VA",
        "0VB",
        "1VA",
        "1VB",
        "3VA",
        "3VB",
        "4VA",
        "4VB",
        "6VA",
        "6VB",
        "WVA",
        "WVB",
        "VVA",
        "VVB",
        "UVA",
        "UVB",
        "QVA",
        "QVB",
        "0vA",
        "0vB",
        "1vA",
        "1vB",
        "3vA",
        "3vB",
        "4vA",
        "4vB",
        "6vA",
        "6vB",
        "WvA",
        "WvB",
        "VvA",
        "VvB",
        "UvA",
        "UvB",
        "QvA",
        "QvB",
        # Galacturonic Acid Nomenculature
        "0OA",
        "0OB",
        "1OA",
        "1OB",
        "2OA",
        "2OB",
        "3OA",
        "3OB",
        "4OA",
        "4OB",
        "ZOA",
        "ZOB",
        "YOA",
        "YOB",
        "WOA",
        "WOB",
        "TOA",
        "TOB",
        "0oA",
        "0oB",
        "1oA",
        "1oB",
        "2oA",
        "2oB",
        "3oA",
        "3oB",
        "4oA",
        "4oB",
        "ZoA",
        "ZoB",
        "YoA",
        "YoB",
        "WoA",
        "WoB",
        "ToA",
        "ToB",
        # Gulose Nomenculature
        "0KA",
        "0KB",
        "1KA",
        "1KB",
        "2KA",
        "2KB",
        "3KA",
        "3KB",
        "4KA",
        "4KB",
        "6KA",
        "6KB",
        "ZKA",
        "ZKB",
        "YKA",
        "YKB",
        "XKA",
        "XKB",
        "WKA",
        "WKB",
        "VKA",
        "VKB",
        "UKA",
        "UKB",
        "TKA",
        "TKB",
        "SKA",
        "SKB",
        "RKA",
        "RKB",
        "QKA",
        "QKB",
        "PKA",
        "PKB",
        "0kA",
        "0kB",
        "1kA",
        "1kB",
        "2kA",
        "2kB",
        "3kA",
        "3kB",
        "4kA",
        "4kB",
        "6kA",
        "6kB",
        "ZkA",
        "ZkB",
        "YkA",
        "YkB",
        "XkA",
        "XkB",
        "WkA",
        "WkB",
        "VkA",
        "VkB",
        "UkA",
        "UkB",
        "TkA",
        "TkB",
        "SkA",
        "SkB",
        "RkA",
        "RkB",
        "QkA",
        "QkB",
        "PkA",
        "PkB",
        # Altrose Nomenculature
        "0EA",
        "0EB",
        "1EA",
        "1EB",
        "2EA",
        "2EB",
        "3EA",
        "3EB",
        "4EA",
        "4EB",
        "6EA",
        "6EB",
        "ZEA",
        "ZEB",
        "YEA",
        "YEB",
        "XEA",
        "XEB",
        "WEA",
        "WEB",
        "VEA",
        "VEB",
        "UEA",
        "UEB",
        "TEA",
        "TEB",
        "SEA",
        "SEB",
        "REA",
        "REB",
        "QEA",
        "QEB",
        "PEA",
        "PEB",
        "0eA",
        "0eB",
        "1eA",
        "1eB",
        "2eA",
        "2eB",
        "3eA",
        "3eB",
        "4eA",
        "4eB",
        "6eA",
        "6eB",
        "ZeA",
        "ZeB",
        "YeA",
        "YeB",
        "XeA",
        "XeB",
        "WeA",
        "WeB",
        "VeA",
        "VeB",
        "UeA",
        "UeB",
        "TeA",
        "TeB",
        "SeA",
        "SeB",
        "ReA",
        "ReB",
        "QeA",
        "QeB",
        "PeA",
        "PeB",
        # Allose Nomenculature
        "0NA",
        "0NB",
        "1NA",
        "1NB",
        "2NA",
        "2NB",
        "3NA",
        "3NB",
        "4NA",
        "4NB",
        "6NA",
        "6NB",
        "ZNA",
        "ZNB",
        "YNA",
        "YNB",
        "XNA",
        "XNB",
        "WNA",
        "WNB",
        "VNA",
        "VNB",
        "UNA",
        "UNB",
        "TNA",
        "TNB",
        "SNA",
        "SNB",
        "RNA",
        "RNB",
        "QNA",
        "QNB",
        "PNA",
        "PNB",
        "0nA",
        "0nB",
        "1nA",
        "1nB",
        "2nA",
        "2nB",
        "3nA",
        "3nB",
        "4nA",
        "4nB",
        "6nA",
        "6nB",
        "ZnA",
        "ZnB",
        "YnA",
        "YnB",
        "XnA",
        "XnB",
        "WnA",
        "WnB",
        "VnA",
        "VnB",
        "UnA",
        "UnB",
        "TnA",
        "TnB",
        "SnA",
        "SnB",
        "RnA",
        "RnB",
        "QnA",
        "QnB",
        "PnA",
        "PnB",
        # Talose Nomenculature
        "0TA",
        "0TB",
        "1TA",
        "1TB",
        "2TA",
        "2TB",
        "3TA",
        "3TB",
        "4TA",
        "4TB",
        "6TA",
        "6TB",
        "ZTA",
        "ZTB",
        "YTA",
        "YTB",
        "XTA",
        "XTB",
        "WTA",
        "WTB",
        "VTA",
        "VTB",
        "UTA",
        "UTB",
        "TTA",
        "TTB",
        "STA",
        "STB",
        "RTA",
        "RTB",
        "QTA",
        "QTB",
        "PTA",
        "PTB",
        "0tA",
        "0tB",
        "1tA",
        "1tB",
        "2tA",
        "2tB",
        "3tA",
        "3tB",
        "4tA",
        "4tB",
        "6tA",
        "6tB",
        "ZtA",
        "ZtB",
        "YtA",
        "YtB",
        "XtA",
        "XtB",
        "WtA",
        "WtB",
        "VtA",
        "VtB",
        "UtA",
        "UtB",
        "TtA",
        "TtB",
        "StA",
        "StB",
        "RtA",
        "RtB",
        "QtA",
        "QtB",
        "PtA",
        "PtB",
        # Iduronic Acid Nomenculature
        "0UA",
        "0UB",
        "1UA",
        "1UB",
        "2UA",
        "2UB",
        "3UA",
        "3UB",
        "4UA",
        "4UB",
        "ZUA",
        "ZUB",
        "YUA",
        "YUB",
        "WUA",
        "WUB",
        "TUA",
        "TUB",
        "0uA",
        "0uB",
        "1uA",
        "1uB",
        "2uA",
        "2uB",
        "3uA",
        "3uB",
        "4uA",
        "4uB",
        "ZuA",
        "ZuB",
        "YuA",
        "YuB",
        "WuA",
        "WuB",
        "TuA",
        "TuB",
        "YuAP",
        # Fucose Nomenculature
        "0FA",
        "0FB",
        "1FA",
        "1FB",
        "2FA",
        "2FB",
        "3FA",
        "3FB",
        "4FA",
        "4FB",
        "ZFA",
        "ZFB",
        "YFA",
        "YFB",
        "WFA",
        "WFB",
        "TFA",
        "TFB",
        "0fA",
        "0fB",
        "1fA",
        "1fB",
        "2fA",
        "2fB",
        "3fA",
        "3fB",
        "4fA",
        "4fB",
        "ZfA",
        "ZfB",
        "YfA",
        "YfB",
        "WfA",
        "WfB",
        "TfA",
        "TfB",
        # Lyxose Nomenculature
        "0DA",
        "0DB",
        "1DA",
        "1DB",
        "2DA",
        "2DB",
        "3DA",
        "3DB",
        "4DA",
        "4DB",
        "ZDA",
        "ZDB",
        "YDA",
        "YDB",
        "WDA",
        "WDB",
        "TDA",
        "TDB",
        "0DD",
        "0DU",
        "1DD",
        "1DU",
        "2DD",
        "2DU",
        "3DD",
        "3DU",
        "5DD",
        "5DU",
        "ZDD",
        "ZDU",
        "0dA",
        "0dB",
        "1dA",
        "1dB",
        "2dA",
        "2dB",
        "3dA",
        "3dB",
        "4dA",
        "4dB",
        "ZdA",
        "ZdB",
        "YdA",
        "YdB",
        "WdA",
        "WdB",
        "TdA",
        "TdB",
        "0dD",
        "0dU",
        "1dD",
        "1dU",
        "2dD",
        "2dU",
        "3dD",
        "3dU",
        "5dD",
        "5dU",
        "ZdD",
        "ZdU",
        # Abequose Nomenculature
        "0AE",
        "2AE",
        "4AE",
        "YGa",
        "0AF",
        "2AF",
        "4AF",
        "YAF",
        # Xylose Nomenculature
        "0XA",
        "0XB",
        "1XA",
        "1XB",
        "2XA",
        "2XB",
        "3XA",
        "3XB",
        "4XA",
        "4XB",
        "ZXA",
        "ZXB",
        "YXA",
        "YXB",
        "WXA",
        "WXB",
        "TXA",
        "TXB",
        "0XD",
        "0XU",
        "1XD",
        "1XU",
        "2XD",
        "2XU",
        "3XD",
        "3XU",
        "5XD",
        "5XU",
        "ZXD",
        "ZXU",
        "0xA",
        "0xB",
        "1xA",
        "1xB",
        "2xA",
        "2xB",
        "3xA",
        "3xB",
        "4xA",
        "4xB",
        "ZxA",
        "ZxB",
        "YxA",
        "YxB",
        "WxA",
        "WxB",
        "TxA",
        "TxB",
        "0xD",
        "0xU",
        "1xD",
        "1xU",
        "2xD",
        "2xU",
        "3xD",
        "3xU",
        "5xD",
        "5xU",
        "ZxD",
        "ZxU",
        # Ribose Nomenculature
        "0RA",
        "0RB",
        "1RA",
        "1RB",
        "2RA",
        "2RB",
        "3RA",
        "3RB",
        "4RA",
        "4RB",
        "ZRA",
        "ZRB",
        "YRA",
        "YRB",
        "WRA",
        "WRB",
        "TRA",
        "TRB",
        "0RD",
        "0RU",
        "1RD",
        "1RU",
        "2RD",
        "2RU",
        "3RD",
        "3RU",
        "5RD",
        "5RU",
        "ZRD",
        "ZRU",
        "0rA",
        "0rB",
        "1rA",
        "1rB",
        "2rA",
        "2rB",
        "3rA",
        "3rB",
        "4rA",
        "4rB",
        "ZrA",
        "ZrB",
        "YrA",
        "YrB",
        "WrA",
        "WrB",
        "TrA",
        "TrB",
        "0rD",
        "0rU",
        "1rD",
        "1rU",
        "2rD",
        "2rU",
        "3rD",
        "3rU",
        "5rD",
        "5rU",
        "ZrD",
        "ZrU",
        # Bacillosamine Nomenculature
        "0BC",
        "3BC",
        "0bC",
        "3bC",
        # Fructose Nomenculature
        "0CA",
        "0CB",
        "1CA",
        "1CB",
        "2CA",
        "2CB",
        "3CA",
        "3CB",
        "4CA",
        "4CB",
        "5CA",
        "5CB",
        "WCA",
        "WCB",
        "0CD",
        "0CU",
        "1CD",
        "1CU",
        "2CD",
        "2CU",
        "3CD",
        "3CU",
        "4CD",
        "4CU",
        "6CD",
        "6CU",
        "WCD",
        "WCU",
        "VCD",
        "VCU",
        "UCD",
        "UCU",
        "QCD",
        "QCU",
        "0cA",
        "0cB",
        "1cA",
        "1cB",
        "2cA",
        "2cB",
        "3cA",
        "3cB",
        "4cA",
        "4cB",
        "5cA",
        "5cB",
        "WcA",
        "WcB",
        "0cD",
        "0cU",
        "1cD",
        "1cU",
        "2cD",
        "2cU",
        "3cD",
        "3cU",
        "4cD",
        "4cU",
        "6cD",
        "6cU",
        "WcD",
        "WcU",
        "VcD",
        "VcU",
        "UcD",
        "UcU",
        "QcD",
        "QcU",
        # Tagatose Nomenculature
        "0JA",
        "0JB",
        "1JA",
        "1JB",
        "2JA",
        "2JB",
        "3JA",
        "3JB",
        "4JA",
        "4JB",
        "5JA",
        "5JB",
        "WJA",
        "WJB",
        "0JD",
        "0JU",
        "1JD",
        "1JU",
        "2JD",
        "2JU",
        "3JD",
        "3JU",
        "4JD",
        "4JU",
        "6JD",
        "6JU",
        "WJD",
        "WJU",
        "VJD",
        "VJU",
        "UJD",
        "UJU",
        "QJD",
        "QJU",
        "0jA",
        "0jB",
        "1jA",
        "1jB",
        "2jA",
        "2jB",
        "3jA",
        "3jB",
        "4jA",
        "4jB",
        "5jA",
        "5jB",
        "WjA",
        "WjB",
        "0jD",
        "0jU",
        "1jD",
        "1jU",
        "2jD",
        "2jU",
        "3jD",
        "3jU",
        "4jD",
        "4jU",
        "6jD",
        "6jU",
        "WjD",
        "WjU",
        "VjD",
        "VjU",
        "UjD",
        "UjU",
        "QjD",
        "QjU",
        # Sorbose Nomenculature
        "0BA",
        "0BB",
        "1BA",
        "1BB",
        "2BA",
        "2BB",
        "3BA",
        "3BB",
        "4BA",
        "4BB",
        "5BA",
        "5BB",
        "WBA",
        "WBB",
        "0BD",
        "0BU",
        "1BD",
        "1BU",
        "2BD",
        "2BU",
        "3BD",
        "3BU",
        "4BD",
        "4BU",
        "6BD",
        "6BU",
        "WBD",
        "WBU",
        "VBD",
        "VBU",
        "UBD",
        "UBU",
        "QBD",
        "QBU",
        "0bA",
        "0bB",
        "1bA",
        "1bB",
        "2bA",
        "2bB",
        "3bA",
        "3bB",
        "4bA",
        "4bB",
        "5bA",
        "5bB",
        "WbA",
        "WbB",
        "0bD",
        "0bU",
        "1bD",
        "1bU",
        "2bD",
        "2bU",
        "3bD",
        "3bU",
        "4bD",
        "4bU",
        "6bD",
        "6bU",
        "WbD",
        "WbU",
        "VbD",
        "VbU",
        "UbD",
        "UbU",
        "QbD",
        "QbU",
        # Psicose Nomenculature
        "0PA",
        "0PB",
        "1PA",
        "1PB",
        "2PA",
        "2PB",
        "3PA",
        "3PB",
        "4PA",
        "4PB",
        "5PA",
        "5PB",
        "WPA",
        "WPB",
        "0PD",
        "0PU",
        "1PD",
        "1PU",
        "2PD",
        "2PU",
        "3PD",
        "3PU",
        "4PD",
        "4PU",
        "6PD",
        "6PU",
        "WPD",
        "WPU",
        "VPD",
        "VPU",
        "UPD",
        "UPU",
        "QPD",
        "QPU",
        "0pA",
        "0pB",
        "1pA",
        "1pB",
        "2pA",
        "2pB",
        "3pA",
        "3pB",
        "4pA",
        "4pB",
        "5pA",
        "5pB",
        "WpA",
        "WpB",
        "0pD",
        "0pU",
        "1pD",
        "1pU",
        "2pD",
        "2pU",
        "3pD",
        "3pU",
        "4pD",
        "4pU",
        "6pD",
        "6pU",
        "WpD",
        "WpU",
        "VpD",
        "VpU",
        "UpD",
        "UpU",
        "QpD",
        "QpU",
        # N-Acetyl Neuraminic Acid Nomenculature
        "0SA",
        "0SB",
        "4SA",
        "4SB",
        "7SA",
        "7SB",
        "8SA",
        "8SB",
        "9SA",
        "9SB",
        "ASA",
        "ASB",
        "BSA",
        "BSB",
        "CSA",
        "CSB",
        "DSA",
        "DSB",
        "ESA",
        "ESB",
        "FSA",
        "FSB",
        "GSA",
        "GSB",
        "HSA",
        "HSB",
        "ISA",
        "ISB",
        "JSA",
        "JSB",
        "KSA",
        "KSB",
        "0sA",
        "0sB",
        "4sA",
        "4sB",
        "7sA",
        "7sB",
        "8sA",
        "8sB",
        "9sA",
        "9sB",
        "AsA",
        "AsB",
        "BsA",
        "BsB",
        "CsA",
        "CsB",
        "DsA",
        "DsB",
        "EsA",
        "EsB",
        "FsA",
        "FsB",
        "GsA",
        "GsB",
        "HsA",
        "HsB",
        "IsA",
        "IsB",
        "JsA",
        "JsB",
        "KsA",
        "KsB",
        # N-Glycolyl Neuraminic Acid Nomenculature
        "0GL",
        "4GL",
        "7GL",
        "8GL",
        "9GL",
        "CGL",
        "DGL",
        "EGL",
        "FGL",
        "GGL",
        "HGL",
        "IGL",
        "JGL",
        "KGL",
        "0gL",
        "4gL",
        "7gL",
        "8gL",
        "9gL",
        "AgL",
        "BgL",
        "CgL",
        "DgL",
        "EgL",
        "FgL",
        "GgL",
        "HgL",
        "IgL",
        "JgL",
        "KgL",
        # Aglycon Nomenculature
        "ROH",
        "OME",
        "TBT",
    }

    def _apply(self, group):
        resname_attr = group.universe._topology.resnames
        # which values in resname attr are in sugar_res?
        matches = [
            ix
            for (nm, ix) in resname_attr.namedict.items()
            if nm in self.sugar_res
        ]
        # index of each atom's resname
        nmidx = resname_attr.nmidx[group.resindices]
        # intersect atom's resname index and matches to sugar_res
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
