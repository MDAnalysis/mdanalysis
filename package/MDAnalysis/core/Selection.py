# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Atom selection Hierarchy --- :mod:`MDAnalysis.core.Selection`
=============================================================

These objects are constructed and applied to the group

In general, Parser.parse() creates a Selection object
from a selection string.

This Selection object is then passed an AtomGroup through its
apply method to apply the Selection to the AtomGroup.

This is all invisible to the user through ag.select_atoms
"""
import six
from six.moves import zip

import collections
import re
import functools
import warnings

import numpy as np
from numpy.lib.utils import deprecate
from Bio.KDTree import KDTree

from MDAnalysis.core import flags
from ..lib import distances
from ..lib.mdamath import triclinic_vectors
from ..exceptions import SelectionError, NoDataError


def unique(ag):
    """Return the unique elements of ag"""
    try:
        return ag.universe.atoms[np.unique(ag.indices)]
    except NoDataError:
        # zero length AG has no Universe
        return ag


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
        return unique(rsel[mask])


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
        return unique(group[:])


class UnarySelection(Selection):
    def __init__(self, parser, tokens):
        sel = parser.parse_expression(self.precedence)
        self.sel = sel


class NotSelection(UnarySelection):
    token = 'not'
    precedence = 5

    def apply(self, group):
        notsel = self.sel.apply(group)
        return unique(group[~np.in1d(group.indices, notsel.indices)])


class GlobalSelection(UnarySelection):
    token = 'global'
    precedence = 5

    def apply(self, group):
        return unique(self.sel.apply(group.universe.atoms))


class ByResSelection(UnarySelection):
    token = 'byres'
    precedence = 1

    def apply(self, group):
        res = self.sel.apply(group)
        unique_res = np.unique(res.resids)
        mask = np.in1d(group.resids, unique_res)

        return unique(group[mask])


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
        # KDTree doesn't support periodic
        if self.periodic:
            self.apply = self._apply_distmat


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
        Limitations: always ignores periodicity
        """
        sel = self.sel.apply(group)
        # All atoms in group that aren't in sel
        sys = group[~np.in1d(group.indices, sel.indices)]

        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(sys.positions)
        found_indices = []
        for atom in sel.positions:
            kdtree.search(atom, self.cutoff)
            found_indices.append(kdtree.get_indices())
        # These are the indices from SYS that were seen when
        # probing with SEL
        unique_idx = np.unique(np.concatenate(found_indices))
        return unique(sys[unique_idx.astype(np.int32)])

    def _apply_distmat(self, group):
        sel = self.sel.apply(group)
        sys = group[~np.in1d(group.indices, sel.indices)]

        box = group.dimensions if self.periodic else None
        dist = distances.distance_array(
            sys.positions, sel.positions, box)

        mask = (dist <= self.cutoff).any(axis=1)

        return unique(sys[mask])


class SphericalLayerSelection(DistanceSelection):
    token = 'sphlayer'
    precedence = 1

    def __init__(self, parser, tokens):
        super(SphericalLayerSelection, self).__init__()
        self.inRadius = float(tokens.popleft())
        self.exRadius = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)

    def _apply_KDTree(self, group):
        """Selection using KDTree but periodic = True not supported.
        """
        sel = self.sel.apply(group)
        ref = sel.center_of_geometry()

        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(group.positions)

        kdtree.search(ref, self.exRadius)
        found_ExtIndices = kdtree.get_indices()
        kdtree.search(ref, self.inRadius)
        found_IntIndices = kdtree.get_indices()
        found_indices = list(set(found_ExtIndices) - set(found_IntIndices))
        return unique(group[found_indices])

    def _apply_distmat(self, group):
        sel = self.sel.apply(group)
        ref = sel.center_of_geometry().reshape(1, 3).astype(np.float32)

        box = group.dimensions if self.periodic else None
        d = distances.distance_array(ref,
                                     group.positions,
                                     box=box)[0]
        mask = d < self.exRadius
        mask &= d > self.inRadius

        return unique(group[mask])


class SphericalZoneSelection(DistanceSelection):
    token = 'sphzone'
    precedence = 1

    def __init__(self, parser, tokens):
        super(SphericalZoneSelection, self).__init__()
        self.cutoff = float(tokens.popleft())
        self.sel = parser.parse_expression(self.precedence)

    def _apply_KDTree(self, group):
        """Selection using KDTree but periodic = True not supported.
        (KDTree routine is ca 15% slower than the distance matrix one)
        """
        sel = self.sel.apply(group)
        ref = sel.center_of_geometry()

        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(group.positions)
        kdtree.search(ref, self.cutoff)
        found_indices = kdtree.get_indices()

        return unique(group[found_indices])

    def _apply_distmat(self, group):
        sel = self.sel.apply(group)
        ref = sel.center_of_geometry().reshape(1, 3).astype(np.float32)

        box = group.dimensions if self.periodic else None
        d = distances.distance_array(ref,
                                     group.positions,
                                     box=box)[0]
        idx = d < self.cutoff
        return unique(group[idx])


class CylindricalSelection(Selection):
    def __init__(self):
        self.periodic = flags['use_periodic_selections']

    def apply(self, group):
        sel = self.sel.apply(group)

        # Calculate vectors between point of interest and our group
        vecs = group.positions - sel.center_of_geometry()

        if self.periodic and not np.any(group.dimensions[:3]==0):
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
                #Triclinic version
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

        return unique(group[mask])


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
        self.ref = np.array([x, y, z])
        self.cutoff = float(tokens.popleft())

    def _apply_KDTree(self, group):
        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(group.positions)
        kdtree.search(self.ref, self.cutoff)
        found_indices = kdtree.get_indices()

        return unique(group[found_indices])

    def _apply_distmat(self, group):
        ref_coor = self.ref[np.newaxis, ...]

        ref_coor = np.asarray(ref_coor, dtype=np.float32)
        box = group.dimensions if self.periodic else None

        dist = distances.distance_array(group.positions, ref_coor, box)
        mask = (dist <= self.cutoff).any(axis=1)
        return unique(group[mask])


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
        return unique(sub)


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
        idx.append(bix[:,0][np.in1d(bix[:,1], grpidx)])
        # right side
        idx.append(bix[:,1][np.in1d(bix[:,0], grpidx)])

        idx = np.union1d(*idx)

        return group.universe.atoms[np.unique(idx)]


class SelgroupSelection(Selection):
    token = 'group'

    def __init__(self, parser, tokens):
        grpname = tokens.popleft()
        try:
            self.grp = parser.selgroups[grpname]
        except KeyError:
            raise ValueError("Failed to find group: {0}".format(grpname))

    def apply(self, group):
        mask = np.in1d(group.indices, self.grp.indices)
        return group[mask]


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
        return unique(self.grp)


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

        return unique(group[mask])


class AtomNameSelection(StringSelection):
    """Select atoms based on 'names' attribute"""
    token = 'name'
    field = 'names'


class AtomTypeSelection(StringSelection):
    """Select atoms based on 'types' attribute"""
    token = 'type'
    field = 'types'


class ResidueNameSelection(StringSelection):
    """Select atoms based on 'resnames' attribute"""
    token = 'resname'
    field = 'resnames'


class SegmentNameSelection(StringSelection):
    """Select atoms based on 'segids' attribute"""
    token = 'segid'
    field = 'segids'


class AltlocSelection(StringSelection):
    """Select atoms based on 'altLoc' attribute"""
    token = 'altloc'
    field = 'altLocs'


class RangeSelection(Selection):
    """Select atoms based on numerical fields

    Allows the use of ':' and '-' to specify a range of values
    For example

      resid 1:10
    """
    def __init__(self, parser, tokens):
        values = grab_not_keywords(tokens)
        if not values:
            raise ValueError("Unexpected token: '{0}'".format(tokens[0]))

        uppers = []
        lowers = []

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
                lower, upper = map(int, selrange.groups())

            lowers.append(lower)
            uppers.append(upper)

        self.lowers = lowers
        self.uppers = uppers

    def _get_vals(self, group):
        return getattr(group, self.field)

    def apply(self, group):
        mask = np.zeros(len(group), dtype=np.bool)
        vals = self._get_vals(group)

        for upper, lower in zip(self.uppers, self.lowers):
            if upper is not None:
                thismask = vals >= lower
                thismask &= vals <= upper
            else:
                thismask = vals == lower

            mask |= thismask
        return unique(group[mask])


class ResidueIDSelection(RangeSelection):
    token = 'resid'
    field = 'resids'


class ResnumSelection(RangeSelection):
    token = 'resnum'
    field = 'resnums'


class ByNumSelection(RangeSelection):
    token = 'bynum'

    def _get_vals(self, group):
        # In this case we'll use 1 indexing since that's what the
        # user will be familiar with
        return group.indices + 1


class ProteinSelection(Selection):
    """Consists of all residues with  recognized residue names.

    Recognized residue names in :attr:`ProteinSelection.prot_res`.

      * from the CHARMM force field::
         awk '/RESI/ {printf "'"'"%s"'"',",$2 }' top_all27_prot_lipid.rtf

      * manually added special CHARMM, OPLS/AA and Amber residue names.

      * still missing: Amber N- and C-terminal residue names

    .. SeeAlso:: :func:`MDAnalysis.lib.util.convert_aa_code`
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
        # Amber: there are also the C-term aas: 'CALA', 'CGLY', 'CSER', ...
        # Amber: there are also the N-term aas: 'NALA', 'NGLY', 'NSER', ...
        'HID', 'HIE', 'HIP', 'ORN', 'DAB', 'LYN', 'HYP', 'CYM', 'CYX', 'ASH',
        'GLH',
        'ACE', 'NME',
    ])
    def __init__(self, parser, tokens):
        pass

    def apply(self, group):
        mask = np.in1d(group.resnames, self.prot_res)
        return unique(group[mask])


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
        return unique(group[mask])


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
        return unique(group[mask])


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
        return unique(group[mask])


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
        return unique(group[mask])


class NucleicSugarSelection(NucleicSelection):
    """Contains all atoms with name C1', C2', C3', C4', O2', O4', O3'.
    """
    token = 'nucleicsugar'
    sug_atoms = np.array(["C1'", "C2'", "C3'", "C4'", "O4'"])

    def apply(self, group):
        mask = np.in1d(group.names, self.sug_atoms)
        mask &= np.in1d(group.resnames, self.nucl_res)
        return unique(group[mask])


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

    def __init__(self, parser, tokens):
        prop = tokens.popleft()
        if prop == "abs":
            self.absolute = True
            prop = tokens.popleft()
        else:
            self.absolute = False
        oper = tokens.popleft()
        self.value = float(tokens.popleft())

        self.prop = prop
        try:
            self.operator = self.ops[oper]
        except KeyError:
            raise ValueError(
                "Invalid operator : '{0}' Use one of : '{1}'"
                "".format(oper, self.ops.keys()))

    def apply(self, group):
        try:
            col = {'x':0, 'y':1, 'z':2}[self.prop]
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

        return unique(group[mask])


class SameSelection(Selection):
    token = 'same'
    precedence = 1

    prop_trans = {
        'fragment': None,
        'x': None,
        'y': None,
        'z': None,
        'residue':'resids',
        'segment':'segids',
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
        if not prop in self.prop_trans:
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
            return unique(group[mask])
        # [xyz] must come before self.prop_trans lookups too!
        try:
            pos_idx = {'x':0, 'y':1, 'z':2}[self.prop]
        except KeyError:
            # The self.prop string was already checked,
            # so don't need error checking here.
            # KeyError at this point is impossible!
            attrname = self.prop_trans[self.prop]
            vals = getattr(res, attrname)
            mask = np.in1d(getattr(group, attrname), vals)

            return unique(group[mask])
        else:
            vals = res.positions[:, pos_idx]
            pos = group.positions[:, pos_idx]

            # isclose only does one value at a time
            mask = np.vstack([np.isclose(pos, v)
                              for v in vals]).any(axis=0)
            return unique(group[mask])


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
        if not self.tokens[0] == None:
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
