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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Atom selection Hierarchy --- :mod:`MDAnalysis.core.selection`
=============================================================

These objects are constructed and applied to the group

In general, Parser.parse() creates a Selection object
from a selection string.

This Selection object is then passed an AtomGroup through its
apply method to apply the Selection to the AtomGroup.

This is all invisible to the user through ag.select_atoms
"""

import re
import numpy as np
from numpy.lib.utils import deprecate
from Bio.KDTree import KDTree
import warnings
import logging

from MDAnalysis.core import flags
from ..lib import distances
from ..lib.mdamath import triclinic_vectors
from ..exceptions import SelectionError


logger = logging.getLogger(__name__)


class AllSelection(object):
    def apply(self, group):
        return group[:]


class NotSelection(object):
    def __init__(self, sel):
        self.sel = sel

    def apply(self, group):
        notsel = self.sel.apply(group)
        return group[~np.in1d(group.indices, notsel.indices)]


class AndSelection(object):
    def __init__(self, lsel, rsel):
        self.rsel = rsel
        self.lsel = lsel

    def apply(self, group):
        rsel = self.rsel.apply(group)
        lsel = self.lsel.apply(group)

        # Mask which lsel indices appear in rsel
        mask = np.in1d(rsel.indices, lsel.indices)
        # and mask rsel according to that
        return rsel[mask]


class OrSelection(object):
    def __init__(self, lsel, rsel):
        self.rsel = rsel
        self.lsel = lsel

    def apply(self, group):
        lsel = self.lsel.apply(group)
        rsel = self.rsel.apply(group)

        # Find unique indices from both these AtomGroups
        # and slice master list using them
        idx = np.unique(np.concatenate([lsel.indices, rsel.indices]))

        return group.universe.atoms[idx]


class GlobalSelection(object):
    def __init__(self, sel):
        self.sel = sel

    def apply(self, group):
        return self.sel.apply(group.universe.atoms)


class DistanceSelection(object):
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
    def __init__(self, sel, cutoff):
        super(AroundSelection, self).__init__()
        self.sel = sel
        self.cutoff = cutoff

    def _apply_KDTree(self, group):
        """KDTree based selection is about 7x faster than distmat
        for typical problems.
        Limitations: always ignores periodicity
        """
        logger.debug("In Around KDTree")
        sel = self.sel.apply(group)
        logger.debug("Reference group is {}".format(sel))
        # All atoms in group that aren't in sel
        sys = group[~np.in1d(group.indices, sel.indices)]
        logger.debug("Other group is {}".format(sys))

        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(sys.positions)
        found_indices = []
        for atom in sel.positions:
            kdtree.search(atom, self.cutoff)
            found_indices.append(kdtree.get_indices())
        logger.debug("Found indices are {}".format(found_indices))
        # These are the indices from SYS that were seen when
        # probing with SEL
        unique_idx = np.unique(np.concatenate(found_indices))
        logger.debug("Unique indices are {}".format(unique_idx))
        return sys[unique_idx]

    def _apply_distmat(self, group):
        logger.debug("In Around Distmat")
        sel = self.sel.apply(group)
        logger.debug("Sel is {}".format(sel))
        sys = group[~np.in1d(group.indices, sel.indices)]
        logger.debug("Sys is {}".format(sys))

        box = group.dimensions if self.periodic else None
        dist = distances.distance_array(
            sys.positions, sel.positions, box)
        logger.debug("dist has shape {}".format(dist.shape))
        logger.debug("dist is {}".format(dist))

        mask = (dist <= self.cutoff).any(axis=1)

        logger.debug("mask has shape {}".format(mask.shape))

        return sys[mask]


class SphericalLayerSelection(DistanceSelection):
    def __init__(self, sel, inRadius, exRadius):
        super(SphericalLayerSelection, self).__init__()
        self.sel = sel
        self.inRadius = inRadius
        self.exRadius = exRadius

    def _apply_KDTree(self, group):
        """Selection using KDTree but periodic = True not supported.
        """
        sel = self.sel.apply(group)
        ref = sel.center_of_geometry()
        sys = group[~np.in1d(group.indices, sel.indices)]

        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(sys.positions)

        kdtree.search(ref, self.exRadius)
        found_ExtIndices = kdtree.get_indices()
        kdtree.search(ref, self.inRadius)
        found_IntIndices = kdtree.get_indices()
        found_indices = list(set(found_ExtIndices) - set(found_IntIndices))
        return sys[found_indices]

    def _apply_distmat(self, group):
        sel = self.sel.apply(group)
        sel_CoG = sel.center_of_geometry().reshape(1, 3).astype(np.float32)

        sys = group[~np.in1d(group.indices, sel.indices)]

        box = group.dimensions if self.periodic else None
        d = distances.distance_array(sel_CoG,
                                     sys.positions,
                                     box=box)[0]
        lmask = d < self.exRadius
        rmask = d > self.inRadius
        mask = lmask & rmask
        return sys[mask]


class SphericalZoneSelection(DistanceSelection):
    def __init__(self, sel, cutoff):
        super(SphericalZoneSelection, self).__init__()
        self.sel = sel
        self.cutoff = cutoff

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

        return sys[found_indices]

    def _apply_distmat(self, group):
        sel = self.sel.apply(group)
        ref = sel.center_of_geometry().reshape(1, 3).astype(np.float32)

        box = group.dimensions if self.periodic else None
        d = distances.distance_array(ref,
                                     group.positions,
                                     box=box)[0]
        idx = d < self.cutoff
        return group[idx]


class CylindricalSelection(object):
    def __init__(self, sel, exRadius, zmax, zmin):
        self.sel = sel
        self.exRadius = exRadius
        self.exRadiusSq = exRadius * exRadius
        self.zmax = zmax
        self.zmin = zmin
        self.periodic = flags['use_periodic_selections']

    def apply(self, group):
        sel = self.sel.apply(group)
        sel_CoG = sel.center_of_geometry()
        coords = group.positions

        if self.periodic and not np.any(group.dimensions[:3]==0):
            box = group.dimensions
            cyl_z_hheight = (self.zmax-self.zmin)/2

            if 2*self.exRadius > box[0]:
                raise NotImplementedError(
                    "The diameter of the cylinder selection ({:.3f}) is larger "
                    "than the unit cell's x dimension ({:.3f}). Can only do "
                    "selections where it is smaller or equal."
                    "".format(2*self.exRadius, box[0]))
            if 2*self.exRadius > box[1]:
                raise NotImplementedError(
                    "The diameter of the cylinder selection ({:.3f}) is larger "
                    "than the unit cell's y dimension ({:.3f}). Can only do "
                    "selections where it is smaller or equal."
                    "".format(2*self.exRadius, box[1]))
            if 2*cyl_z_hheight > box[2]:
                raise NotImplementedError(
                    "The total length of the cylinder selection in z ({:.3f}) "
                    "is larger than the unit cell's z dimension ({:.3f}). Can "
                    "only do selections where it is smaller or equal."
                    "".format(2*cyl_z_hheight, box[2]))

            #how off-center in z is our CoG relative to the cylinder's center
            cyl_center = sel_CoG + [0,0,(self.zmax+self.zmin)/2]
            coords += box[:3] / 2 - cyl_center
            coords = distances.apply_PBC(coords, box=box)

            sel_CoG = box[:3] / 2
            zmin = -cyl_z_hheight
            zmax = cyl_z_hheight
        else:
            zmin = self.zmin
            zmax = self.zmax

        # For performance we first do the selection of the atoms in the
        # rectangular parallelepiped that contains the cylinder.
        lim_min = sel_CoG - [self.exRadius, self.exRadius, -zmin]
        lim_max = sel_CoG + [self.exRadius, self.exRadius, zmax]
        mask_sel = np.all((coords >= lim_min) * (coords <= lim_max), axis=1)
        mask_ndxs = np.where(mask_sel)[0]
        # Now we do the circular part
        xy_vecs = coords[mask_ndxs,:2] - sel_CoG[:2]
        xy_norms = np.sum(xy_vecs**2, axis=1)
        try: # Generic for both 'Layer' and 'Zone' cases
            circ_sel = (xy_norms <= self.exRadiusSq) & (xy_norms >= self.inRadiusSq)
        except AttributeError:
            circ_sel = (xy_norms <= self.exRadiusSq)
        mask_sel[mask_ndxs] = circ_sel
        return group[mask_sel]


class CylindricalZoneSelection(CylindricalSelection):
    def __init__(self, sel, exRadius, zmax, zmin):
        super(CylindricalZoneSelection, self).__init__(
            sel, exRadius, zmax, zmin)


class CylindricalLayerSelection(CylindricalSelection):
    def __init__(self, sel, inRadius, exRadius, zmax, zmin):
        super(CylindricalLayerSelection, self).__init__(
            sel, exRadius, zmax, zmin)
        self.inRadius = inRadius
        self.inRadiusSq = inRadius * inRadius


class PointSelection(DistanceSelection):
    def __init__(self, x, y, z, cutoff):
        super(PointSelection, self).__init__()
        self.ref = np.array([x, y, z])
        self.cutoff = cutoff

    def _apply_KDTree(self, group):
        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(group.positions)
        kdtree.search(self.ref, self.cutoff)
        found_indices = kdtree.get_indices()

        return group[found_indices]

    def _apply_distmat(self, group):
        ref_coor = self.ref[np.newaxis, ...]

        ref_coor = np.asarray(ref_coor, dtype=np.float32)
        box = group.dimensions if self.periodic else None

        dist = distances.distance_array(group.positions, ref_coor, box)
        mask = (dist <= self.cutoff).any(axis=1)
        return group[mask]


class AtomSelection(object):
    def __init__(self, name, resid, segid):
        self.name = name
        self.resid = resid
        self.segid = segid

    def apply(self, group):
        sub = group[group.names == self.name]
        sub = sub[sub.resids == self.resid]
        sub = sub[sub.segids == self.segid]
        return sub


class BondedSelection(object):
    def __init__(self, sel):
        self.sel = sel

    def apply(self, group):
        grp = self.sel.apply(group)
        # Check if we have bonds
        if not group.bonds:
            warnings.warn("Bonded selection has 0 bonds")

        grpidx = grp.indices

        # (n, 2) array of bond indices
        bix = group.bonds.indices

        idx = []
        # left side
        idx.append(bix[:,0][np.in1d(bix[:,1], grpidx)])
        # right side
        idx.append(bix[:,1][np.in1d(bix[:,0], grpidx)])

        idx = np.union1d(*idx)

        return group.universe.atoms[idx]


class SelgroupSelection(object):
    def __init__(self, selgroup):
        self.grp = selgroup

    def apply(self, group):
        idx = np.intersect1d(self.grp.indices, group.indices)
        return group.universe.atoms[idx]


@deprecate(old_name='fullgroup', new_name='global group')
class FullSelgroupSelection(object):
    def __init__(self, selgroup):
        self.grp = selgroup

    def apply(self, group):
        return self.grp


class StringSelection(object):
    def __init__(self, val):
        self.val = val

    def apply(self, group):
        wc_pos = self.val.find('*')
        if wc_pos == -1:  # No wildcard found
            mask = getattr(group, self.field) == self.val
        else:
            values = getattr(group, self.field).astype(np.str_)
            mask = np.char.startswith(values, self.val[:wc_pos])

        return group[mask]


class AtomNameSelection(StringSelection):
    field = 'names'


class AtomTypeSelection(StringSelection):
    field = 'types'


class ResidueNameSelection(StringSelection):
    field = 'resnames'


class SegmentNameSelection(StringSelection):
    field = 'segids'


class AltlocSelection(StringSelection):
    field = 'altLocs'


class ByResSelection(object):
    def __init__(self, sel):
        self.sel = sel

    def apply(self, group):
        logger.debug("In ByResSelection")
        res = self.sel.apply(group)
        logger.debug("Reference group is {}".format(res))

        unique_res = np.unique(res.resindices)
        logger.debug("Ref resindices are {}".format(unique_res))
        mask = np.in1d(group.resindices, unique_res)

        return group[mask]


class RangeSelection(object):
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper


class ResidueIDSelection(RangeSelection):
    def apply(self, group):
        resids = group.resids
        if self.upper is not None:
            lower_mask = resids >= self.lower
            upper_mask = resids <= self.upper
            mask = lower_mask & upper_mask
        else:
            mask = resids == self.lower
        return group[mask]


class ResnumSelection(RangeSelection):
    def apply(self, group):
        resnums = group.resnums
        if self.upper is not None:
            lower_mask = resnums >= self.lower
            upper_mask = resnums <= self.upper
            mask = lower_mask & upper_mask
        else:
            mask = resnums == self.lower
        return group[mask]


class ByNumSelection(RangeSelection):
    def apply(self, group):
        # In this case we'll use 1 indexing since that's what the
        # user will be familiar with
        indices = group.indices + 1
        if self.upper is not None:
            lower_mask = indices >= self.lower
            upper_mask = indices <= self.upper
            mask = lower_mask & upper_mask
        else:
            mask = indices == self.lower
        return group[mask]


class ProteinSelection(object):
    """Consists of all residues with  recognized residue names.

    Recognized residue names in :attr:`ProteinSelection.prot_res`.

      * from the CHARMM force field::
         awk '/RESI/ {printf "'"'"%s"'"',",$2 }' top_all27_prot_lipid.rtf

      * manually added special CHARMM, OPLS/AA and Amber residue names.

      * still missing: Amber N- and C-terminal residue names

    .. SeeAlso:: :func:`MDAnalysis.lib.util.convert_aa_code`
    """
    prot_res = [
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
    ]

    def apply(self, group):
        mask = np.in1d(group.resnames, self.prot_res)
        return group[mask]


class NucleicSelection(object):
    """All atoms in nucleic acid residues with recognized residue names.

    Recognized residue names:

    * from the CHARMM force field ::
        awk '/RESI/ {printf "'"'"%s"'"',",$2 }' top_all27_prot_na.rtf
    * recognized: 'ADE', 'URA', 'CYT', 'GUA', 'THY'
    * recognized (CHARMM in Gromacs): 'DA', 'DU', 'DC', 'DG', 'DT'

    .. versionchanged:: 0.8
       additional Gromacs selections
    """
    nucl_res = [
        'ADE', 'URA', 'CYT', 'GUA', 'THY', 'DA', 'DC', 'DG', 'DT', 'RA',
        'RU', 'RG', 'RC', 'A', 'T', 'U', 'C', 'G']

    def apply(self, group):
        mask = np.in1d(group.resnames, self.nucl_res)
        return group[mask]


class BackboneSelection(ProteinSelection):
    """A BackboneSelection contains all atoms with name 'N', 'CA', 'C', 'O'.

    This excludes OT* on C-termini
    (which are included by, eg VMD's backbone selection).
    """
    bb_atoms = ['N', 'CA', 'C', 'O']

    def apply(self, group):
        namemask = np.in1d(group.names, self.bb_atoms)
        resnamemask = np.in1d(group.resnames, self.prot_res)
        mask = namemask & resnamemask
        return group[mask]


class NucleicBackboneSelection(ProteinSelection):
    """Contains all atoms with name "P", "C5'", C3'", "O3'", "O5'".

    These atoms are only recognized if they are in a residue matched
    by the :class:`NucleicSelection`.
    """
    bb_atoms = ["P", "C5'", "C3'", "O3'", "O5'"]

    def apply(self, group):
        namemask = np.in1d(group.names, self.bb_atoms)
        resnamemask = np.in1d(group.resnames, self.prot_res)
        mask = namemask & resnamemask
        return group[mask]


class BaseSelection(NucleicSelection):
    """Selection of atoms in nucleobases.

    Recognized atom names (from CHARMM):

     'N9', 'N7', 'C8', 'C5', 'C4', 'N3', 'C2', 'N1', 'C6',
     'O6','N2','N6', 'O2','N4','O4','C5M'
    """
    base_atoms = [
        'N9', 'N7', 'C8', 'C5', 'C4', 'N3', 'C2', 'N1', 'C6',
        'O6', 'N2', 'N6',
        'O2', 'N4', 'O4', 'C5M']

    def apply(self, group):
        namemask = np.in1d(group.names, self.base_atoms)
        resnamemask = np.in1d(group.resnames, self.nucl_res)
        mask = namemask & resnamemask
        return group[mask]


class NucleicSugarSelection(NucleicSelection):
    """Contains all atoms with name C1', C2', C3', C4', O2', O4', O3'.
    """
    sug_atoms = ["C1'", "C2'", "C3'", "C4'", "O4'"]

    def apply(self, group):
        namemask = np.in1d(group.names, self.sug_atoms)
        resnamemask = np.in1d(group.resnames, self.nucl_res)
        mask = namemask & resnamemask
        return group[mask]


class PropertySelection(object):
    """Some of the possible properties:
    x, y, z, radius, mass,
    """
    def __init__(self, prop, operator, value, abs_=False):
        self.prop = prop
        self.operator = operator
        self.value = value
        self.abs_ = abs_

    def apply(self, group):
        try:
            col = {'x':0, 'y':1, 'z':2}[self.prop]
        except KeyError:
            pass
        else:
            values = group.positions[:,col]

        if self.prop == 'mass':
            values = group.masses
        elif self.prop == 'charge':
            values = group.charges

        if self.abs_:
            values = np.abs(values)
        mask = self.operator(values, self.value)

        return group[mask]


class SameSelection(object):
    # When adding new keywords here don't forget to also add them to the
    #  case statement under the SAME op, where they are first checked.
    def __init__(self, sel, prop):
        self.sel = sel
        self.prop = prop

    def apply(self, group):
        prop_trans = {
            'residue':'resindices',
            'segment':'segindices',
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

        res = self.sel.apply(group)
        if not res:
            return group[[]]  # empty selection

        # TODO: Fragments!

        try:
            attrname = prop_trans[self.prop]
        except KeyError:
            pass
        else:
            vals = getattr(res, attrname)
            mask = np.in1d(getattr(group, attrname), vals)

            return group[mask]

        try:
            pos_idx = {0:'x', 1:'y', 2:'z'}[self.prop]
        except KeyError:
            pass
        else:
            vals = res.positions[:, pos_idx]
            mask = np.in1d(group.positions[:, pos_idx], vals)

            return group[mask]


class ParseError(Exception):
    pass


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

    #Here are the symbolic tokens that we'll process:
    ALL = 'all'
    NOT = 'not'
    AND = 'and'
    OR = 'or'
    AROUND = 'around'
    SPHLAYER = 'sphlayer'
    SPHZONE = 'sphzone'
    CYLAYER = 'cylayer'
    CYZONE = 'cyzone'
    POINT = 'point'
    BONDED = 'bonded'
    BYRES = 'byres'
    BYNUM = 'bynum'
    PROP = 'prop'
    ATOM = 'atom'
    LPAREN = '('
    RPAREN = ')'
    SEGID = 'segid'
    ALTLOC = 'altloc'
    RESID = 'resid'
    RESNUM = 'resnum'
    RESNAME = 'resname'
    NAME = 'name'
    TYPE = 'type'
    PROTEIN = 'protein'
    NUCLEIC = 'nucleic'
    BB = 'backbone'
    NBB = 'nucleicbackbone'
    BASE = 'nucleicbase'
    SUGAR = 'nucleicsugar'
    SAME = 'same'
    GLOBAL = 'global'
    EOF = 'EOF'
    GT = '>'
    LT = '<'
    GE = '>='
    LE = '<='
    EQ = '=='
    NE = '!='
    SELGROUP = 'group'
    FULLSELGROUP = 'fullgroup'

    classdict = dict([
        (ALL, AllSelection),
        (NOT, NotSelection),
        (AND, AndSelection),
        (BONDED, BondedSelection),
        (OR, OrSelection),
        (SEGID, SegmentNameSelection),
        (RESID, ResidueIDSelection),
        (RESNUM, ResnumSelection),
        (RESNAME, ResidueNameSelection),
        (NAME, AtomNameSelection),
        (ALTLOC, AltlocSelection),
        (TYPE, AtomTypeSelection),
        (BYRES, ByResSelection),
        (BYNUM, ByNumSelection),
        (PROP, PropertySelection),
        (AROUND, AroundSelection),
        (SPHLAYER, SphericalLayerSelection),
        (SPHZONE, SphericalZoneSelection),
        (CYLAYER, CylindricalLayerSelection),
        (CYZONE, CylindricalZoneSelection),
        (POINT, PointSelection),
        (NUCLEIC, NucleicSelection),
        (PROTEIN, ProteinSelection),
        (BB, BackboneSelection),
        (NBB, NucleicBackboneSelection),
        (BASE, BaseSelection),
        (SUGAR, NucleicSugarSelection),
        (ATOM, AtomSelection),
        (SELGROUP, SelgroupSelection),
        (FULLSELGROUP, FullSelgroupSelection),
        (SAME, SameSelection),
        (GLOBAL, GlobalSelection)])
    precedence = dict(
        [(AROUND, 1),
         (SPHLAYER, 1),
         (SPHZONE, 1),
         (CYLAYER, 1),
         (CYZONE, 1),
         (POINT, 1),
         (BYRES, 1),
         (BONDED, 1),
         (SAME, 1),
         (AND, 3),
         (OR, 3),
         (NOT, 5),
         (GLOBAL, 5)])

    # Borg pattern: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66531
    _shared_state = {}

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def _peek_token(self):
        """Looks at the next token in our token stream."""
        return self.tokens[0]

    def _consume_token(self):
        """Pops off the next token in our token stream."""
        return self.tokens.pop(0)

    def _error(self, token):
        """Stops parsing and reports an error."""
        raise ParseError("Parsing error- '{}'\n {} expected"
                         "".format(self.selectstr, token))

    def _expect(self, token):
        if self._peek_token() == token:
            self._consume_token()
        else:
            self._error(token)

    def parse(self, selectstr, selgroups):
        self.selectstr = selectstr
        self.selgroups = selgroups
        self.tokens = selectstr.replace('(', ' ( ').replace(')', ' ) ')
        self.tokens = self.tokens.split() + [self.EOF]
        parsetree = self._parse_expression(0)
        self._expect(self.EOF)
        return parsetree

    def _parse_expression(self, p):
        exp1 = self._parse_subexp()
        while (self._peek_token() in (self.AND, self.OR) and
               self.precedence[self._peek_token()] >= p):  # bin ops
            op = self._consume_token()
            q = 1 + self.precedence[op]
            exp2 = self._parse_expression(q)
            exp1 = self.classdict[op](exp1, exp2)
        return exp1

    def _parse_subexp(self):
        op = self._consume_token()
        if op in (self.NOT, self.BYRES, self.GLOBAL):  # unary operators
            exp = self._parse_expression(self.precedence[op])
            return self.classdict[op](exp)
        elif op == self.AROUND:
            dist = self._consume_token()
            exp = self._parse_expression(self.precedence[op])
            return self.classdict[op](exp, float(dist))
        elif op == self.SPHLAYER:
            inRadius = self._consume_token()
            exRadius = self._consume_token()
            exp = self._parse_expression(self.precedence[op])
            return self.classdict[op](exp, float(inRadius), float(exRadius))
        elif op == self.SPHZONE:
            dist = self._consume_token()
            exp = self._parse_expression(self.precedence[op])
            return self.classdict[op](exp, float(dist))
        elif op == self.CYLAYER:
            inRadius = self._consume_token()
            exRadius = self._consume_token()
            zmax = self._consume_token()
            zmin = self._consume_token()
            exp = self._parse_expression(self.precedence[op])
            return self.classdict[op](exp,
                                      float(inRadius), float(exRadius),
                                      float(zmax), float(zmin))
        elif op == self.CYZONE:
            exRadius = self._consume_token()
            zmax = self._consume_token()
            zmin = self._consume_token()
            exp = self._parse_expression(self.precedence[op])
            return self.classdict[op](exp,
                                      float(exRadius),
                                      float(zmax), float(zmin))
        elif op == self.POINT:
            dist = self._consume_token()
            x = self._consume_token()
            y = self._consume_token()
            z = self._consume_token()
            return self.classdict[op](float(dist), float(x), float(y), float(z))
        elif op == self.BONDED:
            exp = self._parse_expression(self.precedence[op])
            return self.classdict[op](exp)
        elif op == self.LPAREN:
            exp = self._parse_expression(0)
            self._expect(self.RPAREN)
            return exp
        elif op in {self.SEGID, self.RESNAME, self.NAME,
                    self.TYPE, self.ALTLOC}:
            data = self._consume_token()
            if data in {
                    self.LPAREN, self.RPAREN,
                    self.AND, self.OR, self.NOT,
                    self.SEGID, self.RESID, self.RESNAME,
                    self.NAME, self.TYPE}:
                self._error("Identifier")
            return self.classdict[op](data)
        elif op in (self.SELGROUP, self.FULLSELGROUP):
            grpname = self._consume_token()
            return self.classdict[op](self.selgroups[grpname])
        elif op == self.PROTEIN:
            return self.classdict[op]()
        elif op == self.NUCLEIC:
            return self.classdict[op]()
        elif op == self.ALL:
            return self.classdict[op]()
        elif op == self.BB:
            return self.classdict[op]()
        elif op == self.NBB:
            return self.classdict[op]()
        elif op == self.BASE:
            return self.classdict[op]()
        elif op == self.SUGAR:
            return self.classdict[op]()
        elif op in (self.RESID, self.RESNUM, self.BYNUM):
            # can operate on ranges X:Y or X-Y
            data = self._consume_token()
            try:
                lower = int(data)
                upper = None
            except ValueError:
                # check if in appropriate format 'lower:upper' or 'lower-upper'
                selrange = re.match("(\d+)[:-](\d+)",
                                    data)
                if not selrange:
                    self._error(op)
                lower, upper = map(int, selrange.groups())
            return self.classdict[op](lower, upper)
        elif op == self.PROP:
            prop = self._consume_token()
            if prop == "abs":
                abs_ = True
                prop = self._consume_token()
            else:
                abs_ = False
            oper = self._consume_token()
            value = float(self._consume_token())
            ops = dict([
                (self.GT, np.greater),
                (self.LT, np.less),
                (self.GE, np.greater_equal),
                (self.LE, np.less_equal),
                (self.EQ, np.equal),
                (self.NE, np.not_equal),
            ])
            return self.classdict[op](prop, ops[oper], value, abs_)
        elif op == self.ATOM:
            segid = self._consume_token()
            resid = int(self._consume_token())
            name = self._consume_token()
            return self.classdict[op](name, resid, segid)
        elif op == self.SAME:
            prop = self._consume_token()
            self._expect("as")
            if prop in {"name", "type", "resname", "resid", "segid",
                        "mass", "charge", "radius", "bfactor",
                        "resnum", "residue", "segment", "fragment",
                        "x", "y", "z"}:
                exp = self._parse_expression(self.precedence[op])
                return self.classdict[op](exp, prop)
            else:
                self._error(prop)
        else:
            self._error(op)

# The module level instance
Parser = SelectionParser()
