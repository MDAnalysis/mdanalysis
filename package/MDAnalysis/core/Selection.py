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
Atom selection Hierarchy --- :mod:`MDAnalysis.core.Selection`
======================================================================

These objects are constructed and applied to the group

Currently all atom arrays are handled internally as sets, but returned as AtomGroups

"""

import re
import numpy as np
from numpy.lib.utils import deprecate
from Bio.KDTree import KDTree

from .AtomGroup import AtomGroup, Universe
from MDAnalysis.core import flags
from ..lib import distances
from ..lib.mdamath import triclinic_vectors


class Selection:
    def __init__(self):
        # This allows you to build a Selection without tying it to a particular group yet
        # Updatable means every timestep
        self.update = False  # not used at the moment

    def __repr__(self):
        return "<" + self.__class__.__name__ + ">"

    def __and__(self, other):
        return AndSelection(self, other)

    def __or__(self, other):
        return OrSelection(self, other)

    def __invert__(self):
        return NotSelection(self)

    def __hash__(self):
        return hash(repr(self))

    def _apply(self, group):
        # This is an error
        raise NotImplementedError("No _apply function defined for " + repr(self.__class__.__name__))

    def apply(self, group):
        # Cache the result for future use
        # atoms is from Universe
        # returns AtomGroup
        if not (isinstance(group, Universe) or isinstance(group, AtomGroup)):
            raise Exception("Must pass in an AtomGroup or Universe to the Selection")
        # make a set of all the atoms in the group
        # XXX this should be static to all the class members
        Selection._group_atoms = set(group.atoms)
        Selection._group_atoms_list = [
            a for a in Selection._group_atoms
        ]  # need ordered, unique list for back-indexing in Around and Point!
        if not hasattr(group, "coord"):
            Selection.coord = group.universe.coord
        else:
            Selection.coord = group.coord

        if not hasattr(self, "_cache"):
            cache = list(self._apply(group))
            # Decorate/Sort/Undecorate (Schwartzian Transform)
            cache[:] = [(x.index, x) for x in cache]
            cache.sort()
            cache[:] = [val for (key, val) in cache]
            self._cache = AtomGroup(cache)
        return self._cache


class AllSelection(Selection):
    def __init__(self):
        Selection.__init__(self)

    def _apply(self, group):
        return set(group.atoms[:])


class NotSelection(Selection):
    def __init__(self, sel):
        Selection.__init__(self)
        self.sel = sel

    def _apply(self, group):
        notsel = self.sel._apply(group)
        return (set(group.atoms[:]) - notsel)

    def __repr__(self):
        return "<'NotSelection' " + repr(self.sel) + ">"


class AndSelection(Selection):
    def __init__(self, lsel, rsel):
        Selection.__init__(self)
        self.rsel = rsel
        self.lsel = lsel

    def _apply(self, group):
        return self.lsel._apply(group) & self.rsel._apply(group)

    def __repr__(self):
        return "<'AndSelection' " + repr(self.lsel) + "," + repr(self.rsel) + ">"


class OrSelection(Selection):
    def __init__(self, lsel, rsel):
        Selection.__init__(self)
        self.rsel = rsel
        self.lsel = lsel

    def _apply(self, group):
        return self.lsel._apply(group) | self.rsel._apply(group)

    def __repr__(self):
        return "<'OrSelection' " + repr(self.lsel) + "," + repr(self.rsel) + ">"

class GlobalSelection(Selection):
    def __init__(self, sel):
        Selection.__init__(self)
        self.sel = sel

    def _apply(self, group):
        sel = self.sel._apply(group.universe)
        return sel

class AroundSelection(Selection):
    def __init__(self, sel, cutoff, periodic=None):
        Selection.__init__(self)
        self.sel = sel
        self.cutoff = cutoff
        self.sqdist = cutoff * cutoff
        if periodic is None:
            self.periodic = flags['use_periodic_selections']

    def _apply(self, group):
        # make choosing _fast/_slow configurable (while testing)
        if flags['use_KDTree_routines'] in (True, 'fast', 'always'):
            return self._apply_KDTree(group)
        else:
            return self._apply_distmat(group)

    def _apply_KDTree(self, group):
        """KDTree based selection is about 7x faster than distmat for typical problems.
        Limitations: always ignores periodicity
        """
        # group is wrong, should be universe (?!)
        sel_atoms = self.sel._apply(group)
        # list needed for back-indexing
        sys_atoms_list = [a for a in (self._group_atoms - sel_atoms)]
        sel_indices = np.array([a.index for a in sel_atoms], dtype=int)
        sys_indices = np.array([a.index for a in sys_atoms_list], dtype=int)
        sel_coor = Selection.coord[sel_indices]

        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(Selection.coord[sys_indices])
        found_indices = []
        for atom in np.array(sel_coor):
            kdtree.search(atom, self.cutoff)
            found_indices.append(kdtree.get_indices())

        # the list-comprehension here can be understood as a nested loop.
        # for list in found_indices:
        #     for i in list:
        #         yield sys_atoms_list[i]
        # converting found_indices to a np array won't reallt work since
        # each we will find a different number of neighbors for each center in
        # sel_coor.
        res_atoms = [sys_atoms_list[i] for list in found_indices for i in list]
        return set(res_atoms)

    def _apply_distmat(self, group):
        sel_atoms = self.sel._apply(group)  # group is wrong, should be universe (?!)
        sys_atoms_list = [a for a in (self._group_atoms - sel_atoms)]  # list needed for back-indexing
        sel_indices = np.array([a.index for a in sel_atoms], dtype=int)
        sys_indices = np.array([a.index for a in sys_atoms_list], dtype=int)
        sel_coor = Selection.coord[sel_indices]
        sys_coor = Selection.coord[sys_indices]
        if self.periodic:
            box = group.dimensions[:3]  # ignored with KDTree
        else:
            box = None

        dist = distances.distance_array(sys_coor, sel_coor, box)
        res_atoms = [
            sys_atoms_list[i] for i in
            np.any(dist <= self.cutoff, axis=1).nonzero()[0]]  # make list np array and use fancy indexing?
        return set(res_atoms)

    def __repr__(self):
        return "<'AroundSelection' " + repr(self.cutoff) + " around " + repr(self.sel) + ">"


class SphericalLayerSelection(Selection):
    def __init__(self, sel, inRadius, exRadius, periodic=None):
        Selection.__init__(self)
        self.sel = sel
        self.inRadius = inRadius
        self.exRadius = exRadius
        if periodic is None:
            self.periodic = flags['use_periodic_selections']

    def _apply(self, group):
        # make choosing _fast/_slow configurable (while testing)
        if flags['use_KDTree_routines'] in (True, 'fast', 'always'):
            return self._apply_KDTree(group)
        else:
            return self._apply_distmat(group)

    def _apply_KDTree(self, group):
        """Selection using KDTree but periodic = True not supported.
        """
        sys_indices = np.array([a.index for a in self._group_atoms_list])
        sys_coor = Selection.coord[sys_indices]
        # group is wrong, should be universe (?!)
        sel_atoms = self.sel._apply(group)
        sel_CoG = AtomGroup(sel_atoms).center_of_geometry()
        self.ref = np.array((sel_CoG[0], sel_CoG[1], sel_CoG[2]))
        if self.periodic:
            pass  # or warn? -- no periodic functionality with KDTree search
        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(sys_coor)

        kdtree.search(self.ref, self.exRadius)
        found_ExtIndices = kdtree.get_indices()
        kdtree.search(self.ref, self.inRadius)
        found_IntIndices = kdtree.get_indices()
        found_indices = list(set(found_ExtIndices) - set(found_IntIndices))
        res_atoms = [self._group_atoms_list[i] for i in found_indices]
        return set(res_atoms)

    def _apply_distmat(self, group):
        sel_atoms = self.sel._apply(group)  # group is wrong, should be universe (?!)
        sel_CoG = AtomGroup(sel_atoms).center_of_geometry()
        sys_atoms_list = [a for a in (self._group_atoms)]  # list needed for back-indexing
        sys_ag = AtomGroup(sys_atoms_list)
        sel_CoG_str = \
            str("point ") +\
            str(sel_CoG[0]) + " " + str(sel_CoG[1]) + " " + str(sel_CoG[2]) + " " +\
            str(self.exRadius) + " and not point " + str(sel_CoG[0]) + " " + str(sel_CoG[1]) + " " + \
            str(sel_CoG[2]) + " " + str(self.inRadius)
        sel = sys_ag.select_atoms(sel_CoG_str)
        res_atoms = AtomGroup(set(sel))
        if self.periodic:
            box = group.dimensions[:3]  # ignored with KDTree
        else:
            box = None
        return set(res_atoms)

    def __repr__(self):
        return "<'SphericalLayerSelection' inner radius " + repr(self.inRadius) + ", external radius " + repr(
            self.exRadius) + " centered in " + repr(self.sel) + ">"


class SphericalZoneSelection(Selection):
    def __init__(self, sel, cutoff, periodic=None):
        Selection.__init__(self)
        self.sel = sel
        self.cutoff = cutoff
        self.sqdist = cutoff * cutoff
        if periodic is None:
            self.periodic = flags['use_periodic_selections']

    def _apply(self, group):
        # make choosing _fast/_slow configurable (while testing)
        if flags['use_KDTree_routines'] in (True, 'fast', 'always'):
            return self._apply_KDTree(group)
        else:
            return self._apply_distmat(group)

    def _apply_KDTree(self, group):
        """Selection using KDTree but periodic = True not supported.
        (KDTree routine is ca 15% slower than the distance matrix one)
        """
        sys_indices = np.array([a.index for a in self._group_atoms_list])
        sys_coor = Selection.coord[sys_indices]
        sel_atoms = self.sel._apply(group)  # group is wrong, should be universe (?!)
        sel_CoG = AtomGroup(sel_atoms).center_of_geometry()
        self.ref = np.array((sel_CoG[0], sel_CoG[1], sel_CoG[2]))
        if self.periodic:
            pass  # or warn? -- no periodic functionality with KDTree search

        kdtree = KDTree(dim=3, bucket_size=10)
        kdtree.set_coords(sys_coor)
        kdtree.search(self.ref, self.cutoff)
        found_indices = kdtree.get_indices()
        res_atoms = [self._group_atoms_list[i] for i in found_indices]
        return set(res_atoms)

    def _apply_distmat(self, group):
        sel_atoms = self.sel._apply(group)  # group is wrong, should be universe (?!)
        sel_CoG = AtomGroup(sel_atoms).center_of_geometry()
        sys_atoms_list = [a for a in (self._group_atoms)]  # list needed for back-indexing
        sys_ag = AtomGroup(sys_atoms_list)
        sel_CoG_str = str("point ") + str(sel_CoG[0]) + " " + str(sel_CoG[1]) + " " + str(sel_CoG[2]) + " " + str(
            self.cutoff)
        sel = sys_ag.select_atoms(sel_CoG_str)
        res_atoms = AtomGroup(set(sel))
        if self.periodic:
            box = group.dimensions[:3]  # ignored with KDTree
        else:
            box = None
        return set(res_atoms)

    def __repr__(self):
        return "<'SphericalZoneSelection' radius " + repr(self.cutoff) + " centered in " + repr(self.sel) + ">"


class _CylindricalSelection(Selection):
    def __init__(self, sel, exRadius, zmax, zmin, periodic=None):
        Selection.__init__(self)
        self.sel = sel
        self.exRadius = exRadius
        self.exRadiusSq = exRadius * exRadius
        self.zmax = zmax
        self.zmin = zmin
        self.periodic = flags['use_periodic_selections']

    def _apply(self, group):
        #KDTree function not implementable
        return self._apply_distmat(group)

    def _apply_distmat(self, group):
        sel_atoms = self.sel._apply(group)
        sel_CoG = AtomGroup(sel_atoms).center_of_geometry()
        coords = AtomGroup(Selection._group_atoms_list).positions

        if self.periodic and not np.any(Selection.coord.dimensions[:3]==0):
            if not np.allclose(Selection.coord.dimensions[3:],(90.,90.,90.)):
                is_triclinic = True
                box = triclinic_vectors(Selection.coord.dimensions).diagonal()
            else:
                is_triclinic = False
                box = Selection.coord.dimensions[:3]

            cyl_z_hheight = (self.zmax-self.zmin)/2

            if 2*self.exRadius > box[0]:
                raise NotImplementedError("The diameter of the cylinder selection (%.3f) is larger than the unit cell's x dimension (%.3f). Can only do selections where it is smaller or equal." % (2*self.exRadius, box[0]))
            if 2*self.exRadius > box[1]:
                raise NotImplementedError("The diameter of the cylinder selection (%.3f) is larger than the unit cell's y dimension (%.3f). Can only do selections where it is smaller or equal." % (2*self.exRadius, box[1]))
            if 2*cyl_z_hheight > box[2]:
                raise NotImplementedError("The total length of the cylinder selection in z (%.3f) is larger than the unit cell's z dimension (%.3f). Can only do selections where it is smaller or equal." % (2*cyl_z_hheight, box[2]))
            #how off-center in z is our CoG relative to the cylinder's center
            cyl_center = sel_CoG + [0,0,(self.zmax+self.zmin)/2]
            coords += box/2 - cyl_center
            coords = distances.apply_PBC(coords, box=Selection.coord.dimensions)
            if is_triclinic:
                coords = distances.apply_PBC(coords, box=box)
            sel_CoG = box/2
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
            circ_sel = (xy_norms <= self.exRadiusSq) * (xy_norms >= self.inRadiusSq)
        except AttributeError:
            circ_sel = (xy_norms <= self.exRadiusSq)
        mask_sel[mask_ndxs] = circ_sel
        ndxs = np.where(mask_sel)[0]
        res_atoms = set(Selection._group_atoms_list[ndx] for ndx in ndxs)
        return res_atoms

    def __repr__(self):
        return "<'CylindricalSelection' radius " + repr(self.exRadius) + ", zmax " + repr(
            self.zmax) + ", zmin " + repr(self.zmin) + ">"

class CylindricalZoneSelection(_CylindricalSelection):
    def __init__(self, sel, exRadius, zmax, zmin, periodic=None):
        Selection.__init__(self)
        _CylindricalSelection.__init__(self, sel, exRadius, zmax, zmin, periodic)

    def __repr__(self):
        return "<'CylindricalZoneSelection' radius " + repr(self.exRadius) + ", zmax " + repr(
            self.zmax) + ", zmin " + repr(self.zmin) + ">"

class CylindricalLayerSelection(_CylindricalSelection):
    def __init__(self, sel, inRadius, exRadius, zmax, zmin, periodic=None):
        Selection.__init__(self)
        _CylindricalSelection.__init__(self, sel, exRadius, zmax, zmin, periodic)
        self.inRadius = inRadius
        self.inRadiusSq = inRadius * inRadius

    def __repr__(self):
        return "<'CylindricalLayerSelection' inner radius " + repr(self.inRadius) + ", external radius " + repr(
            self.exRadius) + ", zmax " + repr(self.zmax) + ", zmin " + repr(self.zmin) + ">"


class PointSelection(Selection):
    def __init__(self, x, y, z, cutoff, periodic=None):
        Selection.__init__(self)
        self.ref = np.array((float(x), float(y), float(z)))
        self.cutoff = float(cutoff)
        self.cutoffsq = float(cutoff) * float(cutoff)
        if periodic is None:
            self.periodic = flags['use_periodic_selections']

    def _apply(self, group):
        # make choosing _fast/_slow configurable (while testing)
        if flags['use_KDTree_routines'] in ('always',):
            return self._apply_KDTree(group)
        else:
            return self._apply_distmat(group)

    def _apply_KDTree(self, group):
        """Selection using KDTree but periodic = True not supported.
        (KDTree routine is ca 15% slower than the distance matrix one)
        """
        sys_indices = np.array([a.index for a in self._group_atoms_list])
        sys_coor = Selection.coord[sys_indices]
        if self.periodic:
            pass  # or warn? -- no periodic functionality with KDTree search

        CNS = CoordinateNeighborSearch(sys_coor)  # cache the KDTree for this selection/frame?
        found_indices = CNS.search(self.ref, self.cutoff)
        res_atoms = [self._group_atoms_list[i] for i in found_indices]  # make list np array and use fancy indexing?
        return set(res_atoms)

    def _apply_distmat(self, group):
        """Selection that computes all distances."""
        sys_indices = np.array([a.index for a in self._group_atoms_list])
        sys_coor = Selection.coord[sys_indices]
        ref_coor = self.ref[np.newaxis, ...]
        # Fix: Arrarys need to be converted to dtype=float32 to work with distance_array
        sys_coor = np.asarray(sys_coor, dtype=np.float32)
        ref_coor = np.asarray(ref_coor, dtype=np.float32)
        if self.periodic:
            box = group.dimensions[:3]
        else:
            box = None

        dist = distances.distance_array(sys_coor, ref_coor, box)
        res_atoms = [self._group_atoms_list[i] for i in np.any(dist <= self.cutoff, axis=1).nonzero()[0]]
        # make list np array and use fancy indexing?
        return set(res_atoms)

    def __repr__(self):
        return "<'PointSelection' " + repr(self.cutoff) + " Ang around " + repr(self.ref) + ">"


class CompositeSelection(Selection):
    def __init__(self, name=None, type=None, resname=None, resid=None, segid=None):
        Selection.__init__(self)
        self.name = name
        self.type = type
        self.resname = resname
        self.resid = resid
        self.segid = segid

    def _apply(self, group):
        res = []
        for a in group.atoms:
            add = True
            if self.name is not None and a.name != self.name:
                add = False
            if self.type is not None and a.type != self.type:
                add = False
            if self.resname is not None and a.resname != self.resname:
                add = False
            if self.resid is not None and a.resid != self.resid:
                add = False
            if self.segid is not None and a.segid != self.segid:
                add = False
            if add:
                res.append(a)
        return set(res)


class AtomSelection(Selection):
    def __init__(self, name, resid, segid):
        Selection.__init__(self)
        self.name = name
        self.resid = resid
        self.segid = segid

    def _apply(self, group):
        for a in group.atoms:
            if ((a.name == self.name) and (a.resid == self.resid) and (a.segid == self.segid)):
                return set([a])
        return set([])

    def __repr__(self):
        return "<'AtomSelection' " + repr(self.segid) + " " + repr(self.resid) + " " + repr(self.name) + " >"


class SelgroupSelection(Selection):
    def __init__(self, selgroup):
        Selection.__init__(self)
        self._grp = selgroup

    def _apply(self, group):
        common = np.intersect1d(group.atoms.indices, self._grp.atoms.indices)
        res_atoms = [i for i in self._grp if i.index in common]
        return set(res_atoms)

    def __repr__(self):
        return "<" + repr(self.__class__.__name__) + ">"


@deprecate(old_name='fullgroup', new_name='global group')
class FullSelgroupSelection(Selection):
    def __init__(self, selgroup):
        Selection.__init__(self)
        self._grp = selgroup

    def _apply(self, group):
        return set(self._grp._atoms)

    def __repr__(self):
        return "<" + repr(self.__class__.__name__) + ">"


class StringSelection(Selection):
    def __init__(self, field):
        Selection.__init__(self)
        self._field = field

    def _apply(self, group):
        # Look for a wildcard
        value = getattr(self, self._field)
        wc_pos = value.find('*')  # This returns -1, so if it's not in value then use the whole of value
        if wc_pos == -1:
            wc_pos = None
        return set([a for a in group.atoms if getattr(a, self._field)[:wc_pos] == value[:wc_pos]])

    def __repr__(self):
        return "<" + repr(self.__class__.__name__) + ": " + repr(getattr(self, self._field)) + ">"


class AtomNameSelection(StringSelection):
    def __init__(self, name):
        StringSelection.__init__(self, "name")
        self.name = name


class AtomTypeSelection(StringSelection):
    def __init__(self, type):
        StringSelection.__init__(self, "type")
        self.type = type


class ResidueNameSelection(StringSelection):
    def __init__(self, resname):
        StringSelection.__init__(self, "resname")
        self.resname = resname


class SegmentNameSelection(StringSelection):
    def __init__(self, segid):
        StringSelection.__init__(self, "segid")
        self.segid = segid


class AltlocSelection(StringSelection):
    def __init__(self, altLoc):
        StringSelection.__init__(self, "altLoc")
        self.altLoc = altLoc


class ByResSelection(Selection):
    def __init__(self, sel):
        Selection.__init__(self)
        self.sel = sel

    def _apply(self, group):
        res = self.sel._apply(group)
        unique_res = set([(a.resid, a.segid) for a in res])
        sel = []
        for atom in group.atoms:
            if (atom.resid, atom.segid) in unique_res:
                sel.append(atom)
        return set(sel)

    def __repr__(self):
        return "<'ByResSelection'>"


class _RangeSelection(Selection):
    def __init__(self, lower, upper):
        Selection.__init__(self)
        self.lower = lower
        self.upper = upper

    def __repr__(self):
        return "<'" + self.__class__.__name__ + "' " + repr(self.lower) + ":" + repr(self.upper) + " >"


class ResidueIDSelection(_RangeSelection):
    def _apply(self, group):
        if self.upper is not None:
            return set([a for a in group.atoms if (self.lower <= a.resid <= self.upper)])
        else:
            return set([a for a in group.atoms if a.resid == self.lower])


class ResnumSelection(_RangeSelection):
    def _apply(self, group):
        if self.upper is not None:
            return set([a for a in group.atoms if (self.lower <= a.resnum <= self.upper)])
        else:
            return set([a for a in group.atoms if a.resnum == self.lower])


class ByNumSelection(_RangeSelection):
    def _apply(self, group):
        if self.upper is not None:
            # In this case we'll use 1 indexing since that's what the user will be
            # familiar with
            return set([a for a in group.atoms if (self.lower <= (a.index+1) <= self.upper)])
        else:
            return set([a for a in group.atoms if (a.index+1) == self.lower])


class ProteinSelection(Selection):
    """A protein selection consists of all residues with  recognized residue names.

    Recognized residue names in :attr:`ProteinSelection.prot_res`.

      * from the CHARMM force field::
         awk '/RESI/ {printf "'"'"%s"'"',",$2 }' top_all27_prot_lipid.rtf

      * manually added special CHARMM, OPLS/AA and Amber residue names.

      * still missing: Amber N- and C-terminal residue names

    .. SeeAlso:: :func:`MDAnalysis.lib.util.convert_aa_code`
    """
    #: Dictionary of recognized residue names (3- or 4-letter).
    prot_res = dict([(x, None) for x in [
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
    ]])

    def _apply(self, group):
        return set([a for a in group.atoms if a.resname in self.prot_res])

    def __repr__(self):
        return "<'ProteinSelection' >"


class NucleicSelection(Selection):
    """A nucleic selection consists of all atoms in nucleic acid residues with  recognized residue names.

    Recognized residue names:

    * from the CHARMM force field ::
        awk '/RESI/ {printf "'"'"%s"'"',",$2 }' top_all27_prot_na.rtf
    * recognized: 'ADE', 'URA', 'CYT', 'GUA', 'THY'
    * recognized (CHARMM in Gromacs): 'DA', 'DU', 'DC', 'DG', 'DT'

    .. versionchanged:: 0.8
       additional Gromacs selections (see also :class:`NucleicXstalSelection`)
    """
    nucl_res = dict([(x, None) for x in [
        'ADE', 'URA', 'CYT', 'GUA', 'THY', 'DA', 'DC', 'DG', 'DT', 'RA',
        'RU', 'RG', 'RC', 'A', 'T', 'U', 'C', 'G']])

    def _apply(self, group):
        return set([a for a in group.atoms if a.resname in self.nucl_res])

    def __repr__(self):
        return "<'NucleicSelection' >"


class BackboneSelection(ProteinSelection):
    """A BackboneSelection contains all atoms with name 'N', 'CA', 'C', 'O'.

    This excludes OT* on C-termini (which are included by, eg VMD's backbone selection).
    """
    bb_atoms = dict([(x, None) for x in ['N', 'CA', 'C', 'O']])

    def _apply(self, group):
        return set([a for a in group.atoms if (a.name in self.bb_atoms and a.resname in self.prot_res)])

    def __repr__(self):
        return "<'BackboneSelection' >"


class NucleicBackboneSelection(NucleicSelection):
    """A NucleicBackboneSelection contains all atoms with name "P", "C5'", C3'", "O3'", "O5'".

    These atoms are only recognized if they are in a residue matched
    by the :class:`NucleicSelection`.
    """
    bb_atoms = dict([(x, None) for x in ["P", "C5'", "C3'", "O3'", "O5'"]])

    def _apply(self, group):
        return set([a for a in group.atoms if (a.name in self.bb_atoms and a.resname in self.nucl_res)])

    def __repr__(self):
        return "<'NucleicBackboneSelection' >"


class BaseSelection(NucleicSelection):
    """Selection of atoms in nucleobases.

    Recognized atom names (from CHARMM):

     'N9', 'N7', 'C8', 'C5', 'C4', 'N3', 'C2', 'N1', 'C6',
     'O6','N2','N6', 'O2','N4','O4','C5M'
    """
    base_atoms = dict([(x, None) for x in [
        'N9', 'N7', 'C8', 'C5', 'C4', 'N3', 'C2', 'N1', 'C6',
        'O6', 'N2', 'N6',
        'O2', 'N4', 'O4', 'C5M']])

    def _apply(self, group):
        return set([a for a in group.atoms if (a.name in self.base_atoms and a.resname in self.nucl_res)])

    def __repr__(self):
        return "<'BaseSelection' >"


class NucleicSugarSelection(NucleicSelection):
    """A NucleicSugarSelection contains all atoms with name C1', C2', C3', C4', O2', O4', O3'.
    """
    sug_atoms = dict([(x, None) for x in ['C1\'', 'C2\'', 'C3\'', 'C4\'', 'O4\'']])

    def _apply(self, group):
        return set([a for a in group.atoms if (a.name in self.sug_atoms and a.resname in self.nucl_res)])

    def __repr__(self):
        return "<'NucleicSugarSelection' >"


class CASelection(BackboneSelection):
    """Select atoms named CA in protein residues (supposed to be the C-alphas)
    """

    def _apply(self, group):
        return set([a for a in group.atoms if (a.name == "CA" and a.resname in self.prot_res)])

    def __repr__(self):
        return "<'CASelection' >"


class BondedSelection(Selection):
    def __init__(self, sel):
        Selection.__init__(self)
        self.sel = sel

    def _apply(self, group):
        res = self.sel._apply(group)
        # Find all the atoms bonded to each
        sel = []
        for atom in res:
            for b1, b2 in group._bonds:
                if atom.index == b1:
                    sel.append(group.atoms[b2])
                elif atom.index == b2:
                    sel.append(group.atoms[b1])
        return set(sel)

    def __repr__(self):
        return "<'BondedSelection' to " + repr(self.sel) + " >"


class PropertySelection(Selection):
    """Some of the possible properties:
    x, y, z, radius, mass,
    """

    def __init__(self, prop, operator, value, abs=False):
        Selection.__init__(self)
        self.prop = prop
        self.operator = operator
        self.value = value
        self.abs = abs

    def _apply(self, group):
        # For efficiency, get a reference to the actual np position arrays
        if self.prop in ("x", "y", "z"):
            p = getattr(Selection.coord, '_' + self.prop)
            indices = np.array([a.index for a in group.atoms])
            if not self.abs:
                # XXX Hack for difference in np.nonzero between version < 1. and version > 1
                res = np.nonzero(self.operator(p[indices], self.value))
            else:
                res = np.nonzero(self.operator(np.abs(p[indices]), self.value))
            if type(res) == tuple:
                res = res[0]
            result_set = [group.atoms[i] for i in res]
        elif self.prop == "mass":
            result_set = [a for a in group.atoms if self.operator(a.mass, self.value)]
        elif self.prop == "charge":
            result_set = [a for a in group.atoms if self.operator(a.charge, self.value)]
        return set(result_set)

    def __repr__(self):
        if self.abs:
            abs_str = " abs "
        else:
            abs_str = ""
        return "<'PropertySelection' " + abs_str + repr(self.prop) + " " + repr(self.operator.__name__) + " " + repr(
            self.value) + ">"

class SameSelection(Selection):
    # When adding new keywords here don't forget to also add them to the
    #  case statement under the SAME op, where they are first checked.
    def __init__(self, sel, prop):
        Selection.__init__(self)
        self.sel = sel
        self.prop = prop
    def _apply(self, group):
        res = self.sel._apply(group)
        if not res:
            return set([])
        if self.prop in ("residue", "fragment", "segment"):
            atoms = set([])
            for a in res:
                if a not in atoms:
                # This shortcut assumes consistency that all atoms in a residue/fragment/segment
                # belong to the exact same residue/fragment/segment.
                    atoms |= set(getattr(a, self.prop).atoms)
            return Selection._group_atoms & atoms
        elif self.prop in ("name", "type", "resname", "resid", "segid", "mass", "charge", "radius", "bfactor", "resnum"):
            props = [getattr(a, self.prop) for a in res]
            result_set = (a for a in Selection._group_atoms if getattr(a, self.prop) in props)
        elif self.prop in ("x", "y", "z"):
            p = getattr(Selection.coord, "_"+self.prop)
            res_indices = np.array([a.index for a in res])
            sel_indices = np.array([a.index for a in Selection._group_atoms])
            result_set = group.atoms[np.where(np.in1d(p[sel_indices], p[res_indices]))[0]]._atoms
        else:
            self.__error(self.prop, expected=False)
        return set(result_set)
    def __repr__(self):
        return "<'SameSelection' of "+ repr(self.prop)+" >"


class ParseError(Exception):
    pass


class SelectionParser:
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
    BYRES = 'byres'
    BONDED = 'bonded'
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
    NUCLEICXSTAL = 'nucleicxstal'
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
        (ALL, AllSelection), (NOT, NotSelection), (AND, AndSelection), (OR, OrSelection),
        (SEGID, SegmentNameSelection), (RESID, ResidueIDSelection), (RESNUM, ResnumSelection),
        (RESNAME, ResidueNameSelection), (NAME, AtomNameSelection), (ALTLOC, AltlocSelection),
        (TYPE, AtomTypeSelection), (BYRES, ByResSelection),
        (BYNUM, ByNumSelection), (PROP, PropertySelection),
        (AROUND, AroundSelection), (SPHLAYER, SphericalLayerSelection), (SPHZONE, SphericalZoneSelection),
        (CYLAYER, CylindricalLayerSelection), (CYZONE, CylindricalZoneSelection),
        (POINT, PointSelection), (NUCLEIC, NucleicSelection), (PROTEIN, ProteinSelection),
        (BB, BackboneSelection), (NBB, NucleicBackboneSelection),
        (BASE, BaseSelection), (SUGAR, NucleicSugarSelection),
        #(BONDED, BondedSelection), not supported yet, need a better way to walk the bond lists
        (ATOM, AtomSelection), (SELGROUP, SelgroupSelection), (FULLSELGROUP, FullSelgroupSelection),
        (SAME, SameSelection), (GLOBAL, GlobalSelection)])
    associativity = dict([(AND, "left"), (OR, "left")])
    precedence = dict(
        [(AROUND, 1), (SPHLAYER, 1), (SPHZONE, 1), (CYLAYER, 1), (CYZONE, 1), (POINT, 1), (BYRES, 1), (BONDED, 1),
            (SAME, 1), (AND, 3), (OR, 3), (NOT, 5), (GLOBAL, 5)])

    # Borg pattern: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66531
    _shared_state = {}

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __peek_token(self):
        """Looks at the next token in our token stream."""
        return self.tokens[0]

    def __consume_token(self):
        """Pops off the next token in our token stream."""
        return self.tokens.pop(0)

    def __error(self, token, expected=True):
        """Stops parsing and reports an error."""
        if expected:
            raise ParseError("Parsing error- '" + self.selectstr + "'\n" + repr(token) + " expected")
        else:
            raise ParseError("Parsing error- '" + self.selectstr + "'\n" + repr(token) + " unexpected")

    def __expect(self, token):
        if self.__peek_token() == token:
            self.__consume_token()
        else:
            self.__error(token)

    def parse(self, selectstr, selgroups):
        self.selectstr = selectstr
        self.selgroups = selgroups
        self.tokens = selectstr.replace('(', ' ( ').replace(')', ' ) ').split() + [self.EOF]
        parsetree = self.__parse_expression(0)
        self.__expect(self.EOF)
        return parsetree

    def __parse_expression(self, p):
        exp1 = self.__parse_subexp()
        while self.__peek_token() in (self.AND, self.OR) and self.precedence[self.__peek_token()] >= p:  # bin ops
            op = self.__consume_token()
            if self.associativity[op] == "right":
                q = self.precedence[op]
            else:
                q = 1 + self.precedence[op]
            exp2 = self.__parse_expression(q)
            exp1 = self.classdict[op](exp1, exp2)
        return exp1

    def __parse_subexp(self):
        op = self.__consume_token()
        if op in (self.NOT, self.BYRES, self.GLOBAL):  # unary operators
            exp = self.__parse_expression(self.precedence[op])
            return self.classdict[op](exp)
        elif op in (self.AROUND):
            dist = self.__consume_token()
            exp = self.__parse_expression(self.precedence[op])
            return self.classdict[op](exp, float(dist))
        elif op in (self.SPHLAYER):
            inRadius = self.__consume_token()
            exRadius = self.__consume_token()
            exp = self.__parse_expression(self.precedence[op])
            return self.classdict[op](exp, float(inRadius), float(exRadius))
        elif op in (self.SPHZONE):
            dist = self.__consume_token()
            exp = self.__parse_expression(self.precedence[op])
            return self.classdict[op](exp, float(dist))
        elif op in (self.CYLAYER):
            inRadius = self.__consume_token()
            exRadius = self.__consume_token()
            zmax = self.__consume_token()
            zmin = self.__consume_token()
            exp = self.__parse_expression(self.precedence[op])
            return self.classdict[op](exp, float(inRadius), float(exRadius), float(zmax), float(zmin))
        elif op in (self.CYZONE):
            exRadius = self.__consume_token()
            zmax = self.__consume_token()
            zmin = self.__consume_token()
            exp = self.__parse_expression(self.precedence[op])
            return self.classdict[op](exp, float(exRadius), float(zmax), float(zmin))
        elif op in (self.POINT):
            dist = self.__consume_token()
            x = self.__consume_token()
            y = self.__consume_token()
            z = self.__consume_token()
            return self.classdict[op](float(dist), float(x), float(y), float(z))
        elif op == self.BONDED:
            exp = self.__parse_expression(self.precedence[op])
            return self.classdict[op](exp)
        elif op == self.LPAREN:
            exp = self.__parse_expression(0)
            self.__expect(self.RPAREN)
            return exp
        elif op in (self.SEGID, self.RESNAME, self.NAME, self.TYPE, self.ALTLOC):
            data = self.__consume_token()
            if data in (
                    self.LPAREN,
                    self.RPAREN, self.AND, self.OR, self.NOT, self.SEGID, self.RESID, self.RESNAME, self.NAME,
                    self.TYPE):
                self.__error("Identifier")
            return self.classdict[op](data)
        elif op in (self.SELGROUP, self.FULLSELGROUP):
            grpname = self.__consume_token()
            return self.classdict[op](self.selgroups[grpname])
        elif op == self.PROTEIN:
            return self.classdict[op]()
        elif op == self.NUCLEIC:
            return self.classdict[op]()
        elif op == self.NUCLEICXSTAL:
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
        elif op in (self.RESID, self.RESNUM, self.BYNUM):  # can operate on ranges X:Y or X-Y
            data = self.__consume_token()
            try:
                lower = int(data)
                upper = None
            except ValueError:
                selrange = re.match("(\d+)[:-](\d+)",
                                    data)  # check if in appropriate format 'lower:upper' or 'lower-upper'
                if not selrange:
                    self.__error(op)
                lower, upper = map(int, selrange.groups())
            return self.classdict[op](lower, upper)
        elif op == self.PROP:
            prop = self.__consume_token()
            if prop == "abs":
                abs = True
                prop = self.__consume_token()
            else:
                abs = False
            oper = self.__consume_token()
            value = float(self.__consume_token())
            ops = dict([
                (self.GT, np.greater), (self.LT, np.less),
                (self.GE, np.greater_equal), (self.LE, np.less_equal),
                (self.EQ, np.equal), (self.NE, np.not_equal)])
            if oper in ops.keys():
                return self.classdict[op](prop, ops[oper], value, abs)
        elif op == self.ATOM:
            segid = self.__consume_token()
            resid = int(self.__consume_token())
            name = self.__consume_token()
            return self.classdict[op](name, resid, segid)
        elif op == self.SAME:
            prop = self.__consume_token()
            self.__expect("as")
            if prop in ("name", "type", "resname", "resid", "segid", "mass", "charge", "radius", "bfactor",
                        "resnum", "residue", "segment", "fragment", "x", "y", "z"):
                exp = self.__parse_expression(self.precedence[op])
                return self.classdict[op](exp, prop)
            else:
                self.__error(prop, expected=False)
        else:
            self.__error(op, expected=False)

# The module level instance
Parser = SelectionParser()
