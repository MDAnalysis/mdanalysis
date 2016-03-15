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
Topology attribute objects --- :mod:`MDAnalysis.core.topologyattrs'
===================================================================

Common TopologyAttrs used by most topology parsers.

"""
from collections import defaultdict
import itertools
import numpy as np

from . import flags
from ..lib.util import cached
from ..exceptions import NoDataError
from .topologyobjects import TopologyGroup
from . import selection
from . import flags


class TopologyAttr(object):
    """Base class for Topology attributes.

    .. note::   This class is intended to be subclassed, and mostly amounts to a
                skeleton. The methods here should be present in all
                :class:`TopologyAttr` child classes, but by default they raise
                appropriate exceptions.

    Attributes
    ----------
    attrname : str
        the name used for the attribute when attached to a ``Topology`` object
    top : Topology
        handle for the Topology object TopologyAttr is associated with
        
    """
    attrname = 'topologyattrs'
    singular = 'topologyattr'
    top = None

    groupdoc = None
    singledoc = None

    def __init__(self, values):
        self.values = values

    def __len__(self):
        """Length of the TopologyAttr at its intrinsic level."""
        return len(self.values)

    def __getitem__(self, group):
        """Accepts an AtomGroup, ResidueGroup or SegmentGroup"""
        if group.level == 'atom':
            return self.get_atoms(group)
        elif group.level == 'residue':
            return self.get_residues(group)
        elif group.level == 'segment':
            return self.get_segments(group)

    def __setitem__(self, group, values):
        if group.level == 'atom':
            return self.set_atoms(group, values)
        elif group.level == 'residue':
            return self.set_residues(group, values)
        elif group.level == 'segment':
            return self.set_segments(group, values)

    def get_atoms(self, ag):
        """Get atom attributes for a given AtomGroup"""
        # aix = ag.indices
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


## core attributes

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
    target_levels = ['atom']

    def __init__(self):
        pass

    def __len__(self):
        """Length of the TopologyAttr at its intrinsic level."""
        return len(self.top.n_atoms)

    def set_atoms(self, ag, values):
        raise AttributeError("Atom indices are fixed; they cannot be reset")

    def get_atoms(self, ag):
        return ag._ix

    def get_residues(self, rg):
        return self.top.tt.residues2atoms_1d(rg._ix)

    def get_segments(self, sg):
        return self.top.tt.segments2atoms_1d(sg._ix)


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
    target_levels = ['atom', 'residue']

    def __init__(self):
        pass

    def __len__(self):
        """Length of the TopologyAttr at its intrinsic level."""
        return len(self.top.n_residues)

    def get_atoms(self, ag):
        return self.top.tt.atoms2residues(ag._ix)

    def set_atoms(self, ag, values):
        """Set resindex for each atom given. Effectively moves each atom to
        another residue.

        """
        self.top.tt.move_atom(ag._ix, values)

    def get_residues(self, rg):
        return rg._ix

    def set_residues(self, rg, values):
        raise AttributeError("Residue indices are fixed; they cannot be reset")

    def get_segments(self, sg):
        return self.top.tt.segments2residues_1d(sg._ix)


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
    target_levels = ['atom', 'residue', 'segment']

    def __init__(self):
        pass

    def __len__(self):
        """Length of the TopologyAttr at its intrinsic level."""
        return len(self.top.n_segments)

    def get_atoms(self, ag):
        return self.top.tt.atoms2segments(ag._ix)

    def get_residues(self, rg):
        return self.top.tt.residues2segments(rg._ix)

    def get_segments(self, sg):
        return sg._ix

    def set_segments(self, sg, values):
        raise AttributeError("Segment indices are fixed; they cannot be reset")


## atom attributes

class AtomAttr(TopologyAttr):
    """Base class for atom attributes.

    """
    attrname = 'atomattrs'
    singular = 'atomattr'
    target_levels = ['atom']

    def get_atoms(self, ag):
        return self.values[ag._ix]

    def set_atoms(self, ag, values):
        self.values[ag._ix] = values

    def get_residues(self, rg):
        """By default, the values for each atom present in the set of residues
        are returned in a single array. This behavior can be overriden in child
        attributes.

        """
        aix = self.top.tt.residues2atoms_1d(rg._ix)
        return self.values[aix]

    def get_segments(self, sg):
        """By default, the values for each atom present in the set of residues
        are returned in a single array. This behavior can be overriden in child
        attributes.

        """
        aix = self.top.tt.segments2atoms_1d(sg._ix)
        return self.values[aix]


#TODO: update docs to property doc
class Atomids(AtomAttr):
    """ID for each atom.
    """
    attrname = 'ids'
    singular = 'id'


#TODO: update docs to property doc
class Atomnames(AtomAttr):
    """Name for each atom.
    """
    attrname = 'names'
    singular = 'name'
    transplants = defaultdict(list)

    def getattr__(atomgroup, name):
        try:
            return atomgroup._get_named_atom(name)
        except selection.SelectionError:
            raise AttributeError("'{0}' object has no attribute '{1}'".format(
                    atomgroup.__class__.__name__, name))

    transplants['atomgroup'].append(
        ('__getattr__', getattr__))

    transplants['residue'].append(
        ('__getattr__', getattr__))

    def _get_named_atom(group, name):
        """Get all atoms with name *name* in the current AtomGroup.

        For more than one atom it returns a list of :class:`Atom`
        instance. A single :class:`Atom` is returned just as such. If
        no atoms are found, a :exc:`SelectionError` is raised.

        .. versionadded:: 0.9.2
        """
        # There can be more than one atom with the same name
        atomlist = group.atoms.unique[group.atoms.unique.names == name]
        if len(atomlist) == 0:
            raise selection.SelectionError(
                "No atoms with name '{0}'".format(name))
        elif len(atomlist) == 1:
            # XXX: keep this, makes more sense for names
            return atomlist[0]
        else:
            # XXX: but inconsistent (see residues and Issue 47)
            return atomlist

    transplants['atomgroup'].append(
        ('_get_named_atom', _get_named_atom))

    transplants['residue'].append(
        ('_get_named_atom', _get_named_atom))


#TODO: update docs to property doc
class Atomtypes(AtomAttr):
    """Type for each atom"""
    attrname = 'types'
    singular = 'type'


#TODO: update docs to property doc
class Radii(AtomAttr):
    """Radii for each atom"""
    attrname = 'radii'
    singular = 'radius'


class ChainIDs(AtomAttr):
    """ChainID per atom

    Note
    ----
    This is an attribute of the Atom, not Residue or Segment
    """
    attrname = 'chainIDs'
    singular = 'chainID'


class ICodes(AtomAttr):
    """Insertion code for Atoms"""
    attrname = 'icodes'
    singular = 'icode'


class Tempfactors(AtomAttr):
    """Tempfactor for atoms"""
    attrname = 'tempfactors'
    singular = 'tempfactor'


#TODO: need to add cacheing
class Masses(AtomAttr):
    attrname = 'masses'
    singular = 'mass'
    target_levels = ['atom', 'residue', 'segment']
    transplants = defaultdict(list)

    groupdoc = """Mass of each component in the Group.

    If the Group is an AtomGroup, then the masses are for each atom. If the
    Group is a ResidueGroup or SegmentGroup, the masses are for each residue or
    segment, respectively. These are obtained by summation of the member atoms
    for each component.
    """

    singledoc = """Mass of the component."""

    def get_residues(self, rg):
        resatoms = self.top.tt.residues2atoms_2d(rg._ix)

        if isinstance(rg._ix, int):
            # for a single residue
            masses = self.values[resatoms].sum()
        else:
            # for a residuegroup
            masses = np.empty(len(rg))
            for i, row in enumerate(resatoms):
                masses[i] = self.values[row].sum()

        return masses

    def get_segments(self, sg):
        segatoms = self.top.tt.segments2atoms_2d(sg._ix)

        if isinstance(sg._ix, int):
            # for a single segment
            masses = self.values[segatoms].sum()
        else:
            # for a segmentgroup
            masses = np.empty(len(sg))
            for i, row in enumerate(segatoms):
                masses[i] = self.values[row].sum()

        return masses

    def center_of_mass(atomgroup, **kwargs):
        """Center of mass of the AtomGroup.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]

        .. note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
            ``True`` allows the *pbc* flag to be used by default.

        .. versionchanged:: 0.8 Added `pbc` parameter
        """
        pbc = kwargs.pop('pbc', flags['use_pbc'])
        masses = atomgroup.masses

        if pbc:
            positions = atomgroup.pack_into_box(inplace=False)
        else:
            positions = atomgroup.positions

        return np.sum(positions * masses[:, np.newaxis],
                      axis=0) / masses.sum()

    transplants['atomgroup'].append(
        ('center_of_mass', center_of_mass))

    def total_mass(group):
        """Total mass of the Group.
        
        """
        return group.masses.sum()

    transplants['group'].append(
        ('total_mass', total_mass))

    def moment_of_inertia(atomgroup, **kwargs):
        """Tensor moment of inertia relative to center of mass as 3x3 numpy
        array.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]

        .. note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
            ``True`` allows the *pbc* flag to be used by default.

        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', flags['use_pbc'])

        # Convert to local coordinates
        if pbc:
            pos = atomgroup.pack_into_box(inplace=False) - atomgroup.center_of_mass(pbc=True)
        else:
            pos = atomgroup.positions - atomgroup.center_of_mass(pbc=False)

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
        tens = np.zeros((3,3), dtype=np.float64)
        # xx
        tens[0][0] = (masses * (pos[:,1] * pos[:,1] + pos[:,2] * pos[:,2])).sum()
        # xy & yx
        tens[0][1] = tens[1][0] = - (masses * pos[:,0] * pos[:,1]).sum()
        # xz & zx
        tens[0][2] = tens[2][0] = - (masses * pos[:,0] * pos[:,2]).sum()
        # yy
        tens[1][1] = (masses * (pos[:,0] * pos[:,0] + pos[:,2] * pos[:,2])).sum()
        # yz + zy
        tens[1][2] = tens[2][1] = - (masses * pos[:,1] * pos[:,2]).sum()
        # zz
        tens[2][2] = (masses * (pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1])).sum()

        return tens

    transplants['atomgroup'].append(
        ('moment_of_inertia', moment_of_inertia))

    def radius_of_gyration(atomgroup, **kwargs):
        """Radius of gyration.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]

        .. note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
            ``True`` allows the *pbc* flag to be used by default.

        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', flags['use_pbc'])
        masses = atomgroup.masses
        if pbc:
            recenteredpos = atomgroup.pack_into_box(inplace=False) - atomgroup.center_of_mass(pbc=True)
        else:
            recenteredpos = atomgroup.positions - atomgroup.center_of_mass(pbc=False)
        rog_sq = np.sum(masses * np.sum(np.power(recenteredpos, 2), axis=1)) / atomgroup.total_mass()
        return np.sqrt(rog_sq)

    transplants['atomgroup'].append(
        ('radius_of_gyration', radius_of_gyration))

    def shape_parameter(atomgroup, **kwargs):
        """Shape parameter.

        See [Dima2004]_ for background information.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]

        .. note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
            ``True`` allows the *pbc* flag to be used by default.

        .. [Dima2004] Dima, R. I., & Thirumalai, D. (2004). Asymmetry in the
                  shapes of folded and denatured states of proteins. *J
                  Phys Chem B*, 108(21),
                  6564-6570. doi:`10.1021/jp037128y`_

        .. versionadded:: 0.7.7
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', flags['use_pbc'])
        masses = atomgroup.masses
        if pbc:
            recenteredpos = atomgroup.pack_into_box(inplace=False) - atomgroup.center_of_mass(pbc=True)
        else:
            recenteredpos = atomgroup.positions - atomgroup.center_of_mass(pbc=False)
        tensor = np.zeros((3, 3))
        for x in xrange(recenteredpos.shape[0]):
            tensor += masses[x] * np.outer(recenteredpos[x, :],
                                              recenteredpos[x, :])
        tensor /= atomgroup.total_mass()
        eig_vals = np.linalg.eigvalsh(tensor)
        shape = 27.0 * np.prod(eig_vals - np.mean(eig_vals)) / np.power(np.sum(eig_vals), 3)
        return shape

    transplants['atomgroup'].append(
        ('shape_parameter', shape_parameter))

    def asphericity(atomgroup, **kwargs):
        """Asphericity.

        See [Dima2004]_ for background information.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]

        .. note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
            ``True`` allows the *pbc* flag to be used by default.

        .. [Dima2004] Dima, R. I., & Thirumalai, D. (2004). Asymmetry in the
                  shapes of folded and denatured states of proteins. *J
                  Phys Chem B*, 108(21),
                  6564-6570. doi:`10.1021/jp037128y`_

        .. versionadded:: 0.7.7
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        masses = atomgroup.masses
        if pbc:
            recenteredpos = atomgroup.pack_into_box(inplace=False) - atomgroup.center_of_mass(pbc=True)
        else:
            recenteredpos = atomgroup.positions - atomgroup.center_of_mass(pbc=False)
        tensor = np.zeros((3, 3))
        for x in xrange(recenteredpos.shape[0]):
            tensor += masses[x] * np.outer(recenteredpos[x, :],
                                              recenteredpos[x, :])
        tensor /= atomgroup.total_mass()
        eig_vals = np.linalg.eigvalsh(tensor)
        shape = (3.0 / 2.0) * np.sum(np.power(eig_vals - np.mean(eig_vals), 2)) / np.power(
            np.sum(eig_vals), 2)
        return shape

    transplants['atomgroup'].append(
        ('asphericity', asphericity))

    def principal_axes(atomgroup, **kwargs):
        """Calculate the principal axes from the moment of inertia.

        e1,e2,e3 = AtomGroup.principal_axes()

        The eigenvectors are sorted by eigenvalue, i.e. the first one
        corresponds to the highest eigenvalue and is thus the first principal axes.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]

        .. note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
            ``True`` allows the *pbc* flag to be used by default.

        Returns
        -------
        axis_vectors : array
            3 x 3 array with ``v[0]`` as first, ``v[1]`` as second, and
            ``v[2]`` as third eigenvector.

        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        if pbc:
            eigenval, eigenvec = eig(atomgroup.moment_of_inertia(pbc=True))
        else:
            eigenval, eigenvec = eig(atomgroup.moment_of_inertia(pbc=False))
        # Sort
        indices = np.argsort(eigenval)
        # Return transposed in more logical form. See Issue 33.
        return eigenvec[:, indices].T


#TODO: need to add cacheing
#TODO: update docs to property doc
class Charges(AtomAttr):
    attrname = 'charges'
    singular = 'charge'
    target_levels = ['atom', 'residue', 'segment']
    transplants = defaultdict(list)

    def get_residues(self, rg):
        charges = np.empty(len(rg))

        resatoms = self.top.tt.residues2atoms_2d(rg._ix)

        for i, row in enumerate(resatoms):
            charges[i] = self.values[row].sum()

        return charges

    def get_segments(self, sg):
        charges = np.empty(len(sg))

        segatoms = self.top.tt.segments2atoms_2d(sg._ix)

        for i, row in enumerate(segatoms):
            charges[i] = self.values[row].sum()

        return charges

    def total_charge(group):
        """Total charge of the Group.
        
        """
        return group.charges.sum()

    transplants['group'].append(
        ('total_charge', total_charge))


#TODO: update docs to property doc
class Bfactors(AtomAttr):
    """Crystallographic B-factors in A**2 for each atom"""
    attrname = 'bfactors'
    singular = 'bfactor'


#TODO: update docs to property doc
class Occupancies(AtomAttr):
    attrname = 'occupancies'
    singular = 'occupancy'


#TODO: update docs to property doc
class AltLocs(AtomAttr):
    """AltLocs for each atom"""
    attrname = 'altLocs'
    singular = 'altLoc'


## residue attributes

class ResidueAttr(TopologyAttr):
    """Base class for Topology attributes.

    .. note::   This class is intended to be subclassed, and mostly amounts to a
                skeleton. The methods here should be present in all
                :class:`TopologyAttr` child classes, but by default they raise
                appropriate exceptions.

    """
    attrname = 'residueattrs'
    singular = 'residueattr'
    target_levels = ['residue']

    def get_atoms(self, ag):
        rix = self.top.tt.atoms2residues(ag._ix)
        return self.values[rix]

    def get_residues(self, rg):
        return self.values[rg._ix]

    def set_residues(self, rg, values):
        self.values[rg._ix] = values

    def get_segments(self, sg):
        """By default, the values for each residue present in the set of
        segments are returned in a single array. This behavior can be overriden
        in child attributes.

        """
        rix = self.top.tt.segments2residues_1d(sg._ix)
        return self.values[rix]


#TODO: update docs to property doc
class Resids(ResidueAttr):
    """Residue ID"""
    attrname = 'resids'
    singular = 'resid'
    target_levels = ['atom', 'residue']


#TODO: update docs to property doc
class Resnames(ResidueAttr):
    attrname = 'resnames'
    singular = 'resname'
    target_levels = ['atom', 'residue']
    transplants = defaultdict(list)

    def getattr__(residuegroup, resname):
        try:
            return residuegroup._get_named_residue(resname)
        except selection.SelectionError:
            raise AttributeError("'{0}' object has no attribute '{1}'".format(
                    residuegroup.__class__.__name__, resname))

    transplants['residuegroup'].append(('__getattr__', getattr__))

    transplants['segment'].append(('__getattr__', getattr__))

    def _get_named_residue(group, resname):
        """Get all residues with name *resname* in the current ResidueGroup
        or Segment.

        For more than one residue it returns a
        :class:`MDAnalysis.core.groups.ResidueGroup` instance. A single
        :class:`MDAnalysis.core.group.Residue` is returned for a single match. If no
        residues are found, a :exc:`SelectionError` is raised.

        .. versionadded:: 0.9.2
        """
        # There can be more than one residue with the same name
        residues = group.residues.unique[
                group.residues.unique.resnames == resname]
        if len(residues) == 0:
            raise selection.SelectionError(
                "No residues with resname '{0}'".format(resname))
        elif len(residues) == 1:
            # XXX: keep this, makes more sense for names
            return residues[0]
        else:
            # XXX: but inconsistent (see residues and Issue 47)
            return residues

    transplants['residuegroup'].append(
        ('_get_named_residue', _get_named_residue))

    transplants['segment'].append(
        ('_get_named_residue', _get_named_residue))


#TODO: update docs to property doc
class Resnums(ResidueAttr):
    attrname = 'resnums'
    singular = 'resnum'
    target_levels = ['atom', 'residue']


## segment attributes

class SegmentAttr(TopologyAttr):
    """Base class for segment attributes.

    """
    attrname = 'segmentattrs'
    singular = 'segmentattr'
    target_levels = ['segment']

    def get_atoms(self, ag):
        six = self.top.tt.atoms2segments(ag._ix)
        return self.values[six]

    def get_residues(self, rg):
        six = self.top.tt.residues2segments(rg._ix)
        return self.values[six]

    def get_segments(self, sg):
        return self.values[sg._ix]

    def set_segments(self, sg, values):
        self.values[sg._ix] = values


#TODO: update docs to property doc
class Segids(SegmentAttr):
    attrname = 'segids'
    singular = 'segid'
    target_levels = ['atom', 'residue', 'segment']
    transplants = defaultdict(list)

    def getattr__(segmentgroup, segid):
        try:
            return segmentgroup._get_named_segment(segid)
        except selection.SelectionError:
            raise AttributeError("'{0}' object has no attribute '{1}'".format(
                    segmentgroup.__class__.__name__, segid))

    transplants['segmentgroup'].append(
        ('__getattr__', getattr__))

    def _get_named_segment(group, segid):
        """Get all segments with name *segid* in the current SegmentGroup.

        For more than one residue it returns a
        :class:`MDAnalysis.core.groups.SegmentGroup` instance. A single
        :class:`MDAnalysis.core.group.Segment` is returned for a single match. If no
        residues are found, a :exc:`SelectionError` is raised.

        .. versionadded:: 0.9.2
        """
        # There can be more than one segment with the same name
        segments  = group.segments.unique[
                group.segments.unique.segids == segid]
        if len(segments) == 0:
            raise selection.SelectionError(
                "No segments with segid '{0}'".format(segid))
        elif len(segments) == 1:
            # XXX: keep this, makes more sense for names
            return segments[0]
        else:
            # XXX: but inconsistent (see residues and Issue 47)
            return segment 

    transplants['segmentgroup'].append(
        ('_get_named_segment', _get_named_segment))


#TODO: update docs to property doc
class Bonds(AtomAttr):
    """Bonds for atoms"""
    attrname = 'bonds'
    # Singular is the same because one Atom might have
    # many bonds, so still asks for "bonds" in the plural
    singular = 'bonds'
    transplants = defaultdict(list)

    def __init__(self, values, types=None):
        """
        Arguments
        ---------
        values - list of tuples of indices.  Should be zero based
        to match the atom indices

        Eg:  [(0, 1), (1, 2), (2, 3)]
        """
        self.values = values
        self._cache = dict()
        self.order_bonds()
        self._bondtypes = {}
        if types is None:
            types = [None for value in values]
        for value, type in zip(self.values, types):
            self._bondtypes[value] = type

    def __len__(self):
        return len(self._bondDict)

    def order_bonds(self):
        ordered_values = []
        for b in self.values:
            if b[0] > b[-1]:
                b = b[::-1]
            ordered_values.append(b)
        self.values = ordered_values

    @property
    @cached('bd')
    def _bondDict(self):
        """Lazily built mapping of atoms:bonds"""
        bd = defaultdict(list)

        for b in self.values:
            # We always want the first index
            # to be less than the last
            # eg (0, 1) not (1, 0)
            # and (4, 10, 8) not (8, 10, 4)
            if b[0] > b[-1]:
                b = b[::-1]
            for a in b:
                bd[a].append(b)
        return bd

    def get_atoms(self, ag):
        try:
            unique_bonds =  set(itertools.chain(
                *[self._bondDict[a] for a in ag._ix]))
        except TypeError:
            # maybe we got passed an Atom
            unique_bonds = self._bondDict[ag._ix]
        if self._bondtypes is not None:
            _bondtypes = {value: self._bondtypes[value] for value in unique_bonds}
        else:
            _bondtypes = None
        bond_idx = np.array(sorted(unique_bonds))
        return TopologyGroup(bond_idx, ag._u, self.singular[:-1], np.array(_bondtypes.values()))

    def bonded_atoms(self):
        """An AtomGroup of all atoms bonded to this Atom"""
        idx = [b.partner(self).index for b in self.bonds]
        return self._u.atoms[idx]

    transplants['atom'].append(
        ('bonded_atoms', property(bonded_atoms, None, None,
                                  bonded_atoms.__doc__)))

    def fragment(self):
        """The fragment that this Atom is part of

        .. versionadded:: 0.9.0
        """
        return self.universe._fragdict[self]

    def fragments(self):
        """Read-only list of fragments.

        Contains all fragments that any Atom in this AtomGroup is
        part of, the contents of the fragments may extend beyond the
        contents of this AtomGroup.

        .. versionadded 0.9.0
        """
        return tuple(set(a.fragment for a in self))

    transplants['atom'].append(
        ('fragment', property(fragment, None, None,
                              fragment.__doc__)))

    transplants['atomgroup'].append(
        ('fragments', property(fragments, None, None,
                               fragments.__doc__)))


#TODO: update docs to property doc
class Angles(Bonds):
    """Angles for atoms"""
    attrname = 'angles'
    singular = 'angles'


#TODO: update docs to property doc
class Dihedrals(Bonds):
    """Dihedrals for atoms"""
    attrname = 'dihedrals'
    singular = 'dihedrals'


#TODO: update docs to property doc
class Impropers(Bonds):
    attrname = 'impropers'
    singular = 'impropers'
