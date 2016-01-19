# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
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


# TPR parser and tpr support module
# Copyright (c) 2011 Zhuyi Xue
# Released under the  GNU Public Licence, v2

"""
Class definitions for the TPRParser
===================================

"""

from six.moves import range

from collections import namedtuple

TpxHeader = namedtuple(
    "TpxHeader", [
        "number", "ver_str", "precision",
        "fver", "fgen", "file_tag", "natoms", "ngtc", "fep_state", "lamb",
        "bIr", "bTop", "bX", "bV", "bF", "bBox"])
Box = namedtuple("Box", "size rel v")
Mtop = namedtuple("Mtop", "nmoltype moltypes nmolblock")
TPRTopology = namedtuple("TPRTopology", "atoms, bonds, angles, dihe, impr")
Params = namedtuple("Params", "atnr ntypes functype reppow fudgeQQ iparams")
Atom = namedtuple("Atom", ["m", "q", "mB", "qB", "tp", "typeB", "ptype", "resind", "atomnumber"])
Atoms = namedtuple("Atoms", "atoms nr nres type typeB atomnames resnames")
Ilist = namedtuple("Ilist", "nr ik, iatoms")
Molblock = namedtuple("Molblock", [
    "molb_type", "molb_nmol", "molb_natoms_mol",
    "molb_nposres_xA", "molb_nposres_xB"])


class MoleculeKind(object):
    def __init__(self, name, atomkinds, bonds=None, angles=None,
                 dihe=None, impr=None, donors=None, acceptors=None):
        self.name = name  # name of the molecule
        self.atomkinds = atomkinds
        self.bonds = bonds
        self.angles = angles
        self.dihe = dihe
        self.impr = impr
        self.donors = donors
        self.acceptors = acceptors

    def __repr__(self):
        return "Molecule: {0:<20s} #atoms: {1:<10d} #residues: {2:<10d}".format(
            self.name, self.number_of_atoms(), self.number_of_residues())

    def number_of_atoms(self):
        return len(self.atomkinds)

    def number_of_residues(self):
        return len({a.resid for a in self.atomkinds})

    # remap_ method returns [tuple(), tuple(), ..] or []
    # Note: in MDAnalysis, bonds, angles, etc are represented as tuple and not as list
    def remap_bonds(self, atom_start_ndx):
        if self.bonds:
            return [tuple(i + atom_start_ndx for i in b) for b in self.bonds]
        else:
            return []

    def remap_angles(self, atom_start_ndx):
        if self.angles:
            return [tuple(i + atom_start_ndx for i in a) for a in self.angles]
        else:
            return []

    def remap_dihe(self, atom_start_ndx):
        if self.dihe:
            return [tuple(i + atom_start_ndx for i in a) for a in self.dihe]
        else:
            return []

    def remap_impr(self, atom_start_ndx):
        # improper
        if self.impr:
            return [tuple(i + atom_start_ndx for i in a) for a in self.impr]
        else:
            return []


class AtomKind(object):
    def __init__(self, id, name, type, resid, resname, mass, charge):
        # id is only within the scope of a single molecule, not the whole system
        self.id = id
        self.name = name
        self.type = type
        self.resid = resid
        self.resname = resname
        self.mass = mass
        self.charge = charge

    def __repr__(self):
        return \
            ("< AtomKind: id {0:6d}, name {1:5s}, type {2:10s}, resid {3:6d}, resname {4:4s}, mass {5:8.4f}, "
             "charge {6:12.3f} >".format(self.id, self.name, self.type, self.resid,
                                         self.resname, self.mass, self.charge))


class InteractionKind(object):
    def __init__(self, name, long_name, natoms):
        """natoms: number of atoms involved in this type of interaction"""
        self.name = name
        self.long_name = long_name
        self.natoms = natoms

    def process(self, atom_ndx):
        while atom_ndx:
            # format for all info: (type, [atom1, atom2, ...])
            # yield atom_ndx.pop(0), [atom_ndx.pop(0) for i in range(self.natoms)]

            # but currently only [atom1, atom2, ...] is interested
            atom_ndx.pop(0)
            yield [atom_ndx.pop(0) for i in range(self.natoms)]
