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
# doi: 10.25080/majora-629e541a-00e
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
from collections import namedtuple
from ...guesser.tables import Z2SYMB

TpxHeader = namedtuple(
    "TpxHeader", [
        "ver_str", "precision",
        "fver", "fgen", "file_tag", "natoms", "ngtc", "fep_state", "lamb",
        "bIr", "bTop", "bX", "bV", "bF", "bBox", "sizeOfTprBody"])
Box = namedtuple("Box", "size rel v")
Mtop = namedtuple("Mtop", "nmoltype moltypes nmolblock")
Params = namedtuple("Params", "atnr ntypes functype reppow fudgeQQ")
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
        return (
            f"Molecule: {self.name:<20s} "
            f"#atoms: {self.number_of_atoms():<10d} "
            f"#residues: {self.number_of_residues():<10d}"
        )

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
    def __init__(
            self, id, name, type, resid, resname, mass, charge, atomic_number):
        # id is only within the scope of a single molecule, not the whole system
        self.id = id
        self.name = name
        self.type = type
        self.resid = resid
        self.resname = resname
        self.mass = mass
        self.charge = charge
        self.atomic_number = atomic_number

    @property
    def element_symbol(self):
        """
        The symbol of the atom element.

        The symbol corresponding to the atomic number. If the atomic number
        is not recognized, which happens if a particle is not really an
        atom (e.g a coarse-grained particle), an empty string is returned.
        """
        return Z2SYMB.get(self.atomic_number, '')

    def __repr__(self):
        return (
            f"< AtomKind: "
            f"id {self.id:6d}, "
            f"name {self.name:5s}, "
            f"type {self.type:10s}, "
            f"resid {self.resid:6d}, "
            f"resname {self.resname:4s}, "
            f"mass {self.mass:8.4f}, "
            f"charge {6:12.3f} "
            ">"
        )


class InteractionKind(object):
    def __init__(self, name, long_name, natoms):
        """natoms: number of atoms involved in this type of interaction"""
        self.name = name
        self.long_name = long_name
        self.natoms = natoms

    def process(self, atom_ndx):
        # The format for all record is (type, atom1, atom2, ...)
        # but we are only interested in the atoms.
        for cursor in range(0, len(atom_ndx), self.natoms + 1):
            yield atom_ndx[cursor + 1: cursor + 1 + self.natoms]
