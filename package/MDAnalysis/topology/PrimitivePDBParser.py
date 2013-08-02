# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Primitive PDB topology parser
=============================

Use a PDB file to build a minimum internal structure representation (list of atoms).

Reads a PDB file line by line and is not fuzzy about numbering.

.. Warning:: Only cares for atoms and their names; neither
             connectivity nor (partial) charges are deduced. Masses
             are guessed and set to 0 if unknown.
"""

from MDAnalysis.topology.core import guess_atom_type, guess_atom_mass, guess_atom_charge, guess_bonds
import numpy as np
import MDAnalysis.coordinates.PDB


class PDBParseError(Exception):
    """Signifies an error during parsing a PDB file."""
    pass

def parse(filename):
    """Parse atom information from PDB file *filename*.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                 :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.
    """
    structure = {}
    pdb =  MDAnalysis.coordinates.PDB.PrimitivePDBReader(filename)

    __parseatoms_(pdb, structure)
    # TODO: reconstruct bonds from CONECT or guess from distance search
    #       (e.g. like VMD)
    __parsebonds_(filename, pdb, structure, guess_bonds_mode=False)
    return structure


def parse_bonds(filename):
    """Parse atom information from PDB file *filename*.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                 :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.
    """
    structure = {}
    pdb =  MDAnalysis.coordinates.PDB.PrimitivePDBReader(filename)

    __parseatoms_(pdb, structure)
    # TODO: reconstruct bonds from CONECT or guess from distance search
    #       (e.g. like VMD)
    __parsebonds_(filename, pdb, structure, guess_bonds_mode=True)
    return structure

def __parseatoms_(pdb, structure):
    from MDAnalysis.core.AtomGroup import Atom
    attr = "_atoms"  # name of the atoms section
    atoms = []       # list of Atom objects

    # translate list of atoms to MDAnalysis Atom.
    for iatom,atom in enumerate(pdb._atoms):

        # ATOM
        if len(atom.__dict__) == 9:
            atomname = atom.name
            atomtype = atom.element or guess_atom_type(atomname)
            resname = atom.resName
            resid = atom.resSeq
            chain = atom.chainID.strip()
            segid = atom.segID.strip() or chain or "SYSTEM"  # no empty segids (or Universe throws IndexError)
            mass = guess_atom_mass(atomname)
            charge = guess_atom_charge(atomname)
            bfactor = atom.tempFactor
            occupancy = atom.occupancy

            atoms.append(Atom(iatom,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge),\
                              bfactor=bfactor,serial=atom.serial))
        # TER atoms
        elif len(atom.__dict__) == 5:
            pass
            #atoms.append(None)
    structure[attr] = atoms



def __parsebonds_(filename, primitive_pdb_reader, structure, guess_bonds_mode):
    guessed_bonds = set()
    if guess_bonds_mode:
        guessed_bonds = guess_bonds(structure["_atoms"], np.array(primitive_pdb_reader.ts))

    #
    # Mapping between the atom array indicies and atom ids in the original PDB file
    #
    mapping =  dict([(a.serial, i+1) for i, a in  enumerate(structure["_atoms"])])

    bonds = set()
    with open(filename , "r") as filename:
        lines = [(num, line[6:].split()) for num,line in enumerate(filename) if line[:6] == "CONECT"]
        for num, bond in lines:
            atom, atoms = int(bond[0]) , map(int,bond[1:])
            for a in atoms:
                bond = frozenset([mapping[atom], mapping[a] ])
                bonds.add(bond)

    # FIXME by JD: we could use a BondsGroup class perhaps
    structure["_bonds"] = bonds
    structure["_guessed_bonds"] = guessed_bonds
