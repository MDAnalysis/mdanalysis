# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; encoding: utf-8 -*-
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

This topology parser uses a standard PDB file to build a minimum
internal structure representation (list of atoms).

The topology reader reads a PDB file line by line and ignores atom
numbers but only reads residue numbers up to 9,999 correctly. If you
have systems containing at least 10,000 residues then you need to use
a different file format (e.g. the "extended" PDB, *XPDB* format, see
:mod:`~MDAnalysis.topology.ExtendedPDBParser`) that can handle residue
numbers up to 99,999.

.. Note::

   The parser processes atoms and their names. Masses are guessed and set to 0
   if unknown. Partial charges are not set. Bond connectivity can be guessed if
   the ``bonds=True`` keyword is set for
   :class:`~MDAnalysis.core.AtomGroup.Universe`.

.. SeeAlso::

   * :mod:`MDAnalysis.topology.ExtendedPDBParser`
   * :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`
   * :class:`MDAnalysis.core.AtomGroup.Universe`

"""

from MDAnalysis.topology.core import guess_atom_type, guess_atom_mass, guess_atom_charge, guess_bonds
import numpy as np
import MDAnalysis.coordinates.PDB


class PDBParseError(Exception):
    """Signifies an error during parsing a PDB file."""
    pass


class PrimitivePDBParser(object):
    """Parser that obtains a list of atoms from a standard PDB file.

    .. versionadded:: 0.8
    """

    def __init__(self, filename, guess_bonds_mode=False):
        self.PDBReader = MDAnalysis.coordinates.PDB.PrimitivePDBReader
        self.filename = filename
        self.guess_bonds_mode = guess_bonds_mode

    def parse(self):
        """Parse atom information from PDB file *filename*.

        :Returns: MDAnalysis internal *structure* dict

        .. SeeAlso:: The *structure* dict is defined in
                     :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                     :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.
        """
        self.structure = {}
        pdb =  self.PDBReader(self.filename)

        self._parseatoms(pdb)
        # TODO: reconstruct bonds from CONECT or guess from distance search
        #       (e.g. like VMD)
        self._parsebonds(self.filename, pdb)
        return self.structure

    def _parseatoms(self, pdb):
        from MDAnalysis.core.AtomGroup import Atom
        attr = "_atoms"  # name of the atoms section
        atoms = []       # list of Atom objects

        # translate list of atoms to MDAnalysis Atom.
        for iatom,atom in enumerate(pdb._atoms):

            # ATOM
            if len(atom.__dict__) == 10:
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
                altLoc = atom.altLoc
                
                atoms.append(Atom(iatom,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge),\
                                  bfactor=bfactor,serial=atom.serial, altLoc=altLoc))
            # TER atoms
            elif len(atom.__dict__) == 5:
                pass
                #atoms.append(None)
        self.structure[attr] = atoms

    def _parsebonds(self, filename, primitive_pdb_reader):
        guessed_bonds = set()
        if self.guess_bonds_mode:
            guessed_bonds = guess_bonds(self.structure["_atoms"], np.array(primitive_pdb_reader.ts))

        #
        # Mapping between the atom array indicies a.number and atom ids (serial) in the original PDB file
        #
        mapping =  dict((a.serial, a.number) for a in  self.structure["_atoms"])

        bonds = set()
        with open(filename , "r") as filename:
            lines = ((num, line[6:].split()) for num,line in enumerate(filename) if line[:6] == "CONECT")
            for num, bond in lines:
                atom, atoms = int(bond[0]) , map(int,bond[1:])
                for a in atoms:
                    bond = frozenset([mapping[atom], mapping[a] ])
                    bonds.add(bond)

        # FIXME by JD: we could use a BondsGroup class perhaps
        self.structure["_bonds"] = bonds
        self.structure["_guessed_bonds"] = guessed_bonds

# function to keep compatible with the current API; should be cleaned up...
def parse(filename):
    """Parse atom information from PDB file *filename*.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                 :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.

    """
    return PrimitivePDBParser(filename).parse()

def parse_bonds(filename):
    """Parse atom information from PDB file *filename* and guesses bonds.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                 :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.
    """
    return PrimitivePDBParser(filename, guess_bonds_mode=True).parse()
