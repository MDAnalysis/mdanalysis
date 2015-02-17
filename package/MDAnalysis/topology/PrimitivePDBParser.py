# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
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

import numpy as np

from MDAnalysis.core.AtomGroup import Atom
from MDAnalysis.topology.core import (guess_atom_type, guess_atom_mass,
                                      guess_atom_charge, guess_bonds)
import MDAnalysis.coordinates.PDB
from MDAnalysis.core.util import openany
from .base import TopologyReader


class PrimitivePDBParser(TopologyReader):
    """Parser that obtains a list of atoms from a standard PDB file.

    .. versionadded:: 0.8
    """
    def __init__(self, filename, **kwargs):
        super(PrimitivePDBParser, self).__init__(filename, **kwargs)
        self.PDBReader = MDAnalysis.coordinates.PDB.PrimitivePDBReader

    def parse(self):
        """Parse atom information from PDB file *filename*.

        :Returns: MDAnalysis internal *structure* dict

        .. SeeAlso:: The *structure* dict is defined in
                     `MDAnalysis.topology` and the file is
                     read with
                     :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.
        """
        self.structure = {}
        try:
            pdb = self.PDBReader(self.filename)
        except ValueError:
            raise IOError("Failed to open and read PDB file")

        self._parseatoms(pdb)
        self._parsebonds(pdb)
        return self.structure

    def _parseatoms(self, pdb):
        atoms = []

        # translate list of atoms to MDAnalysis Atom.
        for iatom, atom in enumerate(pdb._atoms):
            # ATOM
            if len(atom.__dict__) == 10:
                atomname = atom.name
                atomtype = atom.element or guess_atom_type(atomname)
                resname = atom.resName
                resid = int(atom.resSeq)
                chain = atom.chainID.strip()
                # no empty segids (or Universe throws IndexError)
                segid = atom.segID.strip() or chain or "SYSTEM"
                mass = guess_atom_mass(atomname)
                charge = guess_atom_charge(atomname)
                bfactor = atom.tempFactor
                # occupancy = atom.occupancy
                altLoc = atom.altLoc

                atoms.append(Atom(iatom, atomname, atomtype, resname, resid,
                                  segid, mass, charge,
                                  bfactor=bfactor, serial=atom.serial,
                                  altLoc=altLoc))
            # TER atoms
            #elif len(atom.__dict__) == 5:
            #    pass
            #    #atoms.append(None)
        self.structure["_atoms"] = atoms

    def _parsebonds(self, primitive_pdb_reader):
        if self.guess_bonds_mode:
            guessed_bonds = guess_bonds(self.structure["_atoms"],
                                        np.array(primitive_pdb_reader.ts))
            self.structure["_guessed_bonds"] = guessed_bonds

        # Mapping between the atom array indicies a.number and atom ids
        # (serial) in the original PDB file

        mapping = dict((a.serial, a.number) for a in self.structure["_atoms"])

        bonds = set()
        with openany(self.filename, "r") as fname:
            lines = ((num, line[6:].split()) for num, line in enumerate(fname)
                     if line[:6] == "CONECT")
            for num, bond in lines:
                atom, atoms = int(bond[0]), map(int, bond[1:])
                for a in atoms:
                    bond = tuple([mapping[atom], mapping[a]])
                    bonds.add(bond)

        self.structure["_bonds"] = tuple(bonds)
