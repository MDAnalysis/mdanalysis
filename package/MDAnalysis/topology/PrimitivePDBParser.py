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
   if unknown. Partial charges are not set.

.. SeeAlso::

   * :mod:`MDAnalysis.topology.ExtendedPDBParser`
   * :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`
   * :class:`MDAnalysis.core.AtomGroup.Universe`

Classes
-------

.. autoclass:: PrimitivePDBParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

from ..core.AtomGroup import Atom
from .core import (guess_atom_type, guess_atom_mass,
                   guess_atom_charge)
from ..core.util import openany
from .base import TopologyReader


class PrimitivePDBParser(TopologyReader):
    """Parser that obtains a list of atoms from a standard PDB file.

    .. versionadded:: 0.8
    """
    def __init__(self, filename, **kwargs):
        super(PrimitivePDBParser, self).__init__(filename, **kwargs)
        from ..coordinates.PDB import PrimitivePDBReader
        self.PDBReader = PrimitivePDBReader

    def parse(self):
        """Parse atom information from PDB file *filename*.

        :Returns: MDAnalysis internal *structure* dict

        .. SeeAlso:: The *structure* dict is defined in
                     `MDAnalysis.topology` and the file is
                     read with
                     :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.
        """
        structure = {}
        try:
            pdb = self.PDBReader(self.filename)
        except ValueError:
            raise IOError("Failed to open and read PDB file")

        atoms = self._parseatoms(pdb)
        structure['atoms'] = atoms

        bonds = self._parsebonds(pdb, atoms)
        structure['bonds'] = bonds

        return structure

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
                                  altLoc=altLoc, universe=self._u))
            # TER atoms
            #elif len(atom.__dict__) == 5:
            #    pass
            #    #atoms.append(None)
        return atoms

    def _parsebonds(self, primitive_pdb_reader, atoms):
        # Mapping between the atom array indicies a.number and atom ids
        # (serial) in the original PDB file

        mapping = dict((a.serial, a.number) for a in atoms)

        bonds = set()
        with openany(self.filename, "r") as fname:
            lines = ((num, line[6:].split()) for num, line in enumerate(fname)
                     if line[:6] == "CONECT")
            for num, bond in lines:
                atom, atoms = int(bond[0]), map(int, bond[1:])
                for a in atoms:
                    bond = tuple([mapping[atom], mapping[a]])
                    bonds.add(bond)

        bonds = tuple(bonds)

        return bonds
