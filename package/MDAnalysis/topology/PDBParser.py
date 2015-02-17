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
PDB topology parser
===================

Use a PDB file to build a minimum internal structure representation.

.. Warning:: Only cares for atoms and their names; neither
             connectivity nor (partial) charges are deduced. Masses
             are guessed and set to 0 if unknown.
"""

try:
    # BioPython is overkill but potentially extensible (altLoc etc)
    import Bio.PDB
except ImportError:
    raise ImportError("Bio.PDB from biopython not found. Required for PDB->PSF parser.")

from .base import TopologyReader
from MDAnalysis.core.AtomGroup import Atom
import MDAnalysis.coordinates.pdb.extensions
from MDAnalysis.topology.core import guess_atom_type, guess_atom_mass, guess_atom_charge


class PDBParser(TopologyReader):
    def parse(self):
        """Parse atom information from PDB file *pdbfile*.

        Only reads the list of atoms.

        This functions uses the :class:`Bio.PDB.PDBParser` as used by
        :func:`MDAnalysis.coordinates.pdb.extensions.get_structure`.
        
        :Returns: MDAnalysis internal *structure* dict

        .. SeeAlso:: The *structure* dict is defined in
           `MDAnalysis.topology`.
        """
        structure = {}
        # use Sloppy PDB parser to cope with big PDBs!
        pdb = MDAnalysis.coordinates.pdb.extensions.get_structure(self.filename, "0UNK")

        structure['_atoms'] = self._parseatoms(pdb)

        return structure

    def _parseatoms(self, pdb):
        atoms = []

        # translate Bio.PDB atom objects to MDAnalysis Atom.
        for iatom, atom in enumerate(pdb.get_atoms()):
            residue = atom.parent
            chain_id = residue.parent.id

            atomname = atom.name
            atomtype = guess_atom_type(atomname)
            resname = residue.resname
            resid = int(residue.id[1])
            # no empty segids (or Universe throws IndexError)
            segid = residue.get_segid().strip() or chain_id or "SYSTEM"
            mass = guess_atom_mass(atomname)
            charge = guess_atom_charge(atomname)
            bfactor = atom.bfactor
            # occupancy = atom.occupancy

            atoms.append(Atom(iatom, atomname, atomtype, resname, resid, segid,
                              mass, charge, bfactor=bfactor))

        return atoms
