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
PDB topology parser
===================

Use a PDB file to build a minimum internal structure representation.

.. Warning:: Only cares for atoms and their names; neither
             connectivity nor (partial) charges are deduced. Masses
             are guessed and set to 0 if unknown.
"""

import os.path
try:
    # BioPython is overkill but potentially extensible (altLoc etc)
    import Bio.PDB
except ImportError:
    raise ImportError("Bio.PDB from biopython not found. Required for PDB->PSF parser.")

import MDAnalysis.coordinates.pdb.extensions
from MDAnalysis.topology.core import guess_atom_type, guess_atom_mass, guess_atom_charge

class PDBParseError(Exception):
    """Signifies an error while reading a PDB file."""
    pass

def parse(pdbfile):
    """Parse atom information from PDB file *pdbfile*.

    Only reads the list of atoms.

    This functions uses the :class:`Bio.PDB.PDBParser` as used by
    :func:`MDAnalysis.coordinates.pdb.extensions.get_structure`.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse`.
    """
    root,ext = os.path.splitext(pdbfile)
    if ext.lower() not in ('.pdb', '.ent'):
        raise PDBParseError("%(pdbfile)r is probably not in PDB format (wrong extension).")
    structure = {}
    # use Sloppy PDB parser to cope with big PDBs!
    pdb =  MDAnalysis.coordinates.pdb.extensions.get_structure(pdbfile,"0UNK")

    __parseatoms_(pdb, structure)
    # TODO: reconstruct bonds from CONECT or guess from distance search
    #       (e.g. like VMD)
    return structure

def __parseatoms_(pdb, structure):
    from MDAnalysis.core.AtomGroup import Atom
    attr = "_atoms"  # name of the atoms section
    atoms = []       # list of Atom objects

    # translate Bio.PDB atom objects to MDAnalysis Atom.
    for iatom,atom in enumerate(pdb.get_atoms()):
        residue = atom.parent
        chain_id = residue.parent.id

        atomname = atom.name
        atomtype = guess_atom_type(atomname)
        resname = residue.resname
        resid = residue.id[1]
        segid = residue.get_segid().strip() or chain_id or "SYSTEM"  # no empty segids (or Universe throws IndexError)
        mass = guess_atom_mass(atomname)
        charge = guess_atom_charge(atomname)
        bfactor = atom.bfactor
        occupancy = atom.occupancy

        atoms.append(Atom(iatom,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge),
                          bfactor=bfactor))

    structure[attr] = atoms

