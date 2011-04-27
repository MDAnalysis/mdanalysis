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
#     Molecular Dynamics Simulations. J. Comput. Chem. (2011),
#     doi:10.1002/jcc.21787
#

"""
PDBQT topology parser
=============================

Use a PDBQT file to build a minimum internal structure representation (list of atoms).

Reads a PDBQT file line by line and is not fuzzy about numbering.

This does not support multi-model PDBQT files (yet!)

.. Warning:: Only cares for atoms and their names; neither
             connectivity nor (partial) charges are deduced. Masses
             are guessed and set to 0 if unknown.
"""

import MDAnalysis.coordinates.PDBQT
from MDAnalysis.topology.core import guess_atom_type, guess_atom_mass, guess_atom_charge

class PDBQTParseError(Exception):
    """Signifies an error during parsing a PDBQT file."""
    pass

def parse(filename):
    """Parse atom information from PDBQT file *filename*.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                 :class:`MDAnalysis.coordinates.PDBQT.PDBQTReader`.
    """
    structure = {}
    pdb =  MDAnalysis.coordinates.PDBQT.PDBQTReader(filename)

    __parseatoms_(pdb, structure)
    # TODO: reconstruct bonds from CONECT or guess from distance search
    #       (e.g. like VMD)
    return structure

def __parseatoms_(pdb, structure):
    from MDAnalysis.core.AtomGroup import Atom # Added partialCharge parameter to the Atom classes' constructor
    attr = "_atoms"  # name of the atoms section
    atoms = []       # list of Atom objects

    # translate list of atoms to MDAnalysis Atom.
    for iatom,atom in enumerate(pdb._atoms):
        atomname = atom.name
        atomtype = atom.element or guess_atom_type(atomname)
        resname = atom.resName
        resid = atom.resSeq
        chain = atom.chainID.strip()
        segid = chain or "SYSTEM"  # no empty segids (or Universe throws IndexError)
        mass = guess_atom_mass(atomname)
        charge = guess_atom_charge(atomname)
        bfactor = atom.tempFactor
        occupancy = atom.occupancy
        partialCharge = atom.partialCharge
        element = atom.element
        atoms.append(Atom(iatom,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge),
                          bfactor=bfactor,partialCharge=partialCharge, element=element))

    structure[attr] = atoms
