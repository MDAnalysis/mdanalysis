"""
Primitive PDB topology parser
=============================

Use a PDB file to build a minimum internal structure representation (list of atoms).

Reads a PDB file line by line and is not fuzzy about numbering.

.. Warning:: Only cares for atoms and their names; neither
             connectivity nor (partial) charges are deduced. Masses
             are guessed and set to 0 if unknown.
"""

import MDAnalysis.coordinates.PDB
from MDAnalysis.topology.core import guess_atom_type, guess_atom_mass, guess_atom_charge

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
    return structure

def __parseatoms_(pdb, structure):
    from MDAnalysis.core.AtomGroup import Atom
    attr = "_atoms"  # name of the atoms section
    atoms = []       # list of Atom objects

    # translate list of atoms to MDAnalysis Atom.
    for iatom,atom in enumerate(pdb._atoms):
        atomname = atom.name
        atomtype = atom.element or guess_atom_type(atomname)
        resname = atom.resName
        resid = atom.resSeq
        chain = atom.chainID.strip()
        segid = atom.segID.strip() or chain or "SYSTEM"  # no empty segids (or Universe throws IndexError)
        mass = guess_atom_mass(atomname)
        charge = guess_atom_charge(atomname)

        atoms.append(Atom(iatom,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge)))

    structure[attr] = atoms
