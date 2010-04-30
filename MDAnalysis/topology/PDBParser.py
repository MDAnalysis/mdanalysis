"""Use a PDB file to build a minimum internal structure representation.

.. Warning:: Only cares for atoms and their names; no connectivity is
             deduced, and currently neither (partial) charges nor
             masses are correct.
"""

try:
    # BioPython is overkill but potentially extensible (altLoc etc)
    import Bio.PDB
except ImportError:
    raise ImportError("Bio.PDB from biopython not found. Required for PDB->PSF parser.")

import MDAnalysis.coordinates.pdb.extensions

class PDBParseError(Exception):
    pass

def parse(pdbfile):
    """Parse atom information from PDB file *pdbfile*.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse`.
    """
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
        chain = residue.parent

        atomname = atom.name
        atomtype = guess_atom_type(atomname)
        resname = residue.resname
        resid = residue.id[1]
        segid = residue.get_segid() or "SYSTEM"  # no empty segids (or Universe throws IndexError)
        mass = guess_atom_mass(atomname)
        charge = guess_atom_charge(atomname)

        atoms.append(Atom(iatom,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge)))

    structure[attr] = atoms

def guess_atom_type(atomname):
    # TODO: do something slightly smarter, at least use name/element
    return 0

def guess_atom_mass(atomname):
    # TODO: do something slightly smarter, at least use name/element & dict
    return 12.000

def guess_atom_charge(atomname):
    # TODO: do something slightly smarter, at least use name/element
    return 0.0
