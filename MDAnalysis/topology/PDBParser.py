"""Use a PDB file to build a minimum internal structure representation.

.. Warning:: Only cares for atoms and their names; no connectivity is
             deduced, and currently neither (partial) charges nor
             masses are correct.
"""

import os.path
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
    if atomname[0] == 'N':
    	if atomname[:1] != 'NA':
    		return 14.007
    elif atomname[0] == 'C':
    	if atomname[:2] not in ['CAL','CL ','CLA']:
    		return 12.010
    elif atomname[0] == 'O':
    	return 15.999
    elif atomname[0] == 'S':
    	if atomname[:2] != 'SOD':
		return 32.065
    elif atomname[0] == 'P':
    	return 30.974
    elif atomname[0] == 'H':
    	return 1.008 
    elif atomname[:1] == 'MG':
    	return 24.305
    elif atomname[:2] in ['K  ','POT']:
    	return 39.102
    elif atomname[:1] == 'CL':
    	return 35.450
    elif atomname[:2] in ['NA ','SOD']:
    	return 22.989 
    else:
    	return 0.000

def guess_atom_charge(atomname):
    # TODO: do something slightly smarter, at least use name/element
    return 0.0
