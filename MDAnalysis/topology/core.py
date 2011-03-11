"""
:mod:`MDAnalysis.topology.core` --- Common functions for topology building
==========================================================================

.. function:: build_segments
.. function:: build_bondlist
.. class:: Bond

.. function:: get_parser_for
.. function:: guess_format

.. function:: guess_atom_type
.. function:: guess_atom_mass
.. function:: guess_atom_charge

"""
import os.path
import MDAnalysis.topology

def build_segments(atoms):
    from MDAnalysis.core.AtomGroup import Residue, Segment
    struc = {}
    residues = []
    resatomlist = []
    curr_segname = atoms[0].segid
    curr_resnum = atoms[0].resid
    curr_resname = atoms[0].resname
    for a in atoms:
        if (a.segid == curr_segname):
            if (a.resid == curr_resnum):
                resatomlist.append(a)
            else:
                # New residue
                residues.append(Residue(curr_resname, curr_resnum, resatomlist))
                resatomlist = [a]
                curr_resnum = a.resid
                curr_resname = a.resname
        else:
            # We've come to a new segment
            residues.append(Residue(curr_resname, curr_resnum, resatomlist))
            struc[curr_segname] = Segment(curr_segname, residues)
            residues = []
            resatomlist = [a]
            curr_resnum = a.resid
            curr_resname = a.resname
            curr_segname = a.segid
    # Add the last segment
    residues.append(Residue(curr_resname, curr_resnum, resatomlist))
    struc[curr_segname] = Segment(curr_segname, residues)
    return struc

class Bond(object):
    def __init__(self, a1, a2):
        self.atom1 = a1
        self.atom2 = a2
    def partner(self, atom):
        if atom is self.atom1:
            return self.atom2
        else: return self.atom1
    def length(self):
        bond = self.atom1.pos - self.atom2.pos
        import math
        return math.sqrt((bond[0]**2)+(bond[1]**2)+(bond[2]**2))

def build_bondlists(atoms, bonds):
    for a in atoms:
        a.bonds = []
    for a1, a2 in bonds:
        atom1 = atoms[a1]
        atom2 = atoms[a2]
        b = Bond(atom1, atom2)
        atom1.bonds.append(b)
        atom2.bonds.append(b)

def get_parser_for(filename, permissive=False):
    """Return the appropriate topology parser for *filename*."""
    format = guess_format(filename)
    if permissive:
        return MDAnalysis.topology._topology_parsers_permissive[format]
    return MDAnalysis.topology._topology_parsers[format]

def guess_format(filename):
    """Returns the type of topology file *filename*.

    The current heuristic simply looks at the filename extension but
    more complicated probes could be implemented here or in the
    individual packages (e.g. as static methods).
    """

    # simple extension checking... something more complicated is left
    # for the ambitious
    root, ext = os.path.splitext(filename)
    try:
        if ext.startswith('.'):
            ext = ext[1:]
        format = ext.upper()
    except:
        raise TypeError("Cannot determine topology type for %r" % filename)
    
    if not format in MDAnalysis.topology._topology_parsers:
        raise TypeError("Unknown topology extension %r from %r; only %r are implemented in MDAnalysis." % 
                        (format, filename, MDAnalysis.topology._topology_parsers.keys()))
    return format

# following guess_* used by PDB parser

def guess_atom_type(atomname):
    """Guess atom type from the name.

    Looks in dict to see if element is found, otherwise it uses the first character in the atomname.
    The table comes from CHARMM and AMBER atom types, where the first character is not sufficient to
    determine the atom type.
    This will probably result in some mistakes, but it still better than nothing!
    """
    types = {   'CAL': 'CA',
                'C0': 'CA',
                'CES': 'CS',
                'CLA': 'CL',
                'CLAL': 'CL',
                'FE': 'FE',
                'LIT': 'LI',
                'MG': 'MG',
                'HE': 'HE',
                'NE': 'NE',
                'POT': 'K',
                'SOD': 'NA',
                'ZN': 'ZN',
                'MG': 'MG',
                'CU': 'CU',
                'QK': 'K',
                'QC': 'CE',
                'QL': 'LI',
                'QN': 'NA',
                'QR': 'RB',
                'BC': 'C',
                'AC': 'C',
            }
    if atomname in types:
        return types[atomname]
    else:
        return atomname[0]

def guess_atom_mass(atomname):
    """Guess a mass based on the atom name.

    Masses are hard-coded here and some simple heuristics are used to
    distinguish e.g. CL from C or NA from N. Generally, the first
    letter of the atom name is used to decide.

    .. warning:: Anything not recognized is simply set to 0; if you rely on the
                 masses you might want to double check.
    """
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
    elif atomname[0] == 'H' or (atomname[0] in ('1','2','3','4') and atomname[1] == 'H'):
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
    """Guess atom charge from the name.

    Not implemented; simply returns 0.
    """
    # TODO: do something slightly smarter, at least use name/element
    return 0.0
