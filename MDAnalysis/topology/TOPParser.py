# $Id: TOPParser.py 325 2011-01-22 14:58:19Z Elizabeth Denning $
"""
TOPParser - reads a Amber top file to build the system
======================================================

The format is defined in `PARM parameter/topology file specification`_.

.. Note:: The Amber charge is converted to electron charges as used in
          MDAnalysis and other packages. To get back Amber charges,
          multiply by 18.2223.

.. _`PARM parameter/topology file specification`:
   http://ambermd.org/formats.html#topology
"""

import MDAnalysis.core
import MDAnalysis.core.units
from MDAnalysis.core import util

class TOPParseError(Exception):
    pass

def parse(topfilename):
    # Open and check top validity
    ######  Reading header info POINTERS  #################
    topfile = open(topfilename,'r')
    next_line = skip_line = topfile.next
    header = next_line()
    if header[:3] != "%VE":
        raise TOPParseError("%s is not a valid TOP file" % topfile)
    title = next_line().split()
    if not (title[1] == "TITLE"):
        raise TOPParseError("%s is not a valid TOP file" % topfile) 
    while header[:14] != '%FLAG POINTERS':
	header = next_line()
    header = next_line()
    topremarks = [next_line().strip() for i in range(4)]
    sys_info = []
    for i in topremarks:
 	j = i.split()
	for k in j:
	    sys_info.append(int(k))
    ########################################################

    structure = {}
    final_structure = {}

    def parse_sec(section_info):
        desc, atoms_per, per_line, parsefunc, data_struc, sect_num = section_info
        from math import ceil as c
        # Get the number
        num = sys_info[sect_num]
	if data_struc in ["_resname","_bond"]:
		x = 1
	else:
		header = next_line()
        
        # Now figure out how many lines to read
        numlines = int((float(num)/per_line))
        #print data_struc, numlines
	if parsefunc == __parsebond_:
		parsefunc(next_line, atoms_per, data_struc, final_structure, numlines)
	else:
		parsefunc(next_line, atoms_per, data_struc, structure, numlines)

    sections = [("ATOM_NAME", 1, 20, __parseatoms_, "_name",0),
                ("CHARGE",1, 5, __parsesection_,"_charge",0),
		("MASS",1, 5, __parsesection_,"_mass",0),
		("ATOM_TYPE_INDEX", 1, 10, __parsesectionint_,"_atom_type",0),
		("NUMBER_EXCLUDED_ATOMS", 1, 10, __parseskip_,"_skip",8),
		("NONBONDED_PARM_INDEX", 1, 10, __parseskip_,"_skip",8),
		("RESIDUE_LABEL", 1, 20, __parseatoms_, "_resname",11),
		("RESIDUE_POINTER", 2, 10, __parsesectionint_,"_respoint",11),
		("BOND_FORCE_CONSTANT", 1, 5, __parseskip_,"_skip",8),
		("BOND_EQUIL_VALUE", 1, 5, __parseskip_,"_skip",8),
		("ANGLE_FORCE_CONSTANT", 1, 5, __parseskip_,"_skip",8),
		("ANGLE_EQUIL_VALUE", 1, 5, __parseskip_,"_skip",8),
		("DIHEDRAL_FORCE_CONSTANT", 1, 5, __parseskip_,"_skip",8),
		("DIHEDRAL_PERIODICITY", 1, 5, __parseskip_,"_skip",8),
		("DIHEDRAL_PHASE", 1, 5, __parseskip_,"_skip",8),
		("SOLTY", 1, 5, __parseskip_,"_skip",8),
		("LENNARD_JONES_ACOEF", 1, 5, __parseskip_,"_skip",8),
		("LENNARD_JONES_BCOEF", 1, 5, __parseskip_,"_skip",8)]
		#("BONDS_INC_HYDROGEN", 2, 4, __parsebond_, "_bonds",2),
                #("ANGLES_INC_HYDROGEN", 3, 3, __parsesection_, "_angles"),
                #("DIHEDRALS_INC_HYDROGEN", 4, 2, __parsesection_, "_dihe")]
                #("NIMPHI", 4, 2, __parsesection_, "_impr"),
                #("NDON", 2, 4, __parsesection_,"_donors"),
                #("NACC", 2, 4, __parsesection_,"_acceptors")]

    try:
        for info in sections:
            parse_sec(info)
    except StopIteration:
        raise Exception("The TOP file didn't contain the minimum required section of ATOM_NAME")
    # Completing info respoint to include all atoms in last resid
    structure["_respoint"].append(sys_info[0]+1)
    topfile.close()
    
    atoms = [None,]*sys_info[0]
    from MDAnalysis.core.AtomGroup import Atom
    index = 0
    j = 0 
    for i in range(sys_info[0]):
    	if structure.keys()[0] == "_charge":
		index += 1
		charge = MDAnalysis.core.units.convert(structure.items()[0][1][i], 
                                                       'Amber', MDAnalysis.core.flags['charge_unit'])
		if structure.items()[5][1][j] <= index < structure.items()[5][1][j+1]:
			resid = j + 1
			resname = structure.items()[1][1][j]
		else:
			j += 1
			resid = j + 1
			resname = structure.items()[1][1][j]
		mass = structure.items()[2][1][i]
		atomtype = structure.items()[3][1][i]
		atomname = structure.items()[4][1][i]
		segid = 'SYSTEM'  # does not exist in Amber
	else:
		raise TOPParseError("Something is screwy with this top file")
    	atom_desc = Atom(index-1,atomname,atomtype,resname,resid,segid,mass,charge)
    	atoms[i] = atom_desc
    final_structure["_atoms"] = atoms
    final_structure["_numatoms"] = sys_info[0] 
    return final_structure


import operator

def __parseskip_(lines, atoms_per, attr, structure, numlines):
    section = []
    lines()
    #print lax
    while (lines()[:5] != "%FLAG"):
    	x = 1 #l = lines()

def __parsebond_(lines, atoms_per, attr, structure, numlines):
    section = [] #[None,]*numlines
    for i in xrange(numlines):
        l = lines()
        # Subtract 1 from each number to ensure zero-indexing for the atoms
        f = map(int, l.split())
        fields = [a-1 for a in f]
        for j in range(0, len(fields), atoms_per):
            section.append(tuple(fields[j:j+atoms_per]))
    structure[attr] = section

def __parsesectionint_(lines, atoms_per, attr, structure, numlines):
    section = [] #[None,]*numlines
    y = lines().strip("%FORMAT(")
    y.strip(")")
    x = util.FORTRANReader(y)
    liz = 1
    for i in xrange(numlines+1):
	l = lines()
	# Subtract 1 from each number to ensure zero-indexing for the atoms
	try:
	    for j in range(0, len(x.entries)):
		section.append(int(l[x.entries[j].start:x.entries[j].stop].strip()))
    		liz += 1
	except:
		continue
    structure[attr] = section

def __parsesection_(lines, atoms_per, attr, structure, numlines):
    section = [] #[None,]*numlines
    y = lines().strip("%FORMAT(")
    y.strip(")")
    x = util.FORTRANReader(y)
    for i in xrange(numlines+1):
	l = lines()
	# Subtract 1 from each number to ensure zero-indexing for the atoms
	try:
		for j in range(0, len(x.entries)):
	            section.append(float(l[x.entries[j].start:x.entries[j].stop].strip()))
    	except:
    		continue
    structure[attr] = section

def __parseatoms_(lines, atoms_per, attr, structure, numlines):
    section = [] #[None,]*numlines
    y = lines().strip("%FORMAT(")
    y.strip(")")
    x = util.FORTRANReader(y)
    for i in xrange(numlines+1):
	l = lines()
	# Subtract 1 from each number to ensure zero-indexing for the atoms
	for j in range(0, len(x.entries)):
	    section.append(l[x.entries[j].start:x.entries[j].stop].strip())
    structure[attr] = section
