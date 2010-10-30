# Primitive CRD parser
from __future__ import with_statement

from MDAnalysis.core.AtomGroup import Atom
from MDAnalysis.core.util import FORTRANReader
import core

extformat = FORTRANReader('2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10')
stdformat = FORTRANReader('2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5')

def parse(filename):
    atoms = []
    atom_serial = 0
    with open(filename) as crd:
        for linenum,line in enumerate(crd):
            # reading header
            if line.split()[0] == '*':
                continue 
            elif line.split()[-1] == 'EXT' and bool(int(line.split()[0])) == True: 
                extended = True
                continue 
            elif line.split()[0] == line.split()[-1] and line.split()[0] != '*':
                extended = False 
                continue 
            # anything else should be an atom
            try:
                if extended:
                    # 2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10
                    serial,TotRes,resName,name,x,y,z,chainID,resSeq,tempFactor = extformat.read(line)
                else:
                    # 2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5
                    serial,TotRes,resName,name,x,y,z,chainID,resSeq,tempFactor = stdformat.read(line)
            except:
                print "Check CRD format at line %d: %s" % (linenum, line.rstrip())
                raise #IOError("Check CRD format at line %d: %s" % (linenum, line.rstrip()))

            atomtype = core.guess_atom_type(name)
            mass =  core.guess_atom_mass(name)
            charge = core.guess_atom_charge(name)
            atom_desc = Atom(atom_serial,name,atomtype,resName,TotRes,chainID,mass,charge)
            atoms.append(atom_desc)
            atom_serial += 1

	structure = {}
	structure["_atoms"] = atoms
	# Other attributes are not read since they are not included in .crd files
	other_attrs = ["_bonds" , "_angles" , "_dihe" , "_impr" , "_donors" , "_acceptors"]
	for attr in other_attrs:
		structure[attr] = []
	return structure

