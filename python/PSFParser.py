""" PSFParser - reads a psf file to build the system
"""

class PSFParseError(Exception):
    pass

def parse(psffilename):
    # Open and check psf validity
    psffile = open(psffilename,'r')
    next_line = skip_line = psffile.next
    header = next_line()
    if header[:3] != "PSF":
        raise PSFParseError("%s is not a valid PSF file")
    skip_line()
    title = next_line().split()
    if not (title[1] == "!NTITLE"):
        raise PSFParseError("%s is not a valid PSF file") 
    psfremarks = [next_line() for i in range(int(title[0]))]

    structure = {}

    def parse_sec(section_info):
        desc, atoms_per, per_line, parsefunc, data_struc = section_info
        from math import ceil as c
        header = next_line().split()
        # Get the number
        num = int(header[0])
        sect_type = header[1].strip('!:')
        # Make sure the section type matches the desc
        if not (sect_type == desc): raise PSFParseError("Something is screwy with this psf file")
        # Now figure out how many lines to read
        numlines = int(c(float(num)/per_line))
        # Too bad I don't have generator expressions
        #def repeat(func, num):
        #    for i in xrange(num):
        #        yield func()
        #lines = repeat(next_line, numlines)
        parsefunc(next_line, atoms_per, data_struc, structure, numlines)

    sections = [("NATOM", 1, 1, __parseatoms_, "_atoms"),
                ("NBOND", 2, 4, __parsesection_, "_bonds"),
                ("NTHETA", 3, 3, __parsesection_, "_angles"),
                ("NPHI", 4, 2, __parsesection_, "_dihe"),
                ("NIMPHI", 4, 2, __parsesection_, "_impr")]

    for info in sections:
        skip_line()
        parse_sec(info)

    # Who cares about the rest
    psffile.close()
    return structure

def __parseatoms_(lines, atoms_per, attr, structure, numlines):
    atoms = [None,]*numlines
    from AtomGroup import Atom
    for i in xrange(numlines):
    #for l in lines:
        l = lines()
        fields = l.split()
        # Atom(atomno, atomname, type, resname, resid, segid, mass, charge)
        # We want zero-indexing for atom numbers to make it easy
        atom_desc = Atom(int(fields[0])-1, fields[4], fields[5], fields[3], int(fields[2]), fields[1], float(fields[7]), float(fields[6]))
        atoms[i] = atom_desc

    structure[attr] = atoms

import operator
def __parsesection_(lines, atoms_per, attr, structure, numlines):
    section = [None,]*numlines
    #for l in lines:
    for i in xrange(numlines):
        l = lines()
        # Subtract 1 from each number to ensure zero-indexing for the atoms
        f = map(int, l.split())
        fields = [a-1 for a in f]
        for j in range(0, len(fields), atoms_per):
            section[i] = tuple(fields[j:j+atoms_per])
    structure[attr] = section

def _buildstructure_(atoms):
    from AtomGroup import Residue, Segment
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

def _buildbondlists_(atoms, bonds):
    for a in atoms:
        a.bonds = []
    for a1, a2 in bonds:
        atom1 = atoms[a1]
        atom2 = atoms[a2]
        b = Bond(atom1, atom2)
        atom1.bonds.append(b)
        atom2.bonds.append(b)
