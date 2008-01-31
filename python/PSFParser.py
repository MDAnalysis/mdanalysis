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
        header = next_line()
        while header.strip() == "": header = next_line()
        header = header.split()
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

    try:
        for info in sections:
            skip_line()
            parse_sec(info)
    except StopIteration:
        # Reached the end of the file before we expected
        if not structure.has_key("_atoms"):
            raise Exception("The PSF file didn't contain the minimum required section of NATOM")
    # Who cares about the rest
    psffile.close()
    return structure

def __parseatoms_(lines, atoms_per, attr, structure, numlines):
    """Parses atom section in a Charmm PSF file.

        Currently the only supported PSF format is the standard format
        (which does NOT include the EXT keyword in the header. CHEQ is
        supported in the sense that CHEQ data is simply ignored.
        
        ------ notes ------

        Format from source/psffres.src

        CHEQ:
        II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)

        standard format:
          (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)
          (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)  XPLOR
        expanded format EXT:
          (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)
          (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8,2G14.6) XPLOR
        
        no CHEQ:
        II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I)

        standard format:
          (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
          (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)  XPLOR
        expanded format EXT:
          (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)
          (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8) XPLOR
          """
    atoms = [None,]*numlines
    from AtomGroup import Atom

    # Oli: I don't think that this is the correct OUTPUT format:
    #   psf_atom_format = "   %5d %4s %4d %4s %-4s %-4s %10.6f      %7.4f%s\n"
    # It should be rather something like:
    #   psf_ATOM_format = '%(iatom)8d %(segid)4s %(resid)-4d %(resname)4s '+\
    #                     '%(name)-4s %(type)4s %(charge)-14.6f%(mass)-14.4f%(imove)8d\n'

    # source/psfres.src (CHEQ but not EXTended), see comments above
    #   II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)
    #  (I8,1X,A4, 1X,A4,  1X,A4,  1X,A4,  1X,I4,  1X,2G14.6,     I8,   2G14.6)
    #   0:8   9:13   14:18   19:23   24:28   29:33   34:48 48:62 62:70 70:84 84:98

    for i in xrange(numlines):
        l = lines()
        # 0     1      2      3        4         5       6       7      8
        iatom, segid, resid, resname, atomname, atomtype, charge, mass, imove =  l[:8], l[9:13].strip(), l[14:18], l[19:23].strip(), l[24:28].strip(), l[29:33].strip(), l[34:48], l[48:62], l[62:70]   #  l[70:84], l[84:98] ignore ECH and EHA
        # Atom(atomno, atomname, type, resname, resid, segid, mass, charge)
        # We want zero-indexing for atom numbers to make it easy
        atom_desc = Atom(int(iatom)-1,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge))
        atoms[i] = atom_desc

    structure[attr] = atoms

import operator
def __parsesection_(lines, atoms_per, attr, structure, numlines):
    section = [] #[None,]*numlines
    #for l in lines:
    for i in xrange(numlines):
        l = lines()
        # Subtract 1 from each number to ensure zero-indexing for the atoms
        f = map(int, l.split())
        fields = [a-1 for a in f]
        for j in range(0, len(fields), atoms_per):
            section.append(tuple(fields[j:j+atoms_per]))
    structure[attr] = section

def build_segments(atoms):
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

def build_bondlists(atoms, bonds):
    for a in atoms:
        a.bonds = []
    for a1, a2 in bonds:
        atom1 = atoms[a1]
        atom2 = atoms[a2]
        b = Bond(atom1, atom2)
        atom1.bonds.append(b)
        atom2.bonds.append(b)
