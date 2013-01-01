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

""" LAMMPSParser - reads a lammps data file to build system
"""

import numpy
# from Scientific.IO.FortranFormat import FortranFormat, FortranLine #unused

class Atom(object):
    __slots__ = ("index", "name", "type", "chainid", "charge", "mass", "_positions")
    def __init__(self, index, name, type, chain_id, charge=0, mass=1):
        self.index = index
        self.name = repr(type)
        self.type = type
        self.chainid = chain_id
        self.charge = charge
        self.mass = mass
    def __repr__(self):
        return "<Atom "+repr(self.index+1)+ ": name " + repr(self.type) +" of chain "+repr(self.chainid)+">"
    def __cmp__(self, other):
        return cmp(self.index, other.index)
    def __eq__(self, other):
        return self.index == other.index
    def __hash__(self):
        return hash(self.index)
    def __getattr__(self, attr):
        if attr == 'pos':
            return self._positions[self.index]
        else: super(Atom, self).__getattribute__(attr)
    def __iter__(self):
        pos = self.pos
        return iter((self.index+1, self.chainid, self.type, self.charge, self.mass, pos[0], pos[1], pos[2]))

class LAMMPSParseError(Exception):
    pass

header_keywords= ["atoms","bonds","angles","dihedrals","impropers","atom types","bond types","angle types","dihedral types","improper types","xlo xhi","ylo yhi","zlo zhi"]
connections = dict([["Bonds",("bonds", 3)],["Angles",("angles", 3)],
                    ["Dihedrals",("dihedrals", 4)],["Impropers",("impropers", 2)]])
coeff = dict([["Masses",("atom types", 1)], ["Velocities",("atoms", 3)],
             ["Pair Coeffs",("atom types", 4)],
             ["Bond Coeffs",("bond types", 2)],["Angle Coeffs",("angle types", 4)],
             ["Dihedral Coeffs",("dihedral types", 3)], ["Improper Coeffs",("improper types", 2)]])


def conv_float(l):
    """
    Function to be passed into map or a list comprehension. If the argument is a float it is converted,
    otherwise the original string is passed back
    """
    try:
        n = float(l)
    except ValueError:
        n = l
    return n

class LAMMPSData(object):
    def __init__(self, filename=None):
        self.names = {}
        self.headers = {}
        self.sections = {}
        if filename == None:
            self.title = "LAMMPS data file"
        else:
            # Open and check validity
            file = open(filename, 'r')
            file_iter = file.xreadlines()
            self.title = file_iter.next()
            # Parse headers
            headers = self.headers
            for l in file_iter:
                line = l.strip()
                if len(line) == 0: continue
                found = False
                for keyword in header_keywords:
                    if line.find(keyword) >= 0:
                        found = True
                        values = line.split()
                        if keyword in ("xlo xhi", "ylo yhi", "zlo zhi"):
                            headers[keyword] = (float(values[0]), float(values[1]))
                        else:
                            headers[keyword] = int(values[0])
                if found == False: break
            file.close()
            # Parse sections
            # XXX This is a crappy way to do it
            file = open(filename, 'r')
            file_iter = file.xreadlines()
            # Create coordinate array
            positions = numpy.zeros((headers['atoms'], 3), numpy.Float64)
            sections = self.sections
            for l in file_iter:
                line = l.strip()
                if len(line) == 0: continue
                if coeff.has_key(line):
                    h, numcoeff = coeff[line]
                    # skip line
                    file_iter.next()
                    data = []
                    for i in xrange(headers[h]):
                        fields = file_iter.next().strip().split()
                        data.append(tuple(map(conv_float, fields[1:])))
                    sections[line] = data
                elif connections.has_key(line):
                    h, numfields = connections[line]
                    # skip line
                    file_iter.next()
                    data = []
                    for i in range(headers[h]):
                        fields = file_iter.next().strip().split()
                        data.append(tuple(map(int, fields[1:])))
                    sections[line] = data
                elif line == "Atoms":
                    file_iter.next()
                    data = []
                    for i in xrange(headers["atoms"]):
                        fields = file_iter.next().strip().split()
                        index = int(fields[0])-1
                        a = Atom(index=index, name=fields[2], type=int(fields[2]), chain_id=int(fields[1]), charge=float(fields[3]))
                        a._positions = positions
                        data.append(a)
                        positions[index] = numpy.array([float(fields[4]), float(fields[5]), float(fields[6])])
                    sections[line] = data
                elif line == "Masses":
                    file_iter.next()
                    data = []
                    for i in xrange(headers["atom type"]):
                        fields = file_iter.next().strip().split()
                        print "help"
            self.positions = positions
            file.close()
    def writePSF(self, filename, names=None):
        import string
        # Naveen formatted -- works with MDAnalysis verison 52
        #psf_atom_format = "   %5d %-4s %-4d %-4s %-4s %-4s %10.6f      %7.4f            %1d\n"
        # Liz formatted -- works with MDAnalysis verison 59
        #psf_atom_format = "%8d %4.4s %-4.4s %-4.4s %-4.4s %-4.4s %16.8e %1s %-7.4f %7.7s %s\n"
        # Oli formatted -- works with MDAnalysis verison 81
        psf_atom_format = "%8d %4s %-4s %4s %-4s% 4s %-14.6f%-14.6f%8s\n"
        file = open(filename, 'w')
        file.write("PSF\n\n")
        file.write(string.rjust('0', 8) + ' !NTITLE\n\n')
        file.write(string.rjust(str(len(self.sections["Atoms"])), 8) + ' !NATOM\n')
        #print self.sections["Masses"]
        for i, atom in enumerate(self.sections["Atoms"]):
            if names != None: resname, atomname = names[i]
            else: resname, atomname = 'TEMP', 'XXXX'
            for j, liz in enumerate(self.sections["Masses"]):
                    liz = liz[0]
                    #print j+1, atom.type, liz
                    if j+1 == atom.type: line = [i+1, 'TEMP', str(atom.chainid), resname, atomname, str(atom.type+1), atom.charge, float(liz), 0.]
                    else: continue
            #print line
            file.write(psf_atom_format%tuple(line))

        file.write("\n")
        num_bonds = len(self.sections["Bonds"])
        bond_list = self.sections["Bonds"]
        file.write(string.rjust(str(num_bonds), 8) + ' !NBOND\n')
        for index in range(0, num_bonds, 4):
            try:
                bonds = bond_list[index:index+4]
            except IndexError:
                bonds = bond_list[index:-1]
            bond_line = map(lambda bond: string.rjust(str(bond[1]), 8)+string.rjust(str(bond[2]), 8), bonds)
            file.write(''.join(bond_line)+'\n')
        file.close()
    def writePDB(self, filename):
        import string
        atom_format = "%6s%.5s %4s %4s %.4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n"
        p = self.positions
        file = open(filename, 'w')
        for i, atom in enumerate(self.sections["Atoms"]):
            line = ["ATOM  ", str(i+1), 'XXXX', 'TEMP', str(atom.type+1), p[i,0], p[i,1], p[i,2], 0.0, 0.0, str(atom.type)]
            file.write(atom_format%tuple(line))
