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

"""
AMBER PRMTOP topology parser
============================

Reads a  AMBER top file to build the system. It uses atom types,
partial charges and masses from the PRMTOP file.

The format is defined in `PARM parameter/topology file specification`_.
The reader tries to detect if it is a newer (AMBER 12?) file format
by looking for the flag "ATOMIC_NUMBER".

The parser raises a :exc:`TOPParserError` if it fails to read
the topology file.

.. Note::

   The Amber charge is converted to electron charges as used in
   MDAnalysis and other packages. To get back Amber charges, multiply
   by 18.2223.

.. _`PARM parameter/topology file specification`:
   http://ambermd.org/formats.html#topology

.. versionchanged:: 0.7.6
   parses both amber10 and amber12 formats

"""

import MDAnalysis.core
import MDAnalysis.core.units
from MDAnalysis.core import util

class TOPParseError(Exception):
    """Signifies an error during parsing of a Amber PRMTOP file."""
    pass

def parse(filename):
    """Parse Amber PRMTOP topology file *filename*.

    :Returns: MDAnalysis internal *structure* dict.
    """
    formatversion = 10
    with open(filename) as topfile:
        for line in topfile:
            if line.startswith("%FLAG ATOMIC_NUMBER"):
                formatversion = 12
                break
    # Open and check top validity
    ######  Reading header info POINTERS  #################
    with open(filename,'r') as topfile:
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
            from math import ceil
            # Get the number
            num = sys_info[sect_num]
            if data_struc in ["_resname","_bond"]:
                    pass
            else:
                    header = next_line()

            # Now figure out how many lines to read
            numlines = int(ceil(float(num)/per_line))
            #print data_struc, numlines
            if parsefunc == __parsebond_:
                    parsefunc(next_line, atoms_per, data_struc, final_structure, numlines)
            else:
                    parsefunc(next_line, atoms_per, data_struc, structure, numlines)

        sections = {12:
                        [("ATOM_NAME", 1, 20, __parseatoms_, "_name",0),
                         ("CHARGE",1, 5, __parsesection_,"_charge",0),
                         ("ATOMIC_NUMBER", 1, 10, __parsesectionint_,"_skip",0),
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
                         ("LENNARD_JONES_BCOEF", 1, 5, __parseskip_,"_skip",8),
                         #("BONDS_INC_HYDROGEN", 2, 4, __parsebond_, "_bonds",2),
                         #("ANGLES_INC_HYDROGEN", 3, 3, __parsesection_, "_angles"),
                         #("DIHEDRALS_INC_HYDROGEN", 4, 2, __parsesection_, "_dihe"),
                         #("NIMPHI", 4, 2, __parsesection_, "_impr"),
                         #("NDON", 2, 4, __parsesection_,"_donors"),
                         #("NACC", 2, 4, __parsesection_,"_acceptors"),
                         ],
                    10:
                        [("ATOM_NAME", 1, 20, __parseatoms_, "_name",0),
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
                         ("LENNARD_JONES_BCOEF", 1, 5, __parseskip_,"_skip",8),
                         #("BONDS_INC_HYDROGEN", 2, 4, __parsebond_, "_bonds",2),
                         #("ANGLES_INC_HYDROGEN", 3, 3, __parsesection_, "_angles"),
                         #("DIHEDRALS_INC_HYDROGEN", 4, 2, __parsesection_, "_dihe")]
                         #("NIMPHI", 4, 2, __parsesection_, "_impr"),
                         #("NDON", 2, 4, __parsesection_,"_donors"),
                         #("NACC", 2, 4, __parsesection_,"_acceptors")]
                         ],
                    }
        structure = {}
        try:
            for info in sections[formatversion]:
                 parse_sec(info)
        except StopIteration:
            raise TOPParseError("The TOP file didn't contain the minimum required section of ATOM_NAME")
        # Completing info respoint to include all atoms in last resid
        structure["_respoint"].append(sys_info[0])
        structure["_respoint"][-1] = structure["_respoint"][-1] + 1

    atoms = [None,]*sys_info[0]
    from MDAnalysis.core.AtomGroup import Atom
    index = 0
    j = 0
    for i in range(sys_info[0]):
        index += 1
        charge = MDAnalysis.core.units.convert(structure["_charge"][i],'Amber', MDAnalysis.core.flags['charge_unit'])
        if structure["_respoint"][j] <= index < structure["_respoint"][j+1]:
             resid = j + 1
             resname = structure["_resname"][j]
        else:
             j += 1
             resid = j + 1
             resname = structure["_resname"][j]
        mass = structure["_mass"][i]
        atomtype = structure["_atom_type"][i]
        atomname = structure["_name"][i]
        segid = 'SYSTEM'  # does not exist in Amber

        atom_desc = Atom(index-1,atomname,atomtype,resname,resid,segid,mass,charge)
        atoms[i] = atom_desc
    final_structure["_atoms"] = atoms
    final_structure["_numatoms"] = sys_info[0]
    return final_structure

def __parseskip_(lines, atoms_per, attr, structure, numlines):
    section = []
    ##lines()
    while (lines()[:5] != "%FLAG"):
        pass

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
    for i in xrange(numlines):
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
    for i in xrange(numlines):
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
    for i in xrange(numlines):
        l = lines()
        # Subtract 1 from each number to ensure zero-indexing for the atoms
        for j in range(0, len(x.entries)):
            if l[x.entries[j].start:x.entries[j].stop] != '':
                 #print l[x.entries[j].start:x.entries[j].stop]
                 section.append(l[x.entries[j].start:x.entries[j].stop].strip())
            else:
                 continue
    structure[attr] = section
