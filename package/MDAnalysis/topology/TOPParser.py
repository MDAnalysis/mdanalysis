# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
AMBER PRMTOP topology parser
============================

Reads a  AMBER top file to build the system. It uses atom types,
partial charges and masses from the PRMTOP file.

The format is defined in `PARM parameter/topology file specification`_.
The reader tries to detect if it is a newer (AMBER 12?) file format
by looking for the flag "ATOMIC_NUMBER".

.. Note::

   The Amber charge is converted to electron charges as used in
   MDAnalysis and other packages. To get back Amber charges, multiply
   by 18.2223.

.. _`PARM parameter/topology file specification`:
   http://ambermd.org/formats.html#topology

Classes
-------

.. autoclass:: TOPParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

from math import ceil

from ..core.AtomGroup import Atom
from ..units import convert
from ..lib.util import openany, FORTRANReader
from ..core import flags
from .base import TopologyReader


class TOPParser(TopologyReader):
    """Reads topology information from an AMBER top file.

    It uses atom types, partial charges and masses from the PRMTOP
    file.

    The format is defined in `PARM parameter/topology file
    specification`_.  The reader tries to detect if it is a newer
    (AMBER 12?) file format by looking for the flag "ATOMIC_NUMBER".

    .. _`PARM parameter/topology file specification`:
       http://ambermd.org/formats.html#topology

   .. versionchanged:: 0.7.6
      parses both amber10 and amber12 formats

    """

    def parse(self):
        """Parse Amber PRMTOP topology file *filename*.

        :Returns: MDAnalysis internal *structure* dict.
        """
        formatversion = 10
        with openany(self.filename) as topfile:
            for line in topfile:
                if line.startswith("%FLAG ATOMIC_NUMBER"):
                    formatversion = 12
                    break
        if formatversion == 12:
            sections = [
                ("ATOM_NAME", 1, 20, self._parseatoms, "_name", 0),
                ("CHARGE", 1, 5, self._parsesection, "_charge", 0),
                ("ATOMIC_NUMBER", 1, 10, self._parsesectionint, "_skip", 0),
                ("MASS", 1, 5, self._parsesection, "_mass", 0),
                ("ATOM_TYPE_INDEX", 1, 10, self._parsesectionint, "_atom_type", 0),
                ("NUMBER_EXCLUDED_ATOMS", 1, 10, self._parseskip, "_skip", 8),
                ("NONBONDED_PARM_INDEX", 1, 10, self._parseskip, "_skip", 8),
                ("RESIDUE_LABEL", 1, 20, self._parseatoms, "_resname", 11),
                ("RESIDUE_POINTER", 2, 10, self._parsesectionint, "_respoint", 11),
            ]
                #("BOND_FORCE_CONSTANT", 1, 5, self._parseskip,"_skip",8),
                #("BOND_EQUIL_VALUE", 1, 5, self._parseskip,"_skip",8),
                #("ANGLE_FORCE_CONSTANT", 1, 5, self._parseskip,"_skip",8),
                #("ANGLE_EQUIL_VALUE", 1, 5, self._parseskip,"_skip",8),
                #("DIHEDRAL_FORCE_CONSTANT", 1, 5, self._parseskip,"_skip",8),
                #("DIHEDRAL_PERIODICITY", 1, 5, self._parseskip,"_skip",8),
                #("DIHEDRAL_PHASE", 1, 5, self._parseskip,"_skip",8),
                #("SOLTY", 1, 5, self._parseskip,"_skip",8),
                #("LENNARD_JONES_ACOEF", 1, 5, self._parseskip,"_skip",8),
                #("LENNARD_JONES_BCOEF", 1, 5, self._parseskip,"_skip",8),
                #("BONDS_INC_HYDROGEN", 2, 4, self._parsebond, "_bonds",2),
                #("ANGLES_INC_HYDROGEN", 3, 3, self._parsesection, "_angles"),
                #("DIHEDRALS_INC_HYDROGEN", 4, 2, self._parsesection, "_dihe"),
                #("NIMPHI", 4, 2, self._parsesection, "_impr"),
                #("NDON", 2, 4, self._parsesection,"_donors"),
                #("NACC", 2, 4, self._parsesection,"_acceptors"),
        elif formatversion == 10:
            sections = [
                ("ATOM_NAME", 1, 20, self._parseatoms, "_name", 0),
                ("CHARGE", 1, 5, self._parsesection, "_charge", 0),
                ("MASS", 1, 5, self._parsesection, "_mass", 0),
                ("ATOM_TYPE_INDEX", 1, 10, self._parsesectionint, "_atom_type", 0),
                ("NUMBER_EXCLUDED_ATOMS", 1, 10, self._parseskip, "_skip", 8),
                ("NONBONDED_PARM_INDEX", 1, 10, self._parseskip, "_skip", 8),
                ("RESIDUE_LABEL", 1, 20, self._parseatoms, "_resname", 11),
                ("RESIDUE_POINTER", 2, 10, self._parsesectionint, "_respoint", 11),
            ]
                #("BOND_FORCE_CONSTANT", 1, 5, self._parseskip,"_skip",8),
                #("BOND_EQUIL_VALUE", 1, 5, self._parseskip,"_skip",8),
                #("ANGLE_FORCE_CONSTANT", 1, 5, self._parseskip,"_skip",8),
                #("ANGLE_EQUIL_VALUE", 1, 5, self._parseskip,"_skip",8),
                #("DIHEDRAL_FORCE_CONSTANT", 1, 5, self._parseskip,"_skip",8),
                #("DIHEDRAL_PERIODICITY", 1, 5, self._parseskip,"_skip",8),
                #("DIHEDRAL_PHASE", 1, 5, self._parseskip,"_skip",8),
                #("SOLTY", 1, 5, self._parseskip,"_skip",8),
                #("LENNARD_JONES_ACOEF", 1, 5, self._parseskip,"_skip",8),
                #("LENNARD_JONES_BCOEF", 1, 5, self._parseskip,"_skip",8),
                #("BONDS_INC_HYDROGEN", 2, 4, self._parsebond, "_bonds",2),
                #("ANGLES_INC_HYDROGEN", 3, 3, self._parsesection, "_angles"),
                #("DIHEDRALS_INC_HYDROGEN", 4, 2, self._parsesection, "_dihe")]
                #("NIMPHI", 4, 2, self._parsesection, "_impr"),
                #("NDON", 2, 4, self._parsesection,"_donors"),
                #("NACC", 2, 4, self._parsesection,"_acceptors")]

        # Open and check top validity
        # Reading header info POINTERS
        with openany(self.filename) as topfile:
            next_line = topfile.next
            header = next_line()
            if header[:3] != "%VE":
                raise ValueError("{0} is not a valid TOP file. %VE Missing in header".format(topfile))
            title = next_line().split()
            if not (title[1] == "TITLE"):
                raise ValueError("{0} is not a valid TOP file. 'TITLE' missing in header".format(topfile))
            while header[:14] != '%FLAG POINTERS':
                header = next_line()
            header = next_line()

            topremarks = [next_line().strip() for i in xrange(4)]
            sys_info = [int(k) for i in topremarks for k in i.split()]

            structure = {}
            final_structure = {}

            try:
                for info in sections:
                    self._parse_sec(sys_info, info, next_line,
                                    structure, final_structure)
            except StopIteration:
                raise ValueError("The TOP file didn't contain the minimum"
                                 " required section of ATOM_NAME")
            # Completing info respoint to include all atoms in last resid
            structure["_respoint"].append(sys_info[0])
            structure["_respoint"][-1] = structure["_respoint"][-1] + 1

        atoms = [None, ]*sys_info[0]

        j = 0
        segid = "SYSTEM"
        for i in range(sys_info[0]):
            charge = convert(structure["_charge"][i],
                             'Amber',
                             flags['charge_unit'])
            if structure["_respoint"][j] <= i+1 < structure["_respoint"][j+1]:
                resid = j + 1
                resname = structure["_resname"][j]
            else:
                j += 1
                resid = j + 1
                resname = structure["_resname"][j]
            mass = structure["_mass"][i]
            atomtype = structure["_atom_type"][i]
            atomname = structure["_name"][i]
            #segid = 'SYSTEM'  # does not exist in Amber

            atoms[i] = Atom(i, atomname, atomtype, resname, resid,
                            segid, mass, charge, universe=self._u)
        final_structure["atoms"] = atoms
        final_structure["_numatoms"] = sys_info[0]
        return final_structure

    def _parse_sec(self, sys_info, section_info, next_line, structure, final_structure):
        desc, atoms_per, per_line, parsefunc, data_struc, sect_num = section_info
        # Get the number
        num = sys_info[sect_num]
        if data_struc in ["_resname", "_bond"]:
            pass
        else:
            header = next_line()

        # Now figure out how many lines to read
        numlines = int(ceil(float(num)/per_line))
        #print data_struc, numlines
        if parsefunc == self._parsebond:
            parsefunc(next_line, atoms_per, data_struc, final_structure, numlines)
        else:
            parsefunc(next_line, atoms_per, data_struc, structure, numlines)

    def _parseskip(self, lines, atoms_per, attr, structure, numlines):
        while (lines()[:5] != "%FLAG"):
            pass

    def _parsebond(self, lines, atoms_per, attr, structure, numlines):
        section = []  # [None,]*numlines
        for i in xrange(numlines):
            l = lines()
            # Subtract 1 from each number to ensure zero-indexing for the atoms
            f = map(int, l.split())
            fields = [a-1 for a in f]
            for j in range(0, len(fields), atoms_per):
                section.append(tuple(fields[j:j+atoms_per]))
        structure[attr] = section

    def _parsesectionint(self, lines, atoms_per, attr, structure, numlines):
        section = []  # [None,]*numlines
        y = lines().strip("%FORMAT(")
        y.strip(")")
        x = FORTRANReader(y)
        for i in xrange(numlines):
            l = lines()
            # Subtract 1 from each number to ensure zero-indexing for the atoms
            try:
                for j in xrange(len(x.entries)):
                    section.append(int(l[x.entries[j].start:x.entries[j].stop].strip()))
            except:
                continue
        structure[attr] = section

    def _parsesection(self, lines, atoms_per, attr, structure, numlines):
        section = []  # [None,]*numlines
        y = lines().strip("%FORMAT(")
        y.strip(")")
        x = FORTRANReader(y)
        for i in xrange(numlines):
            l = lines()
            # Subtract 1 from each number to ensure zero-indexing for the atoms
            try:
                for j in range(0, len(x.entries)):
                    section.append(float(l[x.entries[j].start:x.entries[j].stop].strip()))
            except:
                continue
        structure[attr] = section

    def _parseatoms(self, lines, atoms_per, attr, structure, numlines):
        section = []  # [None,]*numlines
        y = lines().strip("%FORMAT(")
        y.strip(")")
        x = FORTRANReader(y)
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
