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
GRO topology parser
===================

Read a list of atoms from a GROMOS/Gromacs GRO coordinate file to build a basic topology.

Atom types and masses are guessed.
"""

from MDAnalysis.core.AtomGroup import Atom
from MDAnalysis.topology.core import guess_atom_type, guess_atom_mass, guess_atom_charge
import numpy

class GROParseError(Exception):
        pass

def parse(filename):
        """Parse GRO file *filename* and return the dict `structure`.

        Only reads the list of atoms.

        :Returns: MDAnalysis internal *structure* dict

        .. SeeAlso:: The *structure* dict is defined in
                     :func:`MDAnalysis.topology.PSFParser.parse`.

        """
        ### Read through .gro file
        atom_iter = 0
        atoms = []
        with open(filename , "r") as grofile:
                for linenum,line in enumerate(grofile):
                        query_atm_line = False
                        try:
                                # Store the various parameters
                                # Not attempting to read velocities
                                resid, resname, name, number = int(line[0:5]) , (line[5:15].split())[0] , (line[5:15].split())[1] , int(line[15:20])
                                # guess based on atom name
                                type = guess_atom_type(name)
                                mass = guess_atom_mass(name)
                                charge = guess_atom_charge(name)
                                segid = "SYSTEM"
                                # ignore coords, as they can be read by coordinates.GRO
                                query_atm_line = True

                        # Not currently doing anything with other lines
                        except (ValueError, IndexError):
                                if linenum == 0:
                                        # Header comment
                                        #hdr_cmt = line
                                        pass
                                elif linenum == 1:
                                        # Header: number of particles
                                        #hdr_np = int(line)
                                        # A bit dodgy; should find a better way
                                        # of locating the box_vectors line
                                        pass
                                else:
                                        #ftr_box = line If the line can't
                                        # otherwise be read properly, then this
                                        # probably indicates a problem with the
                                        # gro line, and an error will be raised
                                        pass
                        except:
                                print "Couldn't read the following line of the .gro file:\n%s" % line
                                raise

                        # If it's an atom line (either with velocities or without) then add it to the list of atoms
                        if query_atm_line == True:
                                # Create an Atom instance
                                # Just use the atom_iter (counting from 0) rather than the number in the .gro file (which wraps at 99999)
                                atom_desc = Atom( atom_iter , name , type , resname , resid , segid , mass , charge)
                                # And add it to the list of atoms
                                atoms.append(atom_desc)
                                atom_iter += 1
        structure = {}
        structure["_atoms"] = atoms
        # Other attributes are not read since they are not included in .gro files
        other_attrs = ["_bonds" , "_angles" , "_dihe" , "_impr" , "_donors" , "_acceptors"]
        for attr in other_attrs:
                structure[attr] = []

        return structure

