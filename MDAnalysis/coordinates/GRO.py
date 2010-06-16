# GRO reader
# Danny .....
"""
GRO file format
===============

Classes to read and write Gromacs_ GRO coordinate files.

.. _Gromacs: http://www.gromacs.org
"""

import numpy

import base
from base import Timestep

class GROReader(base.Reader):
	format = 'GRO'
	units = {'time': 'ps', 'length': 'nm'}

	def __init__(self,grofilename):
		atoms = []

		atom_iter = 0
		grofile = open(grofilename , 'r')
		for linenum,line in enumerate(grofile):
			try:
				# Store the various parameters
				# resid, resname, name, number = int(line[0:5]) , (line[5:15].split())[0] , (line[5:15].split())[1] , int(line[15:20])  Generate number automatically from an iter
				resid, resname, name, number = int(line[0:5]) , (line[5:15].split())[0] , (line[5:15].split())[1] , int(line[15:20])
				# Not attempting to read type, segid, mass or charge
				type , segid , mass , charge = ("?" , 0 , 0 , 0)
				# Do these work with >= 4 decimal points??
				coords = numpy.array( [ float(line[20:28]) , float(line[28:36]) , float(line[36:44]) ] )

				QueryAtomLine = True
				atom_iter += 1

			except ValueError:
				QueryAtomLine = False
			except IndexError:
				QueryAtomLine = False
			except SyntaxError:
				QueryAtomLine = False
			# If the line can't otherwise be read properly then an error will be raised

			if QueryAtomLine == True:
				# Create an Atom instance
				atom_desc = Atom( atom_iter-1 , name , type , resname , resid , segid , mass , charge)
				# And add it to the list of atoms
				atoms.append(atom_desc)

			elif QueryAtomLine == False:
				# Not currently attempting to pick up comment, number of particles, or box info
				pass

		grofile.close()

		structure = {}
		structure["_atoms"] = atoms
		# Other attributes are not read since they are not included in .gro files
		other_attrs = ["_bonds" , "_angles" , "_dihe" , "_impr" , "_donors" , "_acceptors"]
		for attr in other_attrs:
			structure[attr] = []

		return structure


