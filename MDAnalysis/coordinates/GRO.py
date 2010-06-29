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
	units = {'time': None, 'length': 'nm'}

	def __init__(self,grofilename):
		self.grofilename = grofilename
		self.filename = self.grofilename
		coords_list = []

		grofile = open(grofilename , 'r')

		# Read first two lines to get number of atoms
		grofile.readline()
		total_atnums = int(grofile.readline())
		grofile.seek(0)
		for linenum,line in enumerate(grofile):
			# Should work with any precision
			if linenum not in (0,1,total_atnums+2):
				coords_list.append( numpy.array( map( float , line[20:].split()[0:3] ) ) )
			# Unit cell footer
			elif linenum == total_atnums+2:
				unitcell = numpy.array( map( float , line.split() ) )

		grofile.close()

		self.numatoms = len(coords_list)
		coords_list = numpy.array(coords_list)
		self.ts = Timestep(coords_list)
		# ts._unitcell layout is [A, alpha, B, beta, gamma, C]
		# Should be setting the ts.dimensions property, rather than ts._unitcell??
		self.ts._unitcell[0] = unitcell[0]
		self.ts._unitcell[2] = unitcell[1]
		self.ts._unitcell[5] = unitcell[2]
		self.numframes = 1
		self.fixed = 0
		self.skip = 1
		self.periodic = False
		self.delta = 0
		self.skip_timestep = 1
		self.units = {'time': None, 'length': 'nm'}

	def __len__(self):
		return self.numframes
	def __iter__(self):
		yield self.ts  # Just a single frame
		raise StopIteration
	def __getitem__(self, frame):
		if frame != 0:
			raise IndexError('GROReader can only read a single frame at index 0')
		return self.ts
	def _read_next_timestep(self):
		raise Exception, "GROReader can only read a single frame"



