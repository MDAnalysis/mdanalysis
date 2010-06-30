# GRO reader
# Danny .....
"""
GRO file format
===============

Classes to read and write Gromacs_ GRO_ coordinate files.

.. _Gromacs: http://www.gromacs.org
.. _GRO: http://manual.gromacs.org/current/online/gro.html
"""

import numpy

import MDAnalysis.core
import base

class Timestep(base.Timestep):
	@property
	def dimensions(self):
		"""unitcell dimensions (A, B, C, alpha, beta, gamma)

		GRO:
		8.00170   8.00170   5.65806   0.00000   0.00000   0.00000   0.00000   4.00085   4.00085

		PDB:
		CRYST1   80.017   80.017   80.017  60.00  60.00  90.00 P 1           1	
	
		XTC: c.trajectory.ts._unitcell
		array([[ 80.00515747,   0.        ,   0.        ],
		       [  0.        ,  80.00515747,   0.        ],
		       [ 40.00257874,  40.00257874,  56.57218552]], dtype=float32)
		"""
		# unit cell line (from http://manual.gromacs.org/current/online/gro.html)
		# v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
		# 0     1     2      3     4     5     6    7     8
		from MDAnalysis.coordinates.core import _veclength, _angle
		x = self._unitcell[[0,3,4]]
		y = self._unitcell[[5,1,6]]
		z = self._unitcell[[7,8,2]]  # this ordering is correct! (checked it, OB)
		A, B, C = [_veclength(v) for v in x,y,z]
		alpha =  _angle(x,y)
		beta  =  _angle(x,z)
		gamma =  _angle(y,z)
		return numpy.array([A,B,C,alpha,beta,gamma])

class GROReader(base.Reader):
	format = 'GRO'
	units = {'time': None, 'length': 'nm'}
	_Timestep = Timestep

	def __init__(self,grofilename,convert_units=None):
		self.grofilename = grofilename
		self.filename = self.grofilename
		if convert_units is None:
			convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
		self.convert_units = convert_units  # convert length and time to base units

		coords_list = []

		grofile = open(grofilename , 'r')
		try:
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
		finally:
			grofile.close()

		self.numatoms = len(coords_list)
		coords_list = numpy.array(coords_list)
		self.ts = self._Timestep(coords_list)
		# ts._unitcell layout is format dependent; Timestep.dimensions does the conversion
		# behind the scene
		self.ts._unitcell = numpy.zeros(9, dtype=numpy.float32)   # GRO has 9 entries 
		if len(unitcell) == 3:
			# special case: a b c --> (a 0 0) (b 0 0) (c 0 0)
			# see dimensions() below for format (!)
			self.ts._unitcell[[0,1,7]] = unitcell
		elif len(unitcell) == 9:
			self.ts._unitcell[:] = unitcell   # fill all
		else:
			import warnings
			warnings.warn("GRO unitcell has neither 3 nor 9 entries --- might be wrong.")
			self.ts._unitcell[:len(unitcell)] = unitcell   # fill linearly ... not sure about this
		if self.convert_units:
			self.convert_pos_from_native(self.ts._pos)             # in-place !
			self.convert_pos_from_native(self.ts._unitcell)        # in-place ! (all are lengths)
		self.numframes = 1
		self.fixed = 0
		self.skip = 1
		self.periodic = False
		self.delta = 0
		self.skip_timestep = 1

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



