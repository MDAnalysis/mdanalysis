# $Id: TRJ.py 101 2011-01-22 13:19:06Z Elizabeth Denning $
"""
TRJ file format
===============

Classes to read and write Amber_ TRJ_ coordinate files; see the notes on the
`TRJ format`_ which includes a conversion routine for the box.

.. _Amber: http://ambermd.org
.. _TRJ format: http://ambermd.org/formats.html#trajectory
"""

import numpy

import MDAnalysis
import base
import MDAnalysis.core.util as util
#import MDAnalysis.core.AtomGroup._psf as psf
from MDAnalysis.coordinates.core import triclinic_box, triclinic_vectors
from base import Timestep

from copy import deepcopy

class TRJReader(base.Reader):
	format = 'TRJ'
	units = {'time': 'ps', 'length': 'Angstroms'}
	_Timestep = Timestep

	def __init__(self,trjfilename): 
		
		self.filename = trjfilename
		self.__numatoms = None
		self.__numframes = None
		
		self.trjfile = file(self.filename, 'r')
		self.fixed = 0
		self.skip = 1
		self.periodic = False
		self.ts = Timestep(self.numatoms)

		numatoms = self.numatoms 
		x = util.FORTRANReader("10F8.3")
		numlines = numatoms / len(x.entries)
		
		self._read_next_timestep()
		#self.dimensions = triclinic_box(x,y,z)

	def _read_next_timestep(self):
		ts = self.ts
		if self.trjfile is None:
		        self.open_trajectory()
		
		coords_list = []
		numatoms = self.numatoms
		
		# Read first line to get number of atoms
		for linenum,line in enumerate(self.trjfile):
			if len(coords_list) == numatoms * 3:
				counter = 0 
				for i in range(0,numatoms*3,3):
					self.ts._x[counter] = coords_list[i]
					self.ts._y[counter] = coords_list[i+1]
					self.ts._z[counter] = coords_list[i+2]
					counter += 1
				unitcell = numpy.array( map( float, line.split())) 
				self.ts._unitcell = numpy.zeros(3, dtype=numpy.float32)   # TRJ has 3 entries
				coords_list = []
				if len(unitcell) == 3:
					# special case: a b c --> (a 0 0) (b 0 0) (c 0 0)
					self.ts._unitcell[:3] = unitcell
				elif len(unitcell) == 9:
					self.ts._unitcell[:] = unitcell   # fill all
				else:   
					import warnings
					warnings.warn("TRJ unitcell has neither 3 nor 9 entries --- might be wrong.")
					self.ts._unitcell[:len(unitcell)] = unitcell   
				self.ts.frame += 1 
				#self._read_next_frame(ts._x, ts._y, ts._z, ts._unitcell, self.skip)
				return self.ts
			else:
				for j in line.split():
					coords_list.append(float(j))
		
		self.fixed = 0
		self.skip = 1
		self.periodic = False
		self.delta = 0
		self.skip_timestep = 1
		#self.dimensions = triclinic_box
	@property
	def numframes(self):
	        if not self.__numframes is None:   # return cached value
		        return self.__numframes
		try:
			self.__numframes = self._read_trj_numframes(self.filename)
		except IOError:
			return 0
		else:
			return self.__numframes
	
	@property
	def dimensions(self):
	        """unitcell dimensions (A, B, C, alpha, beta, gamma)
		- A, B, C are the lengths of the primitive cell vectors e1, e2, e3
		- alpha = angle(e1, e2)
		- beta = angle(e1, e3)
		- gamma = angle(e2, e3)
		"""
		# Layout of unitcell is [X, Y, Z] with the primitive cell vectors
		x = self.ts._unitcell[0]
		y = self.ts._unitcell[1]
		z = self.ts._unitcell[2]
		return triclinic_box(x,y,z)

	def _read_trj_numatoms(self, filename):
		# this assumes that this is only called once at startup and that the filestream is already open
		# read the first line
		n = self.trjfile.readline() ##_psf.items()[0] #._numatoms

		#self.close_trajectory()
		# need to check type of n
		return int(n)
	
	def _read_trj_numframes(self, filename):
	        self._reopen()
		# the number of lines in the XYZ file will be 2 greater than the number of atoms 
		linesPerFrame = self.numatoms * 3. / 10. 

		counter = 0
		# step through the file (assuming xyzfile has an iterator)
		for i in self.trjfile:
			counter = counter + 1
		self.close_trajectory()
		# need to check this is an integer!
		numframes = int(counter/linesPerFrame)
		return numframes

	@property
	def numatoms(self):
	        if not self.__numatoms is None:   # return cached value
		        return self.__numatoms
		try:
			self.__numatoms = self._read_trj_numatoms(self.filename)
		except IOError:
			return 0
		else:
			return self.__numatoms

	def __del__(self):
		if not self.trjfile is None:
			self.close_trajectory()

	def __len__(self):
		return self.numframes

	def _reopen(self):
	        self.close_trajectory()
		self.open_trajectory()

    	def open_trajectory(self):
		self.trjfile = file(self.filename, 'r')
		self.trjfile.readline()
		# reset ts
		ts = self.ts
		ts.status = 1
		ts.frame = 0
		ts.step = 0
		ts.time = 0
		return self.trjfile

	def close_trajectory(self):
		"""Close trj trajectory file if it was open."""
		if self.trjfile is None:
			return
		self.trjfile.close()
		self.trjfile = None


	def rewind(self):
		self._reopen()
		self.next()

	
	def __iter__(self):
		self._reopen()
		#yield self.ts
		#raise StopIteration
		self.ts.frame = 0  # start at 0 so that the first frame becomes 1
		while True:
			try:
				yield self._read_next_timestep()
			except EOFError:
				self.close_trajectory()
				raise StopIteration
