# $Id: TRJ.py 101 2011-01-22 13:19:06Z Elizabeth Denning $
"""
Amber TRJ file format
=====================

Classes to read Amber_ TRJ coordinate files as defined in `Amber TRJ
format`_. It is also possible to directly read *bzip2* or *gzip*
compressed files.

Amber trajectories are recognised by the suffix '.trj' or '.mdcrd'
(possibly with an additional '.gz' or '.bz2').

.. Note:: Support for Amber is *experimental* and feedback and
   contributions are highly appreciated. Use the `Issue Tracker`_ or
   get in touch on the `MDAnalysis mailinglist`_.

.. _Amber: http://ambermd.org
.. _Amber TRJ format: http://ambermd.org/formats.html#trajectory
.. _Issue Tracker: http://code.google.com/p/mdanalysis/issues/list
.. _MDAnalysis mailinglist: http://groups.google.com/group/mdnalysis-discussion

Units
-----

* lengths in Angstrom
* time in ps (but see below)


Limitations
-----------

* Reads the ASCII `Amber TRJ format`_, but not the `Amber netcdf`_ one
  (yet).

* Periodic boxes are only stored as box lengths A, B, C in an Amber
  trajectory; the reader always assumes that these are orthorhombic
  boxes.

* The trajectory does not contain time information so we simply set
  the time step to 1 ps (or the user could provide it as kwarg *delta*)

* No direct access of frames is implemented, only iteration through
  the trajectory.

* If the trajectory contains exactly *one* atom then it is always
  assumed to be non-periodic (for technical reasons).

.. _Amber netcdf: http://ambermd.org/netcdf/nctraj.html
"""

import numpy

import MDAnalysis
import base
import MDAnalysis.core.util as util

class Timestep(base.Timestep):
	"""Amber trajectory Timestep"""
	@property
	def dimensions(self):
	        """unitcell dimensions (A, B, C, alpha, beta, gamma)
		- A, B, C are the lengths of the primitive cell vectors e1, e2, e3
		- alpha = angle(e1, e2)
		- beta = angle(e1, e3)
		- gamma = angle(e2, e3)

		Note that the Amber trajectory only contains box lengths A,B,C; we assume
		an orthorhombic box.
		"""
		# Layout of unitcell is [A,B,C,90,90,90] with the primitive cell vectors
		return self._unitcell

class TRJReader(base.Reader):
	format = 'TRJ'
	units = {'time': 'ps', 'length': 'Angstroms'}
	_Timestep = Timestep

	# TODO: implement random access via seek
	#       - compute size of frame
	#       - compute seek offset & go
	#       - check that this works for files >2GB

	def __init__(self, trjfilename, numatoms=None, **kwargs): 
		# amber trj REQUIRES the number of atoms from the topology
		if numatoms is None:
			raise ValueError("Amber TRJ reader REQUIRES the numatoms keyword")
		self.filename = trjfilename
		self.__numatoms = numatoms
		self.__numframes = None
		
		self.trjfile = None  # have _read_next_timestep() open it properly!
		self.fixed = 0
		self.skip = 1
		self.skip_timestep = 1   # always 1 for trj at the moment
		self.delta = kwargs.pop("delta", 1.0)   # can set delta manually, default is 1ps
		self.ts = Timestep(self.numatoms)

		# FORMAT(10F8.3)  (X(i), Y(i), Z(i), i=1,NATOM)
		self.default_line_parser = util.FORTRANReader("10F8.3")
		self.lines_per_frame = int(numpy.ceil(3.0 * self.numatoms / len(self.default_line_parser)))
		# The last line per frame might have fewer than 10 
		# We determine right away what parser we need for the last
		# line because it will be the same for all frames.
		last_per_line = 3 * self.numatoms % len(self.default_line_parser)
		self.last_line_parser = util.FORTRANReader("%dF8.3" % last_per_line)

		# FORMAT(10F8.3)  BOX(1), BOX(2), BOX(3)
		# is this always on a separate line??
		self.box_line_parser = util.FORTRANReader("3F8.3")

		# Now check for box
		self.periodic = False
		self._detect_amber_box()

		# open file, read first frame
		self._read_next_timestep()

	def _read_next_timestep(self):
		# FORMAT(10F8.3)  (X(i), Y(i), Z(i), i=1,NATOM)
		ts = self.ts
		if self.trjfile is None:
		        self.open_trajectory()

		# Read coordinat frame:
		#coordinates = numpy.zeros(3*self.numatoms, dtype=numpy.float32)
		_coords = []
		for number,line in enumerate(self.trjfile):
			try:
				_coords.extend(self.default_line_parser.read(line))
			except ValueError:
				# less than 10 entries on the line:
				_coords.extend(self.last_line_parser.read(line))
			if number == self.lines_per_frame - 1:
				break

		# Read box information
		if self.periodic:
			line = self.trjfile.next()
			box = self.box_line_parser.read(line)
			ts._unitcell[:3] = numpy.array(box, dtype=numpy.float32)
			ts._unitcell[3:] = [90.,90.,90.]  # assumed

		# probably slow ... could be optimized by storing the coordinates in X,Y,Z
		# lists or directly filling the array; the array/reshape is not good
		# because it creates an intermediate array
		ts._pos[:] = numpy.array(_coords).reshape(self.numatoms, 3)
		ts.frame += 1
		return ts

	def _detect_amber_box(self):
		"""Detecting a box in a Amber trajectory

		Rewind trajectory and check for potential box data
		after the first frame.

		Set :attr:`TRJReader.periodic` to ``True`` if box was
		found, ``False`` otherwise.

		Only run at the beginning as it *rewinds* the trajctory.

		 - see if there's data after the atoms have been read that looks
		   like::

		     FORMAT(10F8.3)  BOX(1), BOX(2), BOX(3)
		     BOX    : size of periodic box

                 - this WILL fail if we have exactly 1 atom in the trajectory because
		   there's no way to distinguish the coordinates from the box
		   so for 1 atom we always assume no box
		 XXX: needs a Timestep that knows about Amber unitcells!
		"""
		if self.numatoms == 1:
			# for 1 atom we cannot detect the box with the current approach
			self.periodic = False   # see _read_next_timestep()!
			wmsg = "Trajectory contains a single atom: assuming periodic=False"
			warnings.warn(wmsg)
			return False

		self._reopen()
		self.periodic = False      # make sure that only coordinates are read
		self._read_next_timestep()
		ts = self.ts
		line = self.trjfile.next()
		try:
			box = self.box_line_parser.read(line)
			if len(box) == 3:
				ts._unitcell[:3] = box
				ts._unitcell[3:] = [90.,90.,90.]  # assumed
				self.periodic = True
			else:
				raise ValueError  # break into except clause
		except ValueError:
			# 1) line has FEWER than 3 values (box_line_parser.read()
			#    raises ValueError)
			# 2) not exactly 3 values
			# --> no box, but next frame so we need to back up
			self.periodic = False
			ts._unitcell = numpy.zeros(6, numpy.float32)
		self.close()
		return self.periodic

		
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
	
	def _read_trj_numatoms(self, filename):
		raise NotImplementedError("It is not possible to relaibly deduce NATOMS from Amber trj files")
	
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
		self.trjfile, filename = util.anyopen(self.filename, 'r')
		self.header = self.trjfile.readline()  # ignore first line
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
