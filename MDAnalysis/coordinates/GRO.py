# $Id: GRO.py 101 2010-07-22 13:19:06Z Danny.Parton $
"""
GRO file format
===============

Classes to read and write Gromacs_ GRO_ coordinate files; see the notes on the
`GRO format`_ which includes a conversion routine for the box.

.. _Gromacs: http://www.gromacs.org
.. _GRO: http://manual.gromacs.org/current/online/gro.html
.. _GRO format: http://chembytes.wikidot.com/g-grofile
"""

import numpy

import MDAnalysis
import base
import MDAnalysis.core.util as util
from MDAnalysis.coordinates.core import triclinic_box, triclinic_vectors

from copy import deepcopy

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
		x = self._unitcell[[0,3,4]]
		y = self._unitcell[[5,1,6]]
		z = self._unitcell[[7,8,2]]  # this ordering is correct! (checked it, OB)
		return triclinic_box(x,y,z)

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
			# see Timestep.dimensions() above for format (!)
			self.ts._unitcell[:3] = unitcell
		elif len(unitcell) == 9:
			self.ts._unitcell[:] = unitcell   # fill all
		else:   # or maybe raise an error for wrong format??			
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

class GROWriterStandard(base.Writer):
	"""GRO Writer that conforms to the Trajectory API.

	.. Note:: The precision is hard coded to three decimal places.

	.. Warning:: This class will replace GROWriter.
	"""
	
        format = 'GRO'
        units = {'time': None, 'length': 'nm'}

	#: Output format, see http://chembytes.wikidot.com/g-grofile
	fmt_velocities = "%5s%-5s%5s%5s%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
	fmt_no_velocities = "%5s%-5s%5s%5s%8.3f%8.3f%8.3f"

        def __init__(self,filename):
		"""Set up a GROWriter with a precision of 3 decimal places.

		:Arguments:
		   *filename* 
                      output filename
		"""
                self.filename = util.filename(filename,ext='gro')

	def convert_dimensions_to_unitcell(self, ts):
		"""Read dimensions from timestep *ts* and return appropriate unitcell"""
		return self.convert_pos_to_native(triclinic_vectors(ts.dimensions))

	def write(self, selection, frame=None):
		"""Write selection at current trajectory frame to file.
		
		:Arguments:
		  selection
                      MDAnalysis AtomGroup (selection or Universe.atoms)
                :Keywords:
		  frame             
                      optionally move to frame number *frame*

		The GRO format only allows 5 digits for resid and atom
		number. If these number become larger than 99,999 then this
		routine will chop off the leading digits.
		"""
		# write() method that complies with the Trajectory API
		u = selection.universe
		if frame is not None:            
			u.trajectory[frame]  # advance to frame
		else:
			try:
				frame = u.trajectory.ts.frame
			except AttributeError:
				frame = 1   # should catch cases when we are analyzing a single GRO (?)

                output_gro = open(self.filename , 'w')
                # Header
                output_gro.write('Written by MDAnalysis\n')
                output_gro.write("%5d\n" % len(selection))

                # Atom descriptions and coords
		coordinates = selection.coordinates()
		self.convert_pos_to_native(coordinates)   # Convert back to nm from Angstroms, in-place !
		
		for atom_index,atom in enumerate(selection):
			c = coordinates[atom_index]
                        output_line = self.fmt_no_velocities % \
			    (str(atom.resid)[-5:],     # truncate highest digits on overflow
			     atom.resname.strip(), 
			     atom.name.strip(),
			     str(atom.number+1)[-5:],  # number (1-based), truncate highest digits on overflow
			     c[0], c[1], c[2],         # coords - outputted with 3 d.p.
			     )
                        output_gro.write( output_line + '\n' )

                # Footer: box dimensions
		box = self.convert_dimensions_to_unitcell(u.trajectory.ts)
		if numpy.all(u.trajectory.ts.dimensions[3:] == [90.,90.,90.]):
			# orthorhombic cell, only lengths along axes needed in gro
			output_gro.write("%10.5f%10.5f%10.5f\n" % (box[0,0],box[1,1],box[2,2]))
		else:
			# full output
			output_gro.write("%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n" % 
					 (box[0,0],box[1,1],box[2,2],
					  box[0,1],box[0,2],
					  box[1,0],box[1,2],
					  box[2,0],box[2,1]))	
                output_gro.close()


class GROWriter(base.Writer):
        """
The universe option can be either a Universe or an AtomGroup (e.g. a selectAtoms() selection).
The coordinates to be output are taken from the current Timestep (unless a different Timestep is supplied when calling write() )."""
        format = 'GRO'
        units = {'time': None, 'length': 'nm'}
        def __init__(self,filename,universe=None,ndec=3):
		"""Set up a GROWriter.

		:Arguments:
		   *ndec* 
                     the number of decimal places for the
                     coordinates. Currently only works with 3 (default).
		"""
		import warnings
		warnings.warn("Deprecated GROWriter; will be changed in next release to behave in the same "
			      "way as CRDWriter and PrimitivePDBWriter (as detailed in the Trajectory API).",
			      category=DeprecationWarning)

                self.filename = util.filename(filename,ext='gro')
		if isinstance(universe , MDAnalysis.Universe):
                	self.universe = universe
			self.atom_indices = self.universe.atoms.indices()
		# If it's an AtomGroup, then take the atom indices, and then store the Universe which that AtomGroup belongs to
		elif isinstance(universe , MDAnalysis.AtomGroup.AtomGroup):
			self.atom_indices = universe.indices()
			self.universe = universe.universe
		elif universe == None:
			raise Exception, 'Must supply a Universe or AtomGroup'
		else:
			raise Exception, 'Must supply a Universe or AtomGroup'

        def write(self,ts=None,ndec=3):
                """ndec is the number of decimal places for the coordinates. Currently only works with 3 (default).
		If ts=None then we try to get one from the Universe
		"""
                if ts is None:
                        try:
                                ts = self.universe.trajectory.ts
                        except:
                                raise Exception, "Can't find ts in universe.trajectory"
                output_gro = open(self.filename , 'w')
                # Header
                output_gro.write('Written by MDAnalysis\n')
                output_gro.write((' ' * (5 - len(str(len(self.atom_indices))))) + str(len(self.atom_indices)) + '\n')
                # Atom descriptions and coords
		for atom_index in self.atom_indices:
			atom = self.universe.atoms[atom_index]
                        # resid
                        output_line = (' ' * (5 - len(str(atom.resid)))) + str(atom.resid)
                        # resname + atomname
                        output_line += atom.resname + (' ' * (10 - len(atom.resname) - len(atom.name))) + atom.name
                        # number (1-based)
                        output_line += (' ' * (5 - len(str(atom.number+1)))) + str(atom.number+1)
                        # coords - outputted with 3 d.p.
			# These come from a supplied Timestep, if it is supplied, otherwise from the Universe.trajectory.ts
			coords = deepcopy(ts[atom.number])
                	# Convert back to nm from Angstroms
			self.convert_pos_to_native(coords)   # in-place !
                        coords = [ '%.3f' % coords[0] , '%.3f' % coords[1] , '%.3f' % coords[2] ]
                        output_line += (' ' * (4 - len(coords[0].split('.')[0]))) + coords[0]
                        output_line += (' ' * (4 - len(coords[1].split('.')[0]))) + coords[1]
                        output_line += (' ' * (4 - len(coords[2].split('.')[0]))) + coords[2]
                        # Output the line
                        output_gro.write( output_line + '\n' )

                # Footer: box dimensions
		dims = ts.dimensions[:3]
                # Convert back to nm from Angstroms
		self.convert_pos_to_native(dims)   # in-place !
		output_gro.write( (' ' * (4 - len(str(dims[0]).split('.')[0]))) + '%.5f' % dims[0] + (' ' * (4 - len(str(dims[1]).split('.')[0]))) + '%.5f' % dims[1] + (' ' * (4 - len(str(dims[2]).split('.')[0]))) + '%.5f' % dims[2] + '\n')
                output_gro.close()



