# GRO parser
# Danny

from MDAnalysis.core.AtomGroup import Atom
import numpy

class GROParseError(Exception):
	pass

def parse(grofilename):
	### General checks on validity of file
	# Check file ends in .gro
	if grofilename[-4:] != ".gro":
		raise GROParseError("%(grofilename)s file does not end in .gro")

	### Read through .gro file
	atom_iter = 0
	atoms = []
	grofile = open(grofilename , "r")
	for linenum,line in enumerate(grofile):
		query_atm_line = False
		try:
			# Store the various parameters
			# Not attempting to read velocities
			# resid, resname, name, number = int(line[0:5]) , (line[5:15].split())[0] , (line[5:15].split())[1] , int(line[15:20])  Generate number automatically from an iter
			resid, resname, name, number = int(line[0:5]) , (line[5:15].split())[0] , (line[5:15].split())[1] , int(line[15:20])
			# Not attempting to read type, segid, mass or charge
			type , segid , mass , charge = (0 , "SYSTEM" , 12.0 , 0)
			# Or coords, as they can be read by coordinates.GRO
			# Do these work with >= 4 decimal points??
			#coords = numpy.array( [ float(line[20:28]) , float(line[28:36]) , float(line[36:44]) ] )
			query_atm_line = True

		# Not currently doing anything with other lines
		except (ValueError, IndexError, SyntaxError):
			if linenum == 0:
				# Header comment
				#hdr_cmt = line
				pass
			elif linenum == 1:
				# Header: number of particles
				#hdr_np = int(line)
			# A bit dodgy; should find a better way of locating the box_vectors line
				pass
			else:
				#ftr_box = line
		# If the line can't otherwise be read properly, then this probably indicates a problem with the gro line, and an error will be raised
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


	grofile.close()
	structure = {}
	structure["_atoms"] = atoms
	# Other attributes are not read since they are not included in .gro files
	other_attrs = ["_bonds" , "_angles" , "_dihe" , "_impr" , "_donors" , "_acceptors"]
	for attr in other_attrs:
		structure[attr] = []

	return structure

