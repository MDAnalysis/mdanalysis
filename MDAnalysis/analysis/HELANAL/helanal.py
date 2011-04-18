#!/usr/bin/env python
"""
HELANAL algorithm 
=================

Copyright (c) 2009 Benjamin Hall <benjamin.hall@bioch.ox.ac.uk>

Rewritten for clarity.
"""

from __future__ import with_statement

#from CartesianToolkit import *
import numpy

# replace CartesianToolkit with NumPy
# - vecscale(v,a) --> a*c
# - vecadd(a,b) --> a+b
# - vecsub(a,b) --> a-b
try:
	from numpy import rad2deg, deg2rad   # numpy 1.3+
except ImportError:
	def rad2deg(x):             # no need for the numpy out=[] argument 
		return 180.0*x/numpy.pi
	def deg2rad(x):             # no need for the numpy out=[] argument 
		return x*numpy.pi/180.0

def center(coordinates):
	"""Return the geometric center (centroid) of the coordinates.

	Coordinates must be "list of cartesians", i.e. a Nx3 array.
	"""
	return numpy.mean(coordinates, axis=0)

def veclength(v):
	"""Length of vector *v*."""
	# note: this is 3 times faster than numpy.linalg.norm
	return numpy.sqrt(numpy.dot(v,v))

vecscaler = numpy.dot

def vecnorm(a):
	"""Return a/|a|"""
	return a/veclength(a)

def vecangle(a,b):
	"""Angle between two vectors *a* and *b* in degrees.
	
	If one of the lengths is 0 then the angle is returned as 0
	(instead of `nan`).
	"""
	angle = numpy.arccos(numpy.dot(a,b) / (veclength(a)*veclength(b)))
	if numpy.isnan(angle):
		return 0.0
	return rad2deg(angle)

def vecdist(a,b):
	"""Return |a-b|"""
	return veclength(a-b)

veccross = numpy.cross

def wrapangle(angle):
	"""Wrap angle (in radians) to be within -pi < angle =< pi"""
	if angle > numpy.pi:
		angle -= 2*numpy.pi
	elif angle =< -numpy.pi:
		angle += 2*numpy.pi
	return angle
	
	

#import math, numpy, sys
import math,sys,os
try:
	import gmxdump_if as xtcio
except:
	import xtcio
import ndxio


try:
	import psyco
except:
	pass

def helanal_trajectory(xtcfile,pdbfile,ndxfile,start,end,begin,finish,matrix_filename,origin_pdbfile,summary_filename,screw_filename,tilt_filename,fitted_tilt_filename,prefix):
	pdbdata = open(pdbfile)
	[indices,pdblength] = get_backbone_indices(pdbdata)
	pdbdata.close()
	trajectory  = xtcio.read_xtc(xtcfile)
	#origin_pdbfile = "origin.pdb"
	#matrix_filename = "bending_matrix.txt"
	if prefix != None:
		matrix_filename = prefix + matrix_filename
		origin_pdbfile = prefix + origin_pdbfile
		summary_filename = prefix + summary_filename
		screw_filename = prefix + screw_filename
		tilt_filename = prefix + tilt_filename
		fitted_tilt_filename = prefix + fitted_tilt_filename
	backup_file(matrix_filename)
	backup_file(origin_pdbfile)
	backup_file(summary_filename)
	backup_file(screw_filename)
	backup_file(tilt_filename)
	backup_file(fitted_tilt_filename)
	global_height = []
	global_twist = []
	global_rnou = []
	global_bending = []
	global_bending_matrix = []
	global_tilt = []
	global_fitted_tilts = []
	global_screw = []
	frame = 0
	if start != None and end != None:
			print "Analysing from residue", start, "to", end
	elif start != None and end == None:
			print "Analysing from residue", start, "to the C termini"
	elif start == None and end != None:
			print "Analysing from the N termini to", end
	#while(trajectory.next_frame()):
	for cart in trajectory:
		if begin != None:
			if trajectory.time < begin:
				continue
		if finish != None:
			if trajectory.time > finish:
				break
		if pdblength != trajectory.natoms:
			print "Different numbers of atoms in the pdb and xtc files (", pdblength , "vs" , trajectory.natoms ,")! Exiting"
			exit()
		trajectory.convert_coordinates()
		ca_positions = index_cartesians(indices,trajectory.cartesian)
		max_length = len(ca_positions)
		if start != None and end != None:
			#print "Analysing from residue", start, "to", end
			ca_positions = ca_positions[(start-1):(end)]
		elif start != None and end == None:
			#print "Analysing from residue", start, "to the C termini"
			ca_positions = ca_positions[(start-1):]
		elif start == None and end != None:
			#print "Analysing from the N termini to", end
			ca_positions = ca_positions[:(end)]
		if len(global_height) == 0:
			measured_length = len(ca_positions)
			print "Analysing", measured_length, "/", max_length, "residues"
		
		[twist,bending_angles,height,rnou,origins,local_helix_axes,local_screw_angles] = main_loop(ca_positions)

				
		origin_pdb(origins,origin_pdbfile)
		
		#calculate local bending matrix( it is looking at all i, j combinations)
		if len(global_bending_matrix) == 0:
			global_bending_matrix = [ [ [] for item in local_helix_axes] for item in local_helix_axes ]

		for i in range(len(local_helix_axes)):
			for j in range(i+1,len(local_helix_axes)):
				angle = math.acos(vecscaler(local_helix_axes[i],local_helix_axes[j]))*180/math.pi
				global_bending_matrix[i][j].append(angle)
				#global_bending_matrix[j][i].append(angle)
			#global_bending_matrix[i][i].append(0.)

		
		fit_vector, fit_tilt = vector_of_best_fit(origins)
		global_height += height
		global_twist  += twist
		global_rnou   += rnou
		#global_screw.append(local_screw_angles)
		global_fitted_tilts.append(fit_tilt*180/math.pi)
		
		#print out rotations across the helix to a file
		with open(screw_filename,"a") as rot_output:
			print >> rot_output, frame,
			for rotation in local_screw_angles:
				print >> rot_output, rotation,
			print >> rot_output, ""
		
		with open(tilt_filename,"a") as tilt_output:
			print >> tilt_output, frame,
			for tilt in local_helix_axes:
				print >> tilt_output, vecangle(tilt,[0,0,1])*180/math.pi,
			print >> tilt_output, ""
		
		with open(fitted_tilt_filename,"a") as tilt_output:
			print >> tilt_output, frame, fit_tilt*180/math.pi
		
		if len(global_bending) == 0:
			global_bending = [ [] for item in bending_angles ]
			#global_tilt = [ [] for item in local_helix_axes ]
		for store,tmp in zip(global_bending,bending_angles): store.append(tmp)
		#for store,tmp in zip(global_tilt,local_helix_axes): store.append(vecangle(tmp,[0,0,1]))
		
		#simple ticker
		formated_time = "%20.1f" % trajectory.time
		frame += 1
		formated_frame ="%10d" % frame
			
		print '\r',formated_time,' ps' , formated_frame,
		sys.stdout.flush()

	print '\nComplete'
	[twist_mean, twist_sd, twist_abdev] = stats(global_twist)
	[height_mean, height_sd, height_abdev] = stats(global_height)
	[rnou_mean, rnou_sd, rnou_abdev] = stats(global_rnou)
	[ftilt_mean, ftilt_sd, ftilt_abdev] = stats(global_fitted_tilts)
	
	bending_statistics = [ stats(item) for item in global_bending]
	#tilt_statistics =    [ stats(item) for item in global_tilt]

	bending_statistics_matrix = [[stats(col) for col in row] for row in global_bending_matrix]
	mat_output = open(matrix_filename,'w')
	print >> mat_output, "Mean"
	#[ [ print >> mat_output, "%8.3f\t", % col[0] for col in row ] and print '' for row in bending_statistics_matrix ]
	for row in bending_statistics_matrix:
		for col in row:
			formatted_angle = "%6.1f" % col[0]
			print >> mat_output, formatted_angle,
		print >> mat_output, ''
	print >> mat_output, "\nSD"
	#[ [ print >> mat_output, "%8.3f\t", % col[0] for col in row ] and print '' for row in bending_statistics_matrix ]
	for row in bending_statistics_matrix:
		for col in row:
			formatted_angle = "%6.1f" % col[1]
			print >> mat_output, formatted_angle,
		print >> mat_output, ''
	print >> mat_output, "\nABDEV"
	#[ [ print >> mat_output, "%8.3f\t", % col[0] for col in row ] and print '' for row in bending_statistics_matrix ]
	for row in bending_statistics_matrix:
		for col in row:
			formatted_angle = "%6.1f" % col[2]
			print >> mat_output, formatted_angle,
		print >> mat_output, ''
	
	mat_output.close()

	print "Height:", height_mean, "SD", height_sd, "ABDEV", height_abdev, '(nm)'
	print "Twist:", twist_mean, "SD", twist_sd, "ABDEV", twist_abdev
	print "Residues/turn:", rnou_mean, "SD", rnou_sd, "ABDEV", rnou_abdev
	print "Fitted tilt:", ftilt_mean, "SD", ftilt_sd, "ABDEV", ftilt_abdev
	print "Local bending angles:"
	residue_statistics = zip(*bending_statistics)
	measure_names = ["Mean ","SD   ","ABDEV"]
	print "ResID",
	if start == None:
		for item in range(4,len(residue_statistics[0])+4):
			output = "%8d" % item
			print output,
	else:
		for item in range(start+3,len(residue_statistics[0])+start+3):
			output = "%8d" % item
			print output,
	print ""
	for measure,name in zip(residue_statistics,measure_names):
		print name,
		for residue in measure:
			output = "%8.1f" % residue
			print output,
		print ''



	summary_output = open(summary_filename,'w')
	print >> summary_output, "Height:", height_mean, "SD", height_sd, "ABDEV", height_abdev, '(nm)'
	print >> summary_output, "Twist:", twist_mean, "SD", twist_sd, "ABDEV", twist_abdev
	print >> summary_output, "Residues/turn:", rnou_mean, "SD", rnou_sd, "ABDEV", rnou_abdev
	print >> summary_output, "Local bending angles:"
	residue_statistics = zip(*bending_statistics)
	measure_names = ["Mean ","SD   ","ABDEV"]
	print >> summary_output,  "ResID",
	if start == None:
		for item in range(4,len(residue_statistics[0])+4):
			output = "%8d" % item
			print >> summary_output,  output,
	else:
		for item in range(start+3,len(residue_statistics[0])+start+3):
			output = "%8d" % item
			print >> summary_output, output,
	print >> summary_output, ""

	for measure,name in zip(residue_statistics,measure_names):
		print >> summary_output, name,
		for residue in measure:
			output = "%8.1f" % residue
			print >> summary_output, output,
		print >> summary_output, ''

	summary_output.close()

def tilt_correct(number):
	'''Changes an angle (in degrees) so that it is between 0 and 90'''
	if number < 90:
		return number
	else:
		return 180 - number

def backup_file(filename):
	if os.path.exists(filename):
		target_name = "#" + filename
		failure = True
		if not os.path.exists(target_name):
			os.rename(filename,target_name)
			failure = False
		else:
			for i in range(20):
				alt_target_name = target_name + "." + str(i)
				if os.path.exists(alt_target_name):
					continue
				else:
					os.rename(filename,alt_target_name)
					failure = False
					break
		if failure:
			print "Too many backups. Clean up and try again"
			exit()

def stats(some_list):
	if len(some_list) == 0:
		return [0,0,0]
	list_mean = mean(some_list)
	list_sd = sample_sd(some_list,list_mean)
	list_abdev = mean_abs_dev(some_list,list_mean)
	return [list_mean, list_sd, list_abdev]

def helanal_main(pdbfile,start,end):
	pdbdata = open(pdbfile)
	positions = get_backbone_positions(pdbdata)
	max_length = len(positions)

	if start != None and end != None:
			print "Analysing from residue", start, "to", end
			positions = positions[(start-1):(end)]
	elif start != None and end == None:
			print "Analysing from residue", start, "to the C termini"
			positions = positions[(start-1):]
	elif start == None and end != None:
			print "Analysing from the N termini to", end
			positions = positions[:(end)]

	measured_length = len(positions)
	print "Analysing", measured_length, "/", max_length, "residues"
	[twist,bending_angles,height,rnou,origins,local_helix_axes,local_screw_angles] = main_loop(positions)

	#TESTED- origins are correct
	#print current_origin
	#print origins
	
	max_angle = max(bending_angles)
	mean_angle = mean(bending_angles)
	#sd calculated using n-1 to replicate original fortran- assumes a limited sample so uses the sample standard deviation
	sd_angle = sample_sd(bending_angles,mean_angle)
	mean_absolute_deviation_angle = mean_abs_dev(bending_angles,mean_angle)
	#TESTED- stats correct
	#print max_angle, mean_angle, sd_angle, mean_absolute_deviation_angle
	
	#calculate local bending matrix(now it is looking at all i, j combinations)
	for i in local_helix_axes:
		for j in local_helix_axes:
			if i == j:
				angle = 0.
			else:
				angle = math.acos(vecscaler(i,j))*180/math.pi
			string_angle = "%6.0f\t" % angle
			#print string_angle,
		#print ''
		#TESTED- local bending matrix!
	
	#Average helical parameters
	mean_twist = mean(twist)
	sd_twist = sample_sd(twist,mean_twist)
	abdev_twist = mean_abs_dev(twist,mean_twist)
	#TESTED-average twists
	#print mean_twist, sd_twist, abdev_twist
	mean_rnou = mean(rnou)
	sd_rnou = sample_sd(rnou,mean_rnou)
	abdev_rnou = mean_abs_dev(rnou,mean_rnou)
	#TESTED-average residues per turn
	#print mean_rnou, sd_rnou, abdev_rnou
	mean_height = mean(height)
	sd_height = sample_sd(height,mean_height)
	abdev_height = mean_abs_dev(height,mean_height)	
	#TESTED- average rises
	
	print "Height:", mean_height, sd_height, abdev_height
	print "Twist:", mean_twist, sd_twist, abdev_twist
	print "Residues/turn:", mean_rnou, sd_rnou, abdev_rnou
	#print mean_height, sd_height, abdev_height
	print "Local bending angles:"
	for angle in bending_angles:
		output = "%8.1f\t" % angle
		print output,
	print ''
	
	#calculate best fit vector and tilt of said vector
	fit_vector, fit_tilt = vector_of_best_fit(origins)
	print "Best fit tilt =", fit_tilt
	print "Rotation Angles from 1 to n-1"
	for item in local_screw_angles:
		print round(item,1),
	
def origin_pdb(origins,pdbfile):	#special- print origins to pdb, assumes the need to convert nm->A
	with open(pdbfile,'a') as output:
		i=1
		for item in origins:
			tmp = "ATOM    %3d  CA  ALA   %3d    %8.3f%8.3f%8.3f  1.00  0.00" % (i,i,item[0]*10,item[1]*10,item[2]*10)
			print >> output, tmp
			i += 1
		print >> output, "TER\nENDMDL"

def main_loop(positions):
	twist = []
	rnou = []
	height = []
	origins = [[0.,0.,0.] for item in positions[:-2]]
	local_helix_axes = []
	location_rotation_vectors = []
	for i in range(len(positions)-3):
		vec12 = positions[i+1]- positions[i]
		vec23 = positions[i+2] -positions[i+1]
		vec34 = positions[i+3] -positions[i+2]
		
		dv13 = vec12 - vec23
		dv24 = vec23 - vec34
		
		#direction of the local helix axis
		current_uloc = vecnorm(veccross(dv13,dv24))
		local_helix_axes.append(current_uloc)
		
		#TESTED- Axes correct
		#print current_uloc

		dmag = veclength(dv13)
		emag = veclength(dv24)
		
		costheta = vecscaler(dv13,dv24)/(dmag*emag)
		#rnou is the number of residues per turn
		current_twist = math.acos(costheta)
		twist.append(current_twist*180/math.pi)
		rnou.append(2*math.pi/current_twist)
		#radius of local helix cylinder radmag
		
		costheta1 = 1.0 - costheta
		radmag = (dmag*emag)**0.5/(2*costheta1)
		
		#Height of local helix cylinder
		current_height = vecscaler(vec23,current_uloc)
		height.append(current_height)
		#TESTED- Twists etc correct
		#print current_twist*180/math.pi, 2*math.pi/current_twist, height
		
		dv13 = vecnorm(dv13)
		dv24 = vecnorm(dv24)
		
		#record local rotation
		location_rotation_vectors.append(dv13)
		
		rad = [radmag * item for item in dv13]
		current_origin = [(item[0] - item[1]) for item in zip(positions[i+1],rad)]
		origins[i] = current_origin
		
		#TESTED- origins are correct
		#print current_origin
		
		rad = [radmag * item for item in dv24]
		current_origin = [(item[0] - item[1]) for item in zip(positions[i+2],rad)]
		origins[i+1] = current_origin
	#Record final rotation vector
	location_rotation_vectors.append(dv24)
	
	#local bending angles (eg i > i+3, i+3 > i+6)
	
	bending_angles = [0 for item in range(len(local_helix_axes)-3)]
	for axis in range(len(local_helix_axes)-3):
		angle = math.acos(vecscaler(local_helix_axes[axis],local_helix_axes[axis+3]))*180/math.pi
		bending_angles[axis] = angle
		#TESTED- angles are correct
		#print angle
		
	local_screw_angles = []
	#Calculate rotation angles for (+1) to (n-1)
	fit_vector, fit_tilt = vector_of_best_fit(origins)
	for item in location_rotation_vectors:
		local_screw_tmp = rotation_angle(fit_vector,[0,0,1],item)*180/math.pi
		#print local_screw_tmp
		local_screw_angles.append(local_screw_tmp)
		
	return [twist,bending_angles,height,rnou,origins,local_helix_axes,local_screw_angles]

def rotation_angle(helix_vector,axis_vector,rotation_vector):
	reference_vector = veccross(veccross(helix_vector,axis_vector),helix_vector)
	second_reference_vector = veccross(axis_vector,helix_vector)
	screw_angle = vecangle(reference_vector,rotation_vector)
	alt_screw_angle = vecangle(second_reference_vector,rotation_vector)
	updown = veccross(reference_vector,rotation_vector)
	
	if screw_angle > math.pi/4 and screw_angle < 3*math.pi/4:
		pass
	else:
		if screw_angle < math.pi/4 and alt_screw_angle < math.pi/2:
			screw_angle = math.pi/2 - alt_screw_angle
		elif screw_angle < math.pi/4 and alt_screw_angle > math.pi/2:
			screw_angle = alt_screw_angle - math.pi/2
		elif screw_angle > 3* math.pi/4 and alt_screw_angle < math.pi/2:
			screw_angle = math.pi/2 + alt_screw_angle
		elif screw_angle > 3*math.pi/4 and alt_screw_angle > math.pi/2:
			screw_angle = 3*math.pi/2 - alt_screw_angle
		else:
			print "\nBig Screw Up"
	
	if veclength(updown) == 0:
		#vector is at 0 or 180
		print "\nPROBLEM"
		
	helix_dot_rehelix = vecangle(updown,helix_vector)
	
	#if ( helix_dot_rehelix < math.pi/2 and helix_dot_rehelix >= 0 )or helix_dot_rehelix <-math.pi/2:
	if ( helix_dot_rehelix < math.pi/2 and helix_dot_rehelix > -math.pi/2) or (helix_dot_rehelix > 3*math.pi/2):
			screw_angle = 0 - screw_angle
			#print "Same     ", helix_dot_rehelix*180/math.pi
	else:
			#print "Different", helix_dot_rehelix*180/math.pi
			pass
	
	return screw_angle

def vector_of_best_fit(origins):
	origins = numpy.asarray(origins)
	centroids = center(origins)
	M = numpy.matrix(origins - centroids)
	A = M.transpose() * M
	u,s,vh=numpy.linalg.linalg.svd(A)
	vector = vh[0].tolist()[0]
	#Correct vector to face towards first residues
	rough_helix = origins[0] - centroids
	agreement = vecangle(rough_helix,vector)
	if agreement < math.pi/2 and agreement > -math.pi/2:
		pass
	else:
		vector = vh[0]*-1
		vector = vector.tolist()[0]
	best_fit_tilt = vecangle(vector,[0,0,1])
	return vector,best_fit_tilt

def fit(origins):
	#Subroutine to fir plane, circle and line to local helix origins
	#INCOMPLETE
	
	#Not sure exactly what these represent
	[x,y,z] = [0., 0., 0.]
	[x2,y2,z2] = [0., 0., 0.]
	[xy,xz,yz] = [0., 0., 0.]
	
	matp = [ [ 0. for i in range(3)] for i in range(3) ] 
	for item in origins:
		x2 += item[0]**2
		y2 += item[1]**2
		z2 += item[2]**2
		xy += item[0]*item[1]
		xz += item[0]*item[2]
		yz += item[1]*item[2]
		x  += item[0]
		y  += item[1]
		z  += item[2]
	matp[0][0] += x2
	matp[0][1] += xy
	matp[1][0] += xy
	matp[0][2] += xz
	matp[2][0] += xz
	matp[1][1] += y2
	matp[2][2] += z2
	matp[1][2] += yz
	matp[2][1] += yz
	
	matp = numpy.matrix(matp)
	pmat = matp.I
	

def index_cartesians(index_list,frame):
	return [frame[item] for item in index_list]	

from numpy import mean
def sample_sd(a, dummy):
	return numpy.std(a, ddof=1)
def mean_abs_dev(a,mean_a=None):
	if mean_a is None:
		mean_a = mean(a)
	return mean(numpy.fabs(a - mean_a))
		
def get_backbone_positions(pdbdata):
	positions = []
	for line in pdbdata:
		atom_name = line[11:16].strip()
		if atom_name =="CA" or (atom_name[:1] == "B" or atom_name[:2] == "0B"): #Bondini, atomistic or MARTINI
			positions.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
	return positions

def get_backbone_indices(pdbdata):
	indices = []
	current_index = 0
	pdblength = 0
	for line in pdbdata:
		if line[:4] != 'ATOM':
			continue
		pdblength += 1
		atom_name = line[11:16].strip()
		if atom_name =="CA" or ((atom_name[:1] == "B" or atom_name[1:2] == "B") and atom_name[:1] != "C" ): #Bondini, atomistic or MARTINI
			indices.append(current_index)
		current_index +=1
	print "Length", len(indices)
	return indices, pdblength

def getoptions():
	from optparse import OptionParser
	parser = OptionParser()
	
	parser.add_option("-p", "--pdb", dest="pdb_file", default=None,
	         help="PDB file", metavar="FILE")
	parser.add_option("-x", "--xtc", dest="xtc_file", default=None,
	         help="XTC file", metavar="FILE")
	parser.add_option("-n", "--ndx", dest="ndx_file", default=None,
	         help="NDX file", metavar="FILE")
	parser.add_option("-m", "--matrix", dest="mat_file", default="bending_matrix.txt",
	         help="Output file- bending matrix", metavar="FILE")
	parser.add_option("-o", "--origin", dest="ori_file", default="origin.pdb",
	         help="Output file- origin pdb file", metavar="FILE")
	parser.add_option("-r", "--roundup", dest="summary_file", default="summary.txt",
	         help="Output file- all of the basic data", metavar="FILE")
	parser.add_option("-c", "--screw", dest="screw_file", default="screw.xvg",
	         help="Output file- screw rotations of individual residues from 2 to n-1", metavar="FILE")
	parser.add_option("-t", "--tilt", dest="screw_tilt", default="local_tilt.xvg",
	         help="Output file- local tilts of individual residues from 2 to n-1", metavar="FILE")
	parser.add_option("-i", "--fit-tilt", dest="fit_tilt", default="fit_tilt.xvg",
	         help="Output file- tilt of line of best fit applied to origin axes", metavar="FILE")
	parser.add_option("--prefix", dest="prefix", default=None,
	         help="Prefix to add to all output files", metavar="FILE")
	parser.add_option("-s", "--start",
	         dest="start", default=None,type="int",
	         help="Start residue")
	parser.add_option("-e", "--end",
	         dest="end", default=None,type="int",
	         help="End residue")
	parser.add_option("-b", "--begin",
	         dest="begin", default=None,type="int",
	         help="Begin analysing from time (ps)")
	parser.add_option("-f", "--finish",
	         dest="finish", default=None,type="int",
	         help="Stop analysing after time (ps)")
	'''parser.add_option("-v","--verbose",
			 dest="verbose", default=False, action="store_true",
			 help="Generate a file for each residue's information")'''

	(options, args) = parser.parse_args()
	return options
	
	

if __name__ == "__main__":
	import sys
	options = getoptions()
	if options.pdb_file == None:
		print "Needs a PDB file"
		exit()
	if options.xtc_file == None:
		print "No xtc file- working on the pdb alone"
		helanal_main(options.pdb_file,options.start,options.end)
		exit()
	helanal_trajectory(options.xtc_file,options.pdb_file,options.ndx_file,options.start,options.end,options.begin,options.finish,options.mat_file,options.ori_file,options.summary_file,options.screw_file,options.screw_tilt,options.fit_tilt,options.prefix)
