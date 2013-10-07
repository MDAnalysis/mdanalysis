#! /opt/local/bin/python

from MDAnalysis import *
from MDAnalysis.analysis.align import *
import MDAnalysis.KDTree.NeighborSearch as KDNS
from MDAnalysis.core.distances import *
from MDAnalysis.analysis.distances import *
from optparse import OptionParser
import numpy as np
import time

def read_rotamer_weights():

	# read in the rotamer weights
	input = open("rotamer-library/R1A_298K_populations.dat",'r')

	rotamerWeight = [0,0]

	while True:

		line = input.readline()

		if not line: break

		rotamerWeight.append(float(line))
	
	return rotamerWeight

def fit_rotamers(rotamers,protein,site_resid,dcdfile):

	# create an ordered list allowing the rotamer to be fitted onto the backbone of the protein
	fittingSelection = (["name C","name CA","name N"],["protein and name C and resid " + str(site_resid),"protein and name CA and resid " + str(site_resid),"protein and name N and resid " + str(site_resid)])

	# fit the rotamer library onto the protein
	rms_fit_trj(rotamers,protein,select=fittingSelection,mass_weighted=True,filename=dcdfile)

	return dcdfile

def find_clashing_rotamers(fitted_rotamers,protein,site_resid,clash_distance):

	# make a KD tree of the protein neighbouring atoms 
	proteinNotSite = protein.selectAtoms("protein and not name H* and not (resid " + str(site_resid) + " or (resid " + str(site_resid-1) + " and (name C or name O)) or (resid " + str(site_resid+1) + " and name N))")
	proteinNotSiteLookup = KDNS.AtomNeighborSearch(proteinNotSite)

	rotamerSel = fitted_rotamers.selectAtoms("not name H*")

	rotamer_clash = [True]
	rotamer_clash_counter=0

	for rotamer in fitted_rotamers.trajectory:
	
		bumps = proteinNotSiteLookup.search_list(rotamerSel,clash_distance)		
		if bumps:
			rotamer_clash.append(True)
			rotamer_clash_counter+=1
		else:	
			rotamer_clash.append(False)

	return (rotamer_clash, rotamer_clash_counter)


parser = OptionParser()
parser.add_option("--pdb", dest="proteinPDB", help="the protein PDB or GRO file")
parser.add_option("--traj", dest="proteinTRAJ", help ="the protein XTC or DCD file (optional)")
parser.add_option("--resid", type=int, nargs=2, dest="residues", help ="the pair of residues to compute DEER distances")
parser.add_option("--discard", dest="discardFrames", type=int, default=0, help ="discard the first N frames")
parser.add_option("--histogramBins", dest="histogramBins", type=float, nargs=3, default=(0,80,1), help ="the start, end and incr distance for the histogram of the distances in Angstroms")
parser.add_option("--clashDistance", dest="clashDistance", type=float, default=2.2, help ="the distance between heavy-atoms of the label and protein within which they are assumed to be clashing (default 2.2 Angstroms)")
parser.add_option("--output", dest="outputFile", default="output.dat", help ="the path and name of the output histogram file")
parser.add_option("--dcdfilename", dest="dcdFileName", help ="the path and stem of the DCD files of the fitted MTSS rotamers")

(options, args) = parser.parse_args()

# load the reference protein structure
if options.proteinPDB and not options.proteinTRAJ:
	proteinStructure = Universe(options.proteinPDB)
elif options.proteinPDB and options.proteinTRAJ:
	proteinStructure = Universe(options.proteinPDB,options.proteinTRAJ)
else:
	print "protein structure and/or trajectory not correctly specified"
	raise SystemExit

if not options.dcdFileName:
	options.dcdFileName = options.outputFile + "-tmp"

startTime = time.time()

# load the XYZ trajectory of the rotamer library
rotamerLibrary = Universe("rotamer-library/rotamer1_R1A_298K.pdb","rotamer-library/rotamer1_R1A_298K.dcd")

# read in the weights of the different rotamers
rotamerWeights = []
rotamerWeights = read_rotamer_weights()

# setup the main lists
distances = []
weights = []

for protein in proteinStructure.trajectory:

	if protein.frame >= options.discardFrames:

		# define the atoms used to fit the rotamers. Note that an ordered list has to be created as the ordering of C CA N is different in both
		# fit the rotamers onto the protein
		fit_rotamers(rotamerLibrary,proteinStructure,options.residues[0],options.dcdFileName+"-1.dcd")
		rotamersSite1 = Universe("rotamer-library/rotamer1_R1A_298K.pdb",options.dcdFileName+"-1.dcd")
		(rotamer1_clash,rotamer1_clash_total) = find_clashing_rotamers(rotamersSite1,proteinStructure,options.residues[0],options.clashDistance)

		fit_rotamers(rotamerLibrary,proteinStructure,options.residues[1],options.dcdFileName+"-2.dcd")
		rotamersSite2 = Universe("rotamer-library/rotamer1_R1A_298K.pdb",options.dcdFileName+"-2.dcd")
		(rotamer2_clash,rotamer2_clash_total) = find_clashing_rotamers(rotamersSite2,proteinStructure,options.residues[1],options.clashDistance)

		# define the atoms to measure the distances between
		rotamer1nitrogen = rotamersSite1.selectAtoms("name N1")
		rotamer2nitrogen = rotamersSite2.selectAtoms("name N1")
	
		# loop over all the rotamers on the first site
		for rotamer1 in rotamersSite1.trajectory:
		
			# only proceed if it isn't the PDB frame and there is no clash	
			if rotamer1.frame > 1 and not rotamer1_clash[rotamer1.frame]:	

				# loop over all the rotamers on the second site
				for rotamer2 in rotamersSite2.trajectory:

					# only proceed if it isn't the PDB frame and there is no clash			
					if rotamer2.frame > 1 and not rotamer2_clash[rotamer2.frame]:	
				
						# measure and record the distance
						(a,b,distance) = MDAnalysis.analysis.distances.dist(rotamer1nitrogen,rotamer2nitrogen)
				
						distances.append(distance[0])

						# create the weights list
						try:
							weight = rotamerWeights[rotamer1.frame] * rotamerWeights[rotamer2.frame]
						except:
							print "oppps", rotamer1.frame, rotamer2.frame

						weights.append(weight)

# check that at least two distances have been measured				
if len(distances) < 2:
	print "no distances found between the spin pair!"
	raise SystemExit

# create an empty list defining the histogram bin edges
bins=[]
i = options.histogramBins[0]
while i <= options.histogramBins[1]:
	bins.append(i)
	i+=options.histogramBins[2]

# use numpy to histogram the distance data, weighted appropriately
(a,b) = np.histogram(distances,weights=weights,normed=True,bins=bins)

OUTPUT = open(options.outputFile + "-" + str(options.residues[0]) + "-" + str(options.residues[1]) + ".dat",'w')

for (i,j) in enumerate(a):
	print >> OUTPUT, "%6.2f %8.3e" %  (((b[i] + b[i+1])/2), j)

OUTPUT.close()				

print "DONE, elapsed time %6i s" % (int(time.time() - startTime))

	