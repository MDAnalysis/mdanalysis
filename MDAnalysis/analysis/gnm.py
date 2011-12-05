#Analyse a trajectory using elastic network models, following the approach of Hall et al (JACS 2007)
#Ben Hall (benjamin.a.hall@ucl.ac.uk) is to blame
#Copyright 2011; Consider under GPL v2 or later

from __future__ import with_statement

import numpy
from numpy import linalg

from MDAnalysis import MissingDataWarning
from MDAnalysis.core.AtomGroup import AtomGroup
import MDAnalysis.KDTree.NeighborSearch as NS
from MDAnalysis.core.util import norm, angle, parse_residue

import warnings
import logging
logger = logging.getLogger('MDAnalysis.analysis.GNM')

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

def generate_grid(positions,cutoff):
   [x, y, z] = zip(*positions)
   high_x = max(x)
   high_y = max(y)
   high_z = max(z)
   low_x = min(x)
   low_y = min(y)
   low_z = min(z)
   natoms = len(positions)
   #Ok now generate a list with 3 dimensions representing boxes in x, y and z
   grid = [[[[] for i in range(int((high_z - low_z) /cutoff)+1)] for j in range(int((high_y - low_y) /cutoff)+1)] for k in range(int((high_x - low_x) /cutoff)+1)]
   res_positions = []
   for i in range(natoms):
      x_pos = int( (positions[i][0] - low_x)     /cutoff )
      y_pos = int( (positions[i][1] - low_y) /cutoff )
      z_pos = int( (positions[i][2] - low_z) /cutoff )
      grid[x_pos][y_pos][z_pos].append(i)
      res_positions.append([x_pos,y_pos,z_pos])
   return (res_positions,grid,low_x,low_y,low_z)


def generate_kirchoff(positions,cutoff):
   natoms = len(positions)
   matrix = numpy.zeros((natoms,natoms),"float")
   [res_positions,grid,low_x,low_y,low_z] = generate_grid(positions,cutoff)
   icounter = 0
   for icounter in range(natoms):
      #find neighbours from the grid
      neighbour_atoms = []
      for x in (-1,0,1):
         if (res_positions[icounter][0]+x) >= 0 and (res_positions[icounter][0]+x) < len(grid):
            for y in (-1,0,1):
               if (res_positions[icounter][1]+y) >= 0 and (res_positions[icounter][1]+y) < len(grid[0]):
                  for z in (-1,0,1):
                     if (res_positions[icounter][2]+z) >= 0 and (res_positions[icounter][2]+z) < len(grid[0][0]):
                        neighbour_atoms += grid[res_positions[icounter][0]+x][res_positions[icounter][1]+y][res_positions[icounter][2]+z]
      #for jcounter in range(icounter+1,natoms):
      for jcounter in neighbour_atoms:
         if jcounter > icounter and ((positions[icounter][0] -  positions[jcounter][0])**2 + (positions[icounter][1] -  positions[jcounter][1])**2 + (positions[icounter][2] -  positions[jcounter][2])**2) <= cutoff**2:
            matrix[icounter][jcounter] = -1.0
            matrix[jcounter][icounter] = -1.0
            matrix[icounter][icounter] = matrix[icounter][icounter] + 1
            matrix[jcounter][jcounter] = matrix[jcounter][jcounter] + 1
   return (matrix)

def order_list(w):
   ordered = list(w)
   unordered = list(w)
   ordered.sort()
   list_map = {}
   for i in range(len(w)):
      list_map[i] = unordered.index(ordered[i])
   return list_map

def generate_output(w,v,outputobject,time, matrix, nmodes=2, ReportVector=None,counter=0):
   list_map = order_list(w)
   print round(time), w[list_map[1]]
   if ReportVector:
	with open(ReportVector,"a") as oup:
		for item in enumerate(v[list_map[1]]):
			print >> oup, "", counter, time, item[0]+1, w[list_map[1]], item[1]
   outputobject.append((time, w[list_map[1]], v[list_map[1]]))
   #outputobject.append((time, [ w[list_map[i]] for i in range(nmodes) ], [ v[list_map[i]] for i in range(nmodes) ] ))

class GNMAnalysis:
	def __init__(self, universe, selection='name CA', cutoff=7.0, ReportVector=None):
		print "Not recommended by anyone for anything. Frankly, I'm not even sure why you're running it"
		self.u = universe
		self.selection = selection
		self.cutoff = cutoff
		self.results = []  # final result
		self._timesteps = None  # time for each frame
		self.ReportVector=ReportVector

	def run(self, skip=1):
		"""Analyze trajectory and produce timeseries.

		Returns GNM results per frame:

		  results = [(time,eigenvalues,eigenvectors,kirchoff_matrix),(time,eigenvalues,eigenvectors,kirchoff_matrix)... ]


		"""
		logger.info("GNM analysis: starting")
		counter = 0

		self.timeseries = []
		self._timesteps = []

		try:
			self.u.trajectory.time
			def _get_timestep():
				return self.u.trajectory.time
			logger.debug("GNM analysis is recording time step")
		except NotImplementedError:
			# chained reader or xyz(?) cannot do time yet
			def _get_timestep():
				return self.u.trajectory.frame
			logger.warn("GNM analysis is recording frame number instead of time step")

		#Ok, set up the important parameters
		ca = self.u.selectAtoms(self.selection)

		for ts in self.u.trajectory:
			if counter % skip != 0:
				counter += 1
				continue
			counter += 1
			frame = ts.frame
			timestep = _get_timestep()
			self._timesteps.append(timestep)
			ca_positions = ca.coordinates()
			matrix = generate_kirchoff(ca_positions, self.cutoff)
			try:
				[u,w,v] = linalg.svd(matrix)
			except:
				print "\nFrame skip at", timestep,"(SVD failed to converge). Cutoff", self.cutoff
				continue
			#Save the results somewhere useful in some useful format. Usefully.
			generate_output(w,v,self.results,timestep,matrix,ReportVector=self.ReportVector,counter=counter)

