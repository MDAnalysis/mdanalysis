from math import *
import numpy
from MDAnalysis import *
 
universe = Universe(...)
xbox, ybox, zbox = universe.dcd.ts.dimensions[:3]

# preparing 1D z-dimension distribution variables 
zmin = int(floor(-zbox/2)) # Hopefully the system is centered at the origin
zmax = int(ceil(zbox/2))
binwidth = 0.5 # Angstroms
numbins = int((zmax - zmin)/binwidth)
 
# setting the charges for the atoms
charge_dist = None
charges = universe.atoms.charges()

# setting the skip for the traj
universe.dcd.skip = 10

# appending atomic charges to a particular z-dimension (z-bin)
nframes = 0
for ts in universe.dcd:
    print ts.frame
    nframes += 1
    # Recenter z coordinates
    z = universe.atoms.coordinates()[:,2]
    z -= numpy.average(z)
    if charge_dist == None: 
        charge_dist, edges = numpy.histogramdd(z, bins=numbins, range=[(zmin,zmax)], weights=charges)
    else:
        newhist, edges = numpy.histogramdd(z, bins=numbins, range=[(zmin,zmax)], weights=charges)
        charge_dist += newhist
 
# Normalize charge distribution
charge_dist /= nframes
 
# midpoints of bins
edge_1d = edges[0]
zpos = 0.5*(edge_1d[:-1]+edge_1d[1:])
 
# Arbitrary zero reference point
zpos -= min(zpos)
 
# Now integrate charge distribution twice (see derivation in Sachs et al paper)
binwidth = zbox/float(numbins)
binvolume = xbox*ybox*binwidth
 
sum_q_binwidth = numpy.add.accumulate(charge_dist*binwidth)
sum_z_q_binwidth = numpy.add.accumulate(zpos*charge_dist*binwidth)
 
phi = zpos*sum_q_binwidth - sum_z_q_binwidth
 
# Enforce periodic boundary and convert units to mV
pbc_factor = sum_z_q_binwidth[-1]/zbox

# converting charge distribution value to mV units
conv_phi = -180951.262*100
mV = (phi + pbc_factor*midpoint)*conv_phi/binvolume
