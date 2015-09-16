# potential profile along the z direction (membrane system)
#
# Reference:
# JN Sachs, PS Crozier, TB Woolf. J Chem Phys 121 (2004), 10847
#   http://dx.doi.org/10.1063/1.1826056

from math import *
import numpy as np
from MDAnalysis import *

# NOTE: set TOPOLOGY and TRAJECTORY to your data: should be a membrane
#       system with the membrane normal parallel to the z-axis
universe = Universe(TOPLOGY, TRAJECTORY)
xbox, ybox, zbox = universe.trajectory.ts.dimensions[:3]

# preparing 1D z-dimension distribution variables
zmin = int(floor(-zbox / 2))  # Hopefully the system is centered at the origin
zmax = int(ceil(zbox / 2))
binwidth = 0.5  # Angstroms
numbins = int((zmax - zmin) / binwidth)

# setting the charges for the atoms
charge_dist = None
charges = universe.atoms.charges

# setting the skip for the traj
universe.trajectory.skip = 10

# appending atomic charges to a particular z-dimension (z-bin)
nframes = 0
for ts in universe.trajectory:
    print ts.frame
    nframes += 1
    # Recenter z coordinates
    z = universe.atoms.coordinates()[:, 2]
    z -= np.average(z)
    if charge_dist is None:
        charge_dist, edges = np.histogramdd(z, bins=numbins, range=[(zmin, zmax)], weights=charges)
    else:
        newhist, edges = np.histogramdd(z, bins=numbins, range=[(zmin, zmax)], weights=charges)
        charge_dist += newhist

# Normalize charge distribution
charge_dist /= nframes

# midpoints of bins
edge_1d = edges[0]
zpos = 0.5 * (edge_1d[:-1] + edge_1d[1:])

# Arbitrary zero reference point
zpos -= min(zpos)

# Now integrate charge distribution twice
# (see derivation in Sachs et al paper, Eq 6)
binwidth = zbox / float(numbins)
binvolume = xbox * ybox * binwidth

sum_q_binwidth = np.add.accumulate(charge_dist * binwidth)
sum_z_q_binwidth = np.add.accumulate(zpos * charge_dist * binwidth)

phi = zpos * sum_q_binwidth - sum_z_q_binwidth

# Enforce periodic boundary and convert units to mV
pbc_factor = sum_z_q_binwidth[-1] / zbox
midpoint = 0.5 * (zpos[0] + zpos[-1])  # shift by average potential at centre of the box

# converting charge distribution value to mV units
# - MDAnalysis charges are in elementary charge units
# - lengths are in Angstrom (1e-10 m)
epsilon0 = 8.854187817620e-12  # As/Vm   vacuum permittivity
e = 1.602176565e-19  # As      elementary charge
conv_phi = -e / epsilon0 * 1e3 * 1e10  # Vm -> mV Angstrom
mV = (phi + pbc_factor * midpoint) * conv_phi / binvolume
