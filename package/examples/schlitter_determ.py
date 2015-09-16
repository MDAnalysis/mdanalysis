# -*- coding: utf-8 -*-
from MDAnalysis import *
import Numeric
from numpy import linalg

system = AtomGroup.Universe(".psf", ".dcd")
asel = system.select_atoms(' ( name CA ) ')
# print asel

#---------------------------------
#  making the mass matrix
#---------------------------------
kg = asel.masses * (1.6605388e-27)
masses = Numeric.repeat(kg, 3)
mass_matrix = Numeric.identity(len(masses)) * masses

#--------------------------------
#  Preparing to read the CA position at every 5 steps in the traj
#--------------------------------
skip = 5
num_ts = system.dcd.numframes / skip
num_coor = len(asel) * 3
ca_pos = system.dcd.timeseries(asel, skip=skip, format='fac')

#---------------------------------
# converting angstroms to meters and merging the xyz of timeseries
#---------------------------------
ca = (1e-10) * (Numeric.reshape(ca_pos, (num_ts, -1)))
#print "ca", shape(ca)

#---------------------------------
#  making the covariance matrix
#---------------------------------
ca_avg = Numeric.average(ca)
#print "ca_av", shape(ca_avg)
ca2 = ca - ca_avg[Numeric.NewAxis, :]
#print "ca2", shape(ca2)
ca_cov = Numeric.zeros((num_coor, num_coor), Numeric.Float)
for ts in ca2:
    ca_cov += Numeric.outerproduct(ts, ts)
ca_cov /= num_ts
print "ca_cov", shape(ca_cov)
print "mass_matrix", shape(mass_matrix)

#---------------------------------
#  calculating the entropy
#---------------------------------
hplanck_bar = 6.6260755e-34 / (2 * Numeric.pi)  # J*s
k = (1.380658e-23)  # J/K
Avogadro = 6.0221367e23  # /mol
T = 300  # Kelvin
term = (k * T * exp(2) / pow(hplanck_bar, 2))
print "term =", term
determ = linalg.det((term * Numeric.matrixmultiply(ca_cov, mass_matrix)) + identity(len(mass_matrix)))
print "det = ", determ
S_ho = k / 2 * Avogadro * math.log(determ)
print "S_ho=", S_ho
