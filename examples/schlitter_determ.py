from MDAnalysis import *
import Numeric
import numpy
from numpy import *
from numpy import linalg

system = AtomGroup.Universe("/home/denniej0/0_neut/double_system.psf", "/home/denniej0/0_neut/double_2.dcd")
asel = system.selectAtoms(' ( protein or resname HSD and prop z > 0 ) ')
#print asel

#---------------------------------
#  making the mass matrix
#---------------------------------
kg = asel.masses()*(1.6605388e-27)
kg = asel.masses()*(1.66054032e-27)
masses = Numeric.repeat(kg, 3)
mass_matrix = Numeric.identity(len(masses))*masses

skip = 5
num_ts = system._dcd.numframes/skip
num_coor = len(asel)*3
ca_pos = system._dcd.timeseries(asel, skip=skip, format='fac')

#---------------------------------
# angstroms to meters and merging the xyz of timeseries 
#---------------------------------
ca = (1e-10)*(Numeric.reshape(ca_pos, (num_ts, -1)))
#print "ca", shape(ca)

#---------------------------------
#  making the covariance matrix
#---------------------------------
ca_avg = Numeric.average(ca)
#print "ca_av", shape(ca_avg)
ca2 = ca - ca_avg[Numeric.NewAxis,:]
#print "ca2", shape(ca2)
ca_cov = Numeric.zeros((num_coor, num_coor), Numeric.Float)
for ts in ca2:
    ca_cov += Numeric.outerproduct(ts, ts)
ca_cov /= num_ts
print "ca_cov", shape(ca_cov)
print "mass_matrix", shape(mass_matrix)

#---------------------------------
#  calculating the covariance
#---------------------------------
hplanck_bar = 6.6260755e-34/(2*Numeric.pi) # J*s
k =  1.3806580000000001e-23
#k = (1.38e-23)           # J/K
Avogadro = 6.0221367e23   # /mol
T = 300               # Kelvin
term = (k*T*exp(2)/pow(hplanck_bar,2))
print "term =", term
determ = linalg.det((term*Numeric.matrixmultiply(ca_cov,mass_matrix))+identity(len(mass_matrix)))
print "det = ", determ
S_ho = k/2*Avogadro*math.log(determ)
print "S_ho=", S_ho
