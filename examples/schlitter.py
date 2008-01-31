
# This requires a trajectory where the protein has already been C-alpha RMS fitted

import LinearAlgebra
from MDAnalysis import *
import Numeric

#sys = AtomGroup.System("psf file", "dcd file")
system = AtomGroup.System("data/beta_op_adp-vmd.psf", "data/step1.dcd")
asel = system.selectAtoms('segid F and backbone and not type O')

masses = Numeric.repeat(asel.masses(), 3)
mass_matrix = Numeric.sqrt(Numeric.identity(len(masses))*masses)

skip = 2
num_ts = system._dcd.numframes/skip
num_coor = len(asel)*3

#com = Numeric.zeros((num_ts, 3), Numeric.Float)
#for ts in system._dcd:
#    com[ts.frame-1] = system.F.centerOfMass()

ca_pos = system._dcd.timeseries(asel, skip=skip, format='fac')
#ca_pos = ca_pos-com[:, Numeric.NewAxis]
ca = Numeric.reshape(ca_pos, (num_ts, -1))
ca_avg = Numeric.average(ca)
ca2 = ca - ca_avg[Numeric.NewAxis,:]
ca_cov = Numeric.zeros((num_coor, num_coor), Numeric.Float)
for ts in ca2:
    ca_cov += Numeric.outerproduct(ts, ts)
ca_cov /= num_ts
ca_cov1 = Numeric.matrixmultiply(ca_cov, mass_matrix)
del ca_cov
ca_cov2 = Numeric.matrixmultiply(mass_matrix, ca_cov1)
del ca_cov1

N_av = 6.0221367e23
hplanck_bar = 6.6260755e-34/(2*Numeric.pi)
k =  1.3806580000000001e-23
T = 300 # kelvin
eigenv, eigenvec = LinearAlgebra.eigenvectors(ca_cov2)
real = [e.real/100. for e in eigenv]
f = file('eigenval.dat', 'w')
for i, val in enumerate(real):
    f.write(`i+1` + '\t' + `val` + '\n')
f.close()

eigenval = eigenv*1.6605402e-27*1e-20
omega_i = Numeric.sqrt(k*T/(eigenval))

term = (hplanck_bar*omega_i)/(k*T)
summation_terms = (term/(Numeric.exp(term)-1.))-Numeric.log(1.-Numeric.exp(-term))
S_ho = k*N_av*Numeric.sum(summation_terms)
