# -*- coding: utf-8 -*-
import LinearAlgebra
from MDAnalysis import *
import Numeric

system = AtomGroup.System("data/beta_op_adp-vmd.psf", "data/step1.dcd")
asel = system.select_atoms('segid F and backbone and not type O')

masses = Numeric.repeat(asel.masses, 3)
mass_matrix = Numeric.sqrt(Numeric.identity(len(masses)) * masses)

skip = 1
num_ts = system.dcd.numframes / skip
num_coor = len(asel) * 3

ca_pos = system.dcd.timeseries(asel, skip=skip, format='fac')
collection.addTimeseries(Timeseries.CenterOfMass(asel))
ca_com = Numeric.transpose(system.dcd.correl(collection, skip=skip))

# Recenter on center of mass of selection
ca_pos -= ca_com[:, Numeric.NewAxis]

# Remove rotational degrees of freedom
ref = ca_pos[0]
for coor in ca_pos[1:]:
    rotmatrix = rms_fitting.rms_rotation_matrix(coor, ref, kg)
    coor[:] = Numeric.matrixmultiply(coor, rotmatrix).astype(Numeric.Float32)

ca = Numeric.reshape(ca_pos, (num_ts, -1))
ca_avg = Numeric.average(ca)
ca2 = ca - ca_avg[Numeric.NewAxis, :]
ca_cov = Numeric.zeros((num_coor, num_coor), Numeric.Float)
for ts in ca2:
    ca_cov += Numeric.outerproduct(ts, ts)
ca_cov /= num_ts
ca_cov1 = Numeric.matrixmultiply(ca_cov, mass_matrix)
del ca_cov
ca_cov2 = Numeric.matrixmultiply(mass_matrix, ca_cov1)
del ca_cov1

N_av = 6.0221367e23
hplanck_bar = 6.6260755e-34 / (2 * Numeric.pi)
k = 1.3806580000000001e-23
T = 300  # kelvin
eigenv, eigenvec = LinearAlgebra.eigenvectors(ca_cov2)
real = [e.real / 100. for e in eigenv]
f = open('eigenval.dat', 'w')
for i, val in enumerate(real):
    f.write("%i\t%s\n" % (i + 1, val))
f.close()

eigenval = eigenv * 1.6605402e-27 * 1e-20
omega_i = Numeric.sqrt(k * T / (eigenval))

term = (hplanck_bar * omega_i) / (k * T)
summation_terms = (term / (Numeric.exp(term) - 1.)) - Numeric.log(1. - Numeric.exp(-term))
S_ho = k * N_av * Numeric.sum(summation_terms)
