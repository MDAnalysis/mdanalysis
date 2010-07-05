# MDAnalysis example: backbone dihedrals
"""Calculating backbone dihedrals of a protein, using timeseries functionality

Note that you have to modify this script. For instance, change the
selections; it is unlikely that your protein contains a molecule named
'KALP'...
"""


from MDAnalysis import *
from pylab import *
import numpy

dcdfile = ''
psffile = ''
numresidues = 20

universe = Universe(psffile, dcdfile)
protein = universe.selectAtoms("protein")
phi = []
psi = []

for res in range(2, numresidues-1):
    print "Processing residue %d" %res
    phi_sel = universe.selectAtoms("atom KALP %d C"%(res-1), "atom KALP %d N"%res, "atom KALP %d CA"%res, "atom KALP %d C" % res)
    psi_sel = universe.selectAtoms("atom KALP %d N"%res, "atom KALP %d CA"%res, "atom KALP %d C"%res, "atom KALP %d N" % (res+1))
    a = collection.addTimeseries(Timeseries.Dihedral(phi_sel))
    b = collection.addTimeseries(Timeseries.Dihedral(psi_sel))

    phi_data = universe.dcd.correl(a, skip=10)*180./pi
    psi_data = universe.dcd.correl(b, skip=10)*180./pi

    phi_avg = mean(data_phi, axis=1)
    phi_stdev = std(data_phi, axis=1)
    psi_avg = mean(data_psi, axis=1)
    psi_stdev = std(data_psi, axis=1)

    phi.append([phi_avg,phi_stdev])
    psi.append([psi_avg,psi_stdev])

phi = numpy.array(phi)
psi = numpy.array(psi)
res = range(2, kalp+2)
a = errorbar(res, phi[:,0], yerr=phi[:,1], fmt='ro', label="phi")
b = errorbar(res, psi[:,0], yerr=psi[:,1], fmt='bo', label="psi")
legend((a[0], b[0]), ("phi", "psi"), numpoints = 1)
savefig("backbone_dihedrals.png"%kalp)
