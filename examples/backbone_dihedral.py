# MDAnalysis example: backbone dihedrals
"""Calculating backbone dihedrals of a protein, using timeseries functionality

Note that you have to modify this script. For instance, change the
selections; it is unlikely that your protein contains a molecule named
'KALP'...
"""


from MDAnalysis import *
from pylab import *

dcdfile = ''
psffile = ''
numresidues = 20

universe = Universe(psffile, dcdfile)
protein = universe.selectAtoms("protein")

for res in range(2, numresidues-1):
    print "Processing residue %d" %res
    phi_sel = universe.selectAtoms("atom KALP %d C"%(res-1), "atom KALP %d N"%res, "atom KALP %d CA"%res, "atom KALP %d C" % res)
    psi_sel = universe.selectAtoms("atom KALP %d N"%res, "atom KALP %d CA"%res, "atom KALP %d C"%res, "atom KALP %d N" % (res+1))
    collection.addTimeseries(Timeseries.Dihedral(phi_sel))
    collection.addTimeseries(Timeseries.Dihedral(psi_sel))

data = universe.dcd.correl(a, skip=10)*180./pi
avg = mean(data, axis=1)
stdev = std(data, axis=1)
res = range(2, kalp+2)
a = errorbar(res, avg[::2], stdev[::2], fmt='ro', label="phi")
b = errorbar(res, avg[1::2], stdev[1::2], fmt='bo', label="psi")
legend((a[0], b[0]), ("phi", "psi"), numpoints = 1)
savefig("backbone_dihedrals.png"%kalp)
