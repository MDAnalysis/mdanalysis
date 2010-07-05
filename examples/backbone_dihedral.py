# MDAnalysis example: backbone dihedrals
"""Calculating backbone dihedrals of a protein, using timeseries functionality

Note that you have to modify this script. For instance, change the
selections; it is unlikely that your protein contains a molecule named
'KALP'...
"""


from MDAnalysis import *
from pylab import *
from math import pi

dcdfile = ''
psffile = ''
numresidues = 20

universe = Universe(psffile, dcdfile)
protein = universe.selectAtoms("protein")

phi = []
psi = []
for res in range(2, numresidues-1):
    print "Processing residue %d" %res
    #  selection of the atoms involved for the phi for resid '%d' %res
    phi_sel = universe.selectAtoms("atom %d C"%(res-1), "atom %d N"%res, "atom %d CA"%res, "atom %d C" % res)

    #  selection of the atoms involved for the psi for resid '%d' %res
    psi_sel = universe.selectAtoms("atom %d N"%res, "atom %d CA"%res, "atom %d C"%res, "atom %d N" % (res+1))
    
    # definition collecting the timeseries of a dihedral
    a = collection.addTimeseries(Timeseries.Dihedral(phi_sel))
    b = collection.addTimeseries(Timeseries.Dihedral(psi_sel))

    # collection of the timeseries data for every 10 steps in the traj
    data_phi = universe.dcd.correl(a, skip=10)*180./pi
    data_psi = universe.dcd.correl(b, skip=10)*180./pi

    # finding the avg and stdev for each residue
    avg_phi = mean(data_phi)
    stdev_phi = std(data_phi)
    avg_psi = mean(data_psi)
    stdev_psi = std(data_psi)
    phi.append([res,avg_phi,stdev_phi])
    psi.append([res,avg_psi,stdev_psi])

# making an array for phi and psi data
phi = numpy.array(phi)
psi = numpy.array(psi)

# plotting and saving the dihe for each resid
res = range(2, numresidues-1)
a = errorbar(phi[:,0], phi[:,1], phi[:,2], fmt='ro', label="phi")
b = errorbar(psi[:,0], psi[:,1], psi[:,2], fmt='bs', label="psi")
legend((a[0], b[0]), ("phi", "psi"))
savefig("backbone_dihedrals.png")
