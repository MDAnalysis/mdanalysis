from MDAnalysis import *
from pylab import *

def bb_dihe(kalp, file):
    system = AtomGroup.Universe("data/%s.psf"%file, "data/%s.dcd"%file)

    a = Timeseries.TimeseriesCollection()
    for res in range(2, kalp+2):
        print "Processing residue %d" %res
        phi_sel = system.selectAtoms("atom KALP %d C"%(res-1), "atom KALP %d N"%res, "atom KALP %d CA"%res, "atom KALP %d C" % res)
        psi_sel = system.selectAtoms("atom KALP %d N"%res, "atom KALP %d CA"%res, "atom KALP %d C"%res, "atom KALP %d N" % (res+1))
        a.addTimeseries(Timeseries.Dihedral(phi_sel))
        a.addTimeseries(Timeseries.Dihedral(psi_sel))

    data = system._dcd.correl(a, skip=10)*180./pi
    avg = mean(data, axis=1)
    stdev = std(data, axis=1)
    res = range(2, kalp+2)
    a = errorbar(res, avg[::2], stdev[::2], fmt='ro', label="phi")
    b = errorbar(res, avg[1::2], stdev[1::2], fmt='bo', label="psi")
    legend((a[0], b[0]), ("phi", "psi"), numpoints = 1)
    savefig("results/%kalp%d_bb_dihe"%kalp)
    clf()

if __name__ == "__main__":
    for kalp,file in [(16,"kalp16"), (19,"kalp19"), (23,"kalp23.2")]:
        bb_dihe(kalp,file)
