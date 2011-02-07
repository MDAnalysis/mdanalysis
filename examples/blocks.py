#!/usr/bin/env python

"""
MDAnalysis example: Simple blocking analysis
============================================

Calculate the radius of gyration for increasing number of blocks over
a trajectory.

See: H. Flyvbjerg and H. G. Petersen. Error estimates on averages of
correlated data. 91(1):461-466, 1989.

"""

import numpy
import MDAnalysis

def blocked(universe, nblocks, analyze):
    size = universe.trajectory.numframes/nblocks
    blocks = []
    for block in xrange(nblocks):
        a = []
        for ts in u.trajectory[block*size:(block+1)*size]:
            a.append(analyze(universe))
        blocks.append(numpy.average(a))
    blockaverage = numpy.average(blocks)
    blockstd = numpy.std(blocks)

    return nblocks, size, blockaverage, blockstd 

def rgyr(universe):
    return universe.selectAtoms('protein').radiusOfGyration()


if __name__ == "__main__":
    from  MDAnalysis.tests.datafiles import PSF,DCD
    from pylab import errorbar, subplot, xlabel, ylabel
    u = MDAnalysis.Universe(PSF,DCD)
    results = []
    for nblocks in xrange(2,10):
        results.append(blocked(u, nblocks, rgyr))
    r = numpy.array(results)
    
    subplot(211)
    errorbar(r[:,0], r[:,2], yerr=r[:,3])
    xlabel("number of blocks")

    subplot(212)
    errorbar(r[:,1], r[:,2], yerr=r[:,3])
    xlabel("block size")



