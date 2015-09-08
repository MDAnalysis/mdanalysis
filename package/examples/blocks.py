#!/usr/bin/env python
# coding=utf-8

"""
MDAnalysis example: Simple blocking analysis
============================================

Calculate the radius of gyration for increasing number of blocks over
a trajectory.

See: H. Flyvbjerg and H. G. Petersen. Error estimates on averages of
correlated data. 91(1):461-466, 1989.

"""

import numpy as np
import MDAnalysis


def blocked(universe, nblocks, analyze):
    size = universe.trajectory.numframes / nblocks
    blocks = []
    for block in xrange(nblocks):
        a = []
        for ts in u.trajectory[block * size:(block + 1) * size]:
            a.append(analyze(universe))
        blocks.append(np.average(a))
    blockaverage = np.average(blocks)
    blockstd = np.std(blocks)

    return nblocks, size, blockaverage, blockstd


def rgyr(universe):
    return universe.select_atoms('protein').radius_of_gyration()


if __name__ == "__main__":
    from MDAnalysis.tests.datafiles import PSF, DCD

    try:
        import matplotlib

        matplotlib.use('agg')  # no interactive plotting, only save figures
        from pylab import errorbar, subplot, xlabel, ylabel, savefig

        have_matplotlib = True
    except ImportError:
        have_matplotlib = False

    u = MDAnalysis.Universe(PSF, DCD)
    results = []
    for nblocks in xrange(2, 10):
        results.append(blocked(u, nblocks, rgyr))
    r = np.array(results)

    if have_matplotlib:
        subplot(211)
        errorbar(r[:, 0], r[:, 2], yerr=r[:, 3])
        xlabel("number of blocks")
        ylabel(r"$\langle R_{\rm{gyr}} \rangle$ ($\AA$)")

        subplot(212)
        errorbar(r[:, 1], r[:, 2], yerr=r[:, 3])
        xlabel("block size")
        ylabel(r"$\langle R_{\rm{gyr}} \rangle$ ($\AA$)")

        savefig("./figures/blocks.pdf")
        savefig("./figures/blocks.png")

        print "Wrote ./figures/blocks.{pdf,png}" % vars()
