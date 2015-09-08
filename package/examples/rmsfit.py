#!/usr/bin/env python
# Example script, part of MDAnalysis
"""
RMS-fit a trajectory to a reference structure
=============================================

Simple example implementation that shows how to fit by rms and access
the translated and rotated coordinates.

"""

import numpy as np
import MDAnalysis
from MDAnalysis.analysis.align import rms_fit_trj

if __name__ == '__main__':
    from MDAnalysis import Universe
    from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small

    try:
        import matplotlib

        matplotlib.use('agg')
        from pylab import plot, xlabel, ylabel, savefig

        have_matplotlib = True
    except:
        have_matplotlib = False

    ref = Universe(PSF, PDB_small)  # reference structure 4AKE
    trj = Universe(PSF, DCD)  # trajectory of change 1AKE->4AKE

    rms_fit_trj(trj, ref, filename="output/rmsfit.dcd", select="not name H*",
                mass_weighted=True, rmsdfile="output/rmsfit_rmsd.dat")
    print "Fitted trajectory: output/rmsfit.dcd"

    rmsd = np.loadtxt("output/rmsfit_rmsd.dat")
    print "RMSD: output/rmsfit_rmsd.dat"

    if have_matplotlib:
        plot(rmsd, 'k-', linewidth=2)
        xlabel("frame number")
        ylabel(r"mass weighted heavy atom RMSD ($\AA$)")
        savefig("figures/rmsdfit_rmsd.pdf")
        savefig("figures/rmsdfit_rmsd.png")
        print "RMSD graph: figures/rmsdfit_rmsd.{pdf,png}"
