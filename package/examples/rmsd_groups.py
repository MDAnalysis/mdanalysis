#!/usr/bin/env python
# Example script, part of MDAnalysis
"""
RMSD calculation for a group of atoms
=====================================

See docs for MDAnalysis.analysis.rms

"""

import MDAnalysis
from MDAnalysis.tests.datafiles import PSF, DCD, CRD

import MDAnalysis.analysis.rms

u = MDAnalysis.Universe(PSF, DCD)
ref1 = MDAnalysis.Universe(PSF, DCD)  # reference closed AdK (1AKE) (with the default ref_frame=0)

R1 = MDAnalysis.analysis.rms.RMSD(u, ref1,
                                  select="backbone",  # superimpose on whole backbone of the whole protein
                                  groupselections=[
                                      "backbone and (resid 1-29 or resid 60-121 or resid 160-214)",  # CORE
                                      "backbone and resid 122-159",  # LID
                                      "backbone and resid 30-59"],  # NMP
                                  filename="rmsd_all_CORE_LID_NMP_ref1AKE.dat")
R1.run()
R1.save()

ref2 = MDAnalysis.Universe(PSF, CRD)  # reference open AdK (4AKE)
R2 = MDAnalysis.analysis.rms.RMSD(u, ref2,
                                  select="backbone",  # superimpose on whole backbone of the whole protein
                                  groupselections=[
                                      "backbone and (resid 1-29 or resid 60-121 or resid 160-214)",  # CORE
                                      "backbone and resid 122-159",  # LID
                                      "backbone and resid 30-59"],  # NMP
                                  filename="rmsd_all_CORE_LID_NMP_ref4AKE.dat")
R2.run()
R2.save()

# -----------
# plotting
#-----------

import matplotlib.pyplot as plt

rmsd1 = R1.rmsd.T  # transpose makes it easier for plotting
rmsd2 = R2.rmsd.T  # transpose makes it easier for plotting
time = rmsd1[1]

fig = plt.figure(figsize=(8, 4))

ax1 = fig.add_subplot(121)
ax1.plot(time, rmsd1[2], 'k-', label="all")
ax1.plot(time, rmsd1[3], 'k--', label="CORE")
ax1.plot(time, rmsd1[4], 'r--', label="LID")
ax1.plot(time, rmsd1[5], 'b--', label="NMP")
#ax1.legend(loc="best")
ax1.set_xlabel("time (ps)")
ax1.set_ylabel(r"RMSD ($\AA$)")
ax1.set_title("reference: 1AKE (closed)")

ax2 = fig.add_subplot(122)
ax2.plot(time, rmsd2[2], 'k-', label="all")
ax2.plot(time, rmsd2[3], 'k--', label="CORE")
ax2.plot(time, rmsd2[4], 'r--', label="LID")
ax2.plot(time, rmsd2[5], 'b--', label="NMP")
ax2.legend(loc="best")
ax2.set_xlabel("time (ps)")
#ax2.yaxis.set_visible(False)
ax2.set_ylim(ax1.get_ylim())
ax2.set_title("reference: 4AKE (open)")

fig.savefig("rmsd_all_CORE_LID_NMP.pdf")
fig.savefig("rmsd_all_CORE_LID_NMP.svg")
fig.savefig("rmsd_all_CORE_LID_NMP.png")
