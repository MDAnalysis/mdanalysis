#!/usr/bin/python
import sys
import os
import MDAnalysis
import MDAnalysis.analysis.gnm
from MDAnalysis.tests.datafiles import PSF, DCD

u = MDAnalysis.Universe(PSF, DCD)
C = MDAnalysis.analysis.gnm.GNMAnalysis(u, ReportVector="output.txt")

C.run()
output = zip(*C.results)

outputfile = open("eigenvalues.dat", "w")
for item in output[1]:
    print >> outputfile, item
outputfile.close()
