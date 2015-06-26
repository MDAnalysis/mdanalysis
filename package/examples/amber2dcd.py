#!/usr/bin/env python
# coding=utf-8

"""
MDAnalysis example: Convert Amber formatted trajectory to DCD
=============================================================

This example shows how one can use MDAnalysis to convert between
different trajectory formats.

"""

# from MDAnalysis.tests.datafiles import PRMpbc,TRJpbc_bz2
from MDAnalysis.tests.datafiles import PRM, TRJ_bz2
from MDAnalysis import Universe, Writer
from MDAnalysis.lib.util import greedy_splitext

import os.path

topol = PRM  # PRMpbc
intrj = TRJ_bz2  # TRJpbc_bz2
ext = '.dcd'  # output format determined by extension

root, oldext = greedy_splitext(os.path.basename(intrj))
outtrj = root + ext
outpdb = root + '.pdb'

u = Universe(topol, intrj)

# create a writer instance for the output trajectory
w = Writer(outtrj, u.trajectory.numatoms)

# loop through the trajectory and write a frame for every step
for ts in u.trajectory:
    w.write(ts)
    print "Converted frame %d" % ts.frame
w.close_trajectory()
print "Converted %r --> %r" % (intrj, outtrj)

# make a pdb file as a simple 'topology'
u.trajectory.rewind()
u.atoms.write(outpdb)
print "Created %r to be used with the trajectory" % outpdb
