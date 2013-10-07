#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# Oliver Beckstein, 2013, released into the Public Domain

"""
Alternative combined rotamer library as a TRR file
--------------------------------------------------

Combine rotamer structures with probability into a TRR file where the
lambda holds the probability.

To get structure and probability::

   u = MDAnalysis.Universe("rotamer1_R1A_298K.pdb", trr)
   for ts in u.trajectory:
       probability = ts.lmbda
       rotamer = u.atoms.positions  # or ts._pos

"""

import numpy as np
import MDAnalysis
from itertools import izip

def build_trr(trr="./rotamer1_R1A_298K.trr"):
    L = MDAnalysis.Universe("./rotamer1_R1A_298K.pdb", "./rotamer1_R1A_298K.dcd")
    pop = np.loadtxt("R1A_298K_populations.dat")
    with MDAnalysis.Writer(trr, L.atoms.numberOfAtoms()) as W:
        for ts, p in izip(L.trajectory, pop):
            ts.lmbda = p  # add lambda for TRR
            W.write(ts)
    return trr

def test(trr="rotamer1_R1A_298K.trr"):
    u = MDAnalysis.Universe("rotamer1_R1A_298K.pdb", trr)
    pp = np.array([ts.lmbda for ts in u.trajectory])
    ref = np.loadtxt("R1A_298K_populations.dat")
    assert np.abs(pp - ref).max() < 1e-10

if __name__ == "__main__":
    import argparse
    # http://stackoverflow.com/questions/12151306/argparse-way-to-include-default-values-in-help
    # http://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments
    parser = argparse.ArgumentParser(description='Build a TRR-formated rotamer library from the original DCD and dat files.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('filename', nargs="?", default="./rotamer1_R1A_298K.trr",
                        help="Write library to TRR filename")
    args = parser.parse_args()
    trr = build_trr(trr=args.trr)
    print("Created TRR-formatted rotamer library {0}".format(trr))
    print("Testing...")
    try:
        test(trr)
    except AssertionError:
        print("New library contains errors. Do not use and investigate!")
        sys.exit(1)
