#!/usr/bin/env python
# coding=utf-8
"""
:Author: Oliver Beckstein
:Year: 2010
:Copyright: GNU Public License v3

MDAnalysis example: Lipid bilayer position and thickness
========================================================

This example shows how to use the
:class:`MDAnalysis.analysis.leaflet.LeafletFinder` to determine the
centre of a bilayer and its approximate thickness of a approximately planar bilayer.

.. Note::

   This is a very promitive script. Depending on what you want to know
   and what your system looks like you will have to work significantly
   harder to compute appropriate numbers. This approach is not
   generally recommended for quantitative work!

Algorithm
---------

We compute the z coordinates of the centroids (center of geometry) of
the upper and lower leaflet, named *c1z* and *c2z*. Leaflets are
identified through the phosphate headgroup atoms.  The z-coordinate of
the **center of the bilayer** is approximately

  z_mem = (c1z + c2z)/2

The approximate thickness is

  t_mem = abs(c1z - c2z)

"""

usage = """%%prog [options] structure-file

Calculate the approximate thickness and z-coordinate of the center of
a planar lipid bilayer. The bilayer normal is assumed to be roughly
parallel to the z-axis.

The default selection string (--headgroup-selection) is suitable for
coarse-grained bilayers including phospholipids and
cholesterol. Adjust it so that it selects one atom from each head
group.
"""

import numpy as np
import MDAnalysis
from MDAnalysis.analysis.leaflet import LeafletFinder, optimize_cutoff


def get_membrane_parameters(universe, leafletfinder):
    L = leafletfinder

    L0 = L.group(0)
    L1 = L.group(1)

    parameters = {}
    parameters['thickness'] = np.abs((L1.centroid() - L0.centroid())[2])  # z coordinate
    parameters['zmem'] = 0.5 * (L1.centroid() + L0.centroid())[2]  # z coordinate
    return parameters


if __name__ == "__main__":
    import errno
    import os.path
    from optparse import OptionParser

    parser = OptionParser(usage=usage % vars())
    parser.add_option("-t", "--topology", dest="topology",
                      help="optional topology file e.g. PSF; useful when using "
                           "a PDB file as structure, which only provides 3-letter resnames")
    parser.add_option("--pbc", dest="pbc", action="store_true",
                      help="take periodic boundaries into account? "
                           "(Only works for orthorhombic boxes) [%default]")
    parser.add_option("-d", "--cutoff", dest="cutoff", type="float", metavar="DIST",
                      help="look for neighbouring lipids within DIST Angstroem [%default]")
    parser.add_option("-S", "--headgroup-selection", dest="selection", metavar="SELECTION",
                      help="MDAnalysis selection string that selects one atom from "
                           "each headgroup of every lipid ['%default']")
    parser.add_option("--optimize", dest="optimize", action="store_true",
                      help="find cutoff automatically that minimizes the total number "
                           "of unconnected groups but is larger than a single group. "
                           "Can take a while... Ignores --cutoff. [%default]")
    parser.add_option("--max-imbalance", dest="max_imbalance", type="float",
                      metavar="Q",
                      help="When optimizing, discard any solutions for which "
                           "|N0 - N1|/|N0 + N1| > Q (Ni = number of lipids in group i) "
                           "This heuristic picks groups with balanced numbers of lipids. "
                           "[%default]")

    parser.set_defaults(pbc=False, cutoff=15.0, optimize=False,
                        max_imbalance=0.2,
                        # PO-lipids (and CHOL for CG) but not CHARMM K+
                        selection="(name P* or name ROH) and not name POT")
    # TODO: selection should be set based on the identities of lipids.
    # combined with fingerprinting in the future to do this automagically;
    #       Hard coded for now.

    options, args = parser.parse_args()

    try:
        structure = args[0]
    except IndexError:
        raise ValueError("Need structure file (pdb, gro, ...) as input")

    if not os.path.exists(structure):
        raise IOError(errno.ENOENT, "Structure file not found", structure)

    if not options.topology:
        options.topology = structure

    u = MDAnalysis.Universe(options.topology, structure)

    if options.optimize:
        print "# Finding best cutoff (--optimize)..."
        try:
            cutoff, N = optimize_cutoff(u, options.selection, pbc=options.pbc,
                                        max_imbalance=options.max_imbalance)
        except:
            raise RuntimeError("Failed cutoff optimization, try without --optimize")
        print "# Optimized cutoff=%(cutoff).1f A, finding %(N)d disconnected groups" % vars()
    else:
        cutoff = options.cutoff
        print "# Using fixed cutoff=%(cutoff).1f A" % vars()

    LF = LeafletFinder(u, options.selection, cutoff=cutoff, pbc=options.pbc)
    p = get_membrane_parameters(u, LF)

    # show results
    print "#" + 60 * "="
    print "thickness tmem = %(thickness).2f A" % p
    print "center    zmem = %(zmem).2f A" % p
    print "#" + 60 * "="
