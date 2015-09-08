#!/usr/bin/env python
# Example script, part of MDAnalysis
"""
:Author: Oliver Beckstein
:Year: 2010
:Copyright: GNU Public License v3

MDAnalysis example: Lipid bilayer composition
=============================================

Find the leaflets of a membrane system and print a breakdown of the
counts of the different residue names in each leaflet.
"""

usage = """%%prog [options] structure-file

Calculate the lipid composition in the two leaflets of the membrane
stored in *structure-file*.

The default selection string (--headgroup-selection) is suitable for
coarse-grained bilayers including phospholipids and
cholesterol. Adjust it so that it selects one atom from each head
group.
"""

import MDAnalysis
from MDAnalysis.analysis.leaflet import LeafletFinder, optimize_cutoff

import numpy as np

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
        raise IOError(errno.ENOENT, "PQR file not found", structure)

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

    print "# group sizes = %r " % LF.sizes()

    # two leaflets
    def print_line(symbol="-"):
        print "#" + (12 + 5) * symbol

    print_line("=")
    print "#%2s  %5s  %6s" % ("ll", "resn", "count")
    print_line("=")

    for groupindex in xrange(len(LF.components)):
        resnames = [a.resname for a in LF.groups(groupindex)]
        # there CERTAINLY is a better way to count occurrences than this...
        keys = np.unique(resnames)
        for k in keys:
            count = resnames.count(k)
            print " %2d  %5s  %6d" % (groupindex, k, count)
        total = LF.sizes()[groupindex]
        if total > 1:
            print_line()
            print "#%2d  %5s  %6d" % (groupindex, '', total)
        print
