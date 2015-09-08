#!/usr/bin/env python
# Example script, part of MDAnalysis
"""
:Author: Oliver Beckstein
:Year: 2010
:Copyright: GNU Public License v3

MDAnalysis example: Leaflet indentification
===========================================

Simple use of :class:`MDAnalysis.analysis.leaflet.LeafletFinder`:
finds leaflets in a membrane system and write out a selection for
visualization in VMD.

See the source for details and see ``memanal-composition.py`` for a
more sophisticated example.

LeafletFinder Algorithm
-----------------------
  1. build a graph of all phosphate distances < cutoff
  2. identify the largest connected subgraphs
  3. analyse first and second largest graph (assume they are the
     leaflets, anything else are stray micelles or individually
     floating lipids)

One could identify group 1 and 2 as upper and lower leaflet by
comparing the median of the centres of masses; this is left as an
exercise...
"""
from MDAnalysis.analysis.leaflet import LeafletFinder


if __name__ == "__main__":
    import sys

    try:
        PDB, selection = sys.argv[1:3]
    except ValueError:
        print "usage: leaflet.py PDB SELECTION"
        sys.exit(1)
    print "PDB=%(PDB)r selection=%(selection)r" % vars()
    L = LeafletFinder(PDB, selection)
    print "Number of lipids in leaflets: %r" % L.sizes()
    macrovmd = PDB + ".vmd"
    L.write_vmd(macrovmd)
    print "Load macros for vmd from file %r" % macrovmd
