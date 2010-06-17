#!/usr/bin/env python
# Example script, part of MDAnalysis
"""
:Author: Oliver Beckstein
:Year: 2010
:Copyright: GNU Public License v3

Leaflet indentification
=======================

Algorithm:
  1. build a graph of all phosphate distances < cutoff
  2. identify the largest connected subgraphs
  3. analyse first and second largest graph and identify as upper and lower
     leaflet by comparing the median of the centres of masses
"""
from __future__ import with_statement

import numpy
import MDAnalysis
import networkx as NX

class LeafletFinder(object):
    """Identify atoms in the same leaflet of a lipid bilayer.

    Algorithm:

      1. build a graph of all phosphate distances < cutoff
      2. identify the largest connected subgraphs
      3. analyse first and second largest graph and identify as upper and lower
         leaflet by comparing the median of the centres of masses

    The components of the graph are stored in the list
    :attr:`LeafletFinder.components`; the atoms in each component are numbered
    consecutively, starting at 0. To obtain the atoms in the input structure
    use :meth:`LeafletFinder.atoms`::

       L = LeafletFinder(PDB, 'name P*')
       leaflet0 = L.atoms(0)
       leaflet1 = L.atoms(1)
    """

    def __init__(self, universe, selectionstring, cutoff=15.0, pbc=False):
        """Initialize from a *universe* or pdb file.

        :Arguments:
             *universe*
                 :class:`MDAnalysis.Universe` or a PDB file name.
             *selectionstring*
                 :meth:`MDAnalysis.Universe.selectAtoms` selection
                 for atoms that define the lipid head groups, e.g.
                 "name PO4" or "name P*"
        :Keywords:
             *cutoff*
                 head group-defining atoms within a distance of *cutoff*
                 Angstroem are deemed to be in the same leaflet [15.0]
             *pbc*
                 take periodic boundary conditions into account (only works
                 for orthorhombic boxes) [``False``]
        """
        if type(universe) is str:
            universe = MDAnalysis.Universe(universe)
        self.universe = universe
        self.selectionstring = selectionstring
        self.selection = universe.selectAtoms(selectionstring)
        self.pbc = pbc
        self._init_graph(cutoff)

    def _init_graph(self, cutoff):
        self.cutoff = cutoff
        self.graph = self._get_graph()
        self.components = self._get_components()

    # The last two lines of code in _get_graph() and the single line in
    # _get_components() is all that's needed to make the leaflet detection
    # work.

    def _get_graph(self):
        """Build graph from adjacency matrix at the given cutoff."""
        # could use self_distance_array to speed up but then need to deal with the sparse indexing
        if self.pbc:
            box = self.universe.trajectory.ts.dimensions
        else:
            box = None
        coord = self.selection.coordinates()        
        adj =  (MDAnalysis.distances.distance_array(coord,coord,box=box) < self.cutoff)
        return NX.Graph(adj)
        
    def _get_components(self):
        """Return connected components (as sorted numpy arrays), sorted by size."""
        return [numpy.sort(component) for component in NX.connected_components(self.graph)]

    def update(self, cutoff=None):
        """Update components, possibly with a different *cutoff*"""
        if cutoff is None:
            cutoff = self.cutoff
        self._init_graph(cutoff)

    def sizes(self):
        """Dict of component index with size of component."""
        return dict(((idx,len(component)) for idx,component in enumerate(self.components)))

    def atoms(self, component_index):
        """Return a :class:`MDAnalysis.AtomGroup.AtomGroup` for *component_index*."""
        # maybe cache this?
        return MDAnalysis.AtomGroup.AtomGroup(
            [self.selection[i] for i in self.components[component_index]])

    def write_vmd(self, filename, numcols=8):
        """Write VMD atomselect macros to *filename*."""
        # should be in MDAnalysis as part of AtomGroup.write()
        with open(filename, 'w') as vmd:
            vmd.write("# leaflets based on selection=%(selectionstring)r cutoff=%(cutoff)f\n" % vars(self))
            for i in xrange(len(self.components)):                
                name = "leaflet_%d" % (i+1)
                vmd.write("atomselect macro %(name)s {index " % vars())
                for iatom, atom in enumerate(self.atoms(i)):
                    index = atom.number  # VMD index is 0-based (as is MDAnalysis?)
                    vmd.write(" %(index)d" % vars())
                    if (iatom+1) % numcols == 0:
                        vmd.write("\\\n\t")
                vmd.write(" }\n")

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
    macrovmd = PDB+".vmd"
    L.write_vmd(macrovmd)
    print "Load macros for vmd from file %r" % macrovmd
    


