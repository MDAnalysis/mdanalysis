# lipid membrane leaflet analysis
# Part of MDAnalysis http://mdanalysis.googlecode.com
# Copyright (c) 2006-2010 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# Released under the GNU Public Licence, v2

"""
:Author: Oliver Beckstein
:Year: 2010
:Copyright: GNU Public License v3

:mod:`MDAnalysis.analysis.leaflet` --- Leaflet indentification
==============================================================

Algorithm:
  1. build a graph of all phosphate distances < cutoff
  2. identify the largest connected subgraphs
  3. analyse first and second largest graph, which correspond to the leaflets

You could identify the upper and lower leaflet of a planar membrane by
comparing the median of the centres of masses, or for a vesicle by
comparing distances from the centre of geometry, but neither is
implemented at the moment.


.. autoclass:: LeafletFinder
   :members:

"""
from __future__ import with_statement

import numpy
import MDAnalysis
import networkx as NX

class LeafletFinder(object):
    """Identify atoms in the same leaflet of a lipid bilayer.

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

    def atoms_iter(self):
        for component_index in xrange(len(self.components)):
            yield self.atoms(component_index)

    def write_selection(self, filename, **kwargs):
        """Write selections for the leaflets to *filename*.

        The format is typically determined by the extension of *filename*
        (e.g. "vmd", "pml", or "ndx" for VMD, PyMol, or Gromacs).

        See :class:`MDAnalysis.selections.base.SelectionWriter` for all
        options.
        """
        import MDAnalysis.selections        
        SelectionWriter = MDAnalysis.selections.get_writer(filename, kwargs.pop('format',None))
        writer = SelectionWriter(
            filename, mode=kwargs.pop('mode', 'wa'), 
            preamble="leaflets based on selection=%(selectionstring)r cutoff=%(cutoff)f\n" % vars(self), 
            **kwargs)
        for i, ag in enumerate(self.atoms_iter()):
            name = "leaflet_%d" % (i+1)
            writer.write(ag, name=name)


