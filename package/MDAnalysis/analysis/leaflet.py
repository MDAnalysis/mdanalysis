# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Leaflet identification --- :mod:`MDAnalysis.analysis.leaflet`
==============================================================

:Author: Oliver Beckstein
:Year: 2010
:Copyright: GNU Public License v3

Algorithm:
  1. build a graph of all phosphate distances < cutoff
  2. identify the largest connected subgraphs
  3. analyse first and second largest graph, which correspond to the leaflets

One can use this information to identify

* the upper and lower leaflet of a *planar membrane* by comparing the
  the :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.center_of_geometry` of
  the leaflet groups, or

* the outer and inner leaflet of a *vesicle* by comparing histograms
  of distances from the centre of geometry (or possibly simply the
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.radius_of_gyration`).

See example scripts in the ``examples/`` directory on how to use
:class:`LeafletFinder`. The function :func:`optimize_cutoff` implements a
(slow) heuristic method to find the best cut off for the LeafletFinder
algorithm.

.. autoclass:: LeafletFinder
   :members:

.. autofunction:: optimize_cutoff

"""

import numpy as np
import MDAnalysis
import networkx as NX
import distances
import warnings


class LeafletFinder(object):
    """Identify atoms in the same leaflet of a lipid bilayer.

    The components of the graph are stored in the list
    :attr:`LeafletFinder.components`; the atoms in each component are numbered
    consecutively, starting at 0. To obtain the atoms in the input structure
    use :meth:`LeafletFinder.groups`::

       L = LeafletFinder(PDB, 'name P*')
       leaflet0 = L.groups(0)
       leaflet1 = L.groups(1)

    The residues can be accessed through the standard MDAnalysis mechanism::

       leaflet0.residues

    provides a :class:`~MDAnalysis.core.AtomGroup.ResidueGroup`
    instance. Similarly, all atoms in the first leaflet are then ::

       leaflet0.residues.atoms
    """

    def __init__(self, universe, selectionstring, cutoff=15.0, pbc=False, sparse=None):
        """Initialize from a *universe* or pdb file.

        :Arguments:
             *universe*
                 :class:`MDAnalysis.Universe` or a PDB file name.
             *selection*
                 :class:`MDAnalysis.core.AtomGroup.AtomGroup` or a
                 :meth:`MDAnalysis.Universe.select_atoms` selection string
                 for atoms that define the lipid head groups, e.g.
                 universe.atoms.PO4 or "name PO4" or "name P*"
        :Keywords:
             *cutoff*
                 head group-defining atoms within a distance of *cutoff*
                 Angstroms are deemed to be in the same leaflet [15.0]
             *pbc*
                 take periodic boundary conditions into account (only works
                 for orthorhombic boxes) [``False``]
             *sparse*
                 ``None``: use fastest possible routine; ``True``: use slow
                 sparse matrix implementation (for large systems); ``False``:
                 use fast :func:`~MDAnalysis.analysis.distances.distance_array`
                 implementation [``None``].
        """
        universe = MDAnalysis.as_Universe(universe)
        self.universe = universe
        self.selectionstring = selectionstring
        if type(self.selectionstring) == MDAnalysis.core.AtomGroup.AtomGroup:
            self.selection = self.selectionstring
        else:
            self.selection = universe.select_atoms(self.selectionstring)
        self.pbc = pbc
        self.sparse = sparse
        self._init_graph(cutoff)

    def _init_graph(self, cutoff):
        self.cutoff = cutoff
        self.graph = self._get_graph()
        self.components = self._get_components()

    # The last two calls in _get_graph() and the single line in
    # _get_components() are all that are needed to make the leaflet
    # detection work.

    def _get_graph(self):
        """Build graph from adjacency matrix at the given cutoff.
        Automatically select between high and low memory usage versions of
        contact_matrix."""
        # could use self_distance_array to speed up but then need to deal with the sparse indexing
        if self.pbc:
            box = self.universe.trajectory.ts.dimensions
        else:
            box = None
        coord = self.selection.positions
        if self.sparse is False:
            # only try distance array
            try:
                adj = distances.contact_matrix(coord, cutoff=self.cutoff, returntype="numpy", box=box)
            except ValueError:
                warnings.warn('N x N matrix too big, use sparse=True or sparse=None', category=UserWarning,
                              stacklevel=2)
                raise
        elif self.sparse is True:
            # only try sparse
            adj = distances.contact_matrix(coord, cutoff=self.cutoff, returntype="sparse", box=box)
        else:
            # use distance_array and fall back to sparse matrix
            try:
                # this works for small-ish systems and depends on system memory
                adj = distances.contact_matrix(coord, cutoff=self.cutoff, returntype="numpy", box=box)
            except ValueError:
                # but use a sparse matrix method for larger systems for memory reasons
                warnings.warn(
                    'N x N matrix too big - switching to sparse matrix method (works fine, but is currently rather '
                    'slow)',
                    category=UserWarning, stacklevel=2)
                adj = distances.contact_matrix(coord, cutoff=self.cutoff, returntype="sparse", box=box)
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
        return dict(((idx, len(component)) for idx, component in enumerate(self.components)))

    def groups(self, component_index=None):
        """Return a :class:`MDAnalysis.core.AtomGroup.AtomGroup` for *component_index*.

        If no argument is supplied, then a list of all leaflet groups is returned.

        .. SeeAlso:: :meth:`LeafletFinder.group` and :meth:`LeafletFinder.groups_iter`
        """
        if component_index is None:
            return list(self.groups_iter())
        else:
            return self.group(component_index)

    def group(self, component_index):
        """Return a :class:`MDAnalysis.core.AtomGroup.AtomGroup` for *component_index*."""
        # maybe cache this?
        return MDAnalysis.core.AtomGroup.AtomGroup(
            [self.selection[i] for i in self.components[component_index]])

    def groups_iter(self):
        """Iterator over all leaflet :meth:`groups`"""
        for component_index in xrange(len(self.components)):
            yield self.group(component_index)

    def write_selection(self, filename, **kwargs):
        """Write selections for the leaflets to *filename*.

        The format is typically determined by the extension of *filename*
        (e.g. "vmd", "pml", or "ndx" for VMD, PyMol, or Gromacs).

        See :class:`MDAnalysis.selections.base.SelectionWriter` for all
        options.
        """
        import MDAnalysis.selections

        SelectionWriter = MDAnalysis.selections.get_writer(filename, kwargs.pop('format', None))
        writer = SelectionWriter(
            filename, mode=kwargs.pop('mode', 'wa'),
            preamble="leaflets based on selection=%(selectionstring)r cutoff=%(cutoff)f\n" % vars(self),
            **kwargs)
        for i, ag in enumerate(self.groups_iter()):
            name = "leaflet_%d" % (i + 1)
            writer.write(ag, name=name)

    def __repr__(self):
        return "<LeafletFinder(%r, cutoff=%.1f A) with %d atoms in %d groups>" % \
               (self.selectionstring, self.cutoff, self.selection.n_atoms,
               len(self.components))


def optimize_cutoff(universe, selection, dmin=10.0, dmax=20.0, step=0.5,
                    max_imbalance=0.2, **kwargs):
    """Find cutoff that minimizes number of disconnected groups.

    Applies heuristics to find best groups:

    1. at least two groups (assumes that there are at least 2 leaflets)
    2. reject any solutions for which:

              `|N0 - N1|/|N0 + N1| > *max_imbalance*`

       Ni = number of lipids in group i. This heuristic picks groups with
       balanced numbers of lipids.

    :Arguments:
      *universe*
          :class:`MDAnalysis.Universe` instance
      *selection*
          selection string as used for :class:`LeafletFinder`
      *dmin*, *dmax*, *step*
          scan cutoffs from *dmin* to *dmax* at stepsize*step (in Angstroms)
      *max_imbalance*
          tuning parameter for the balancing heuristic (2) [0.2]
      *kwargs*
          other arguments for  :class:`LeafletFinder`

    :Returns: ``(cutoff,N)`` optimum cutoff and number of groups found
    :Raises: can die in various ways if really no appropriate number of groups
             can be found; needs to be made more robust
    """
    kwargs.pop('cutoff', None)  # not used, so we filter it
    _sizes = []
    for cutoff in np.arange(dmin, dmax, step):
        LF = LeafletFinder(universe, selection, cutoff=cutoff, **kwargs)
        # heuristic:
        #  1) N > 1
        #  2) no imbalance between large groups:
        sizes = LF.sizes()
        if len(sizes) < 2:
            continue
        n0 = float(sizes[0])  # sizes of two biggest groups ...
        n1 = float(sizes[1])  # ... assumed to be the leaflets
        imbalance = np.abs(n0 - n1) / (n0 + n1)
        # print "sizes: %(sizes)r; imbalance=%(imbalance)f" % vars()
        if imbalance > max_imbalance:
            continue
        _sizes.append((cutoff, len(LF.sizes())))
    results = np.rec.fromrecords(_sizes, names="cutoff,N")
    del _sizes
    results.sort(order=["N", "cutoff"])  # sort ascending by N, then cutoff
    return results[0]  # (cutoff,N) with N>1 and shortest cutoff
