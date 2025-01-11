# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Leaflet identification --- :mod:`MDAnalysis.analysis.leaflet`
==============================================================

This module implements the *LeafletFinder* algorithm, described in
[Michaud-Agrawal2011]_. It can identify the lipids in a bilayer of
arbitrary shape and topology, including planar and undulating bilayers
under periodic boundary conditions or vesicles.

One can use this information to identify

* the upper and lower leaflet of a *planar membrane* by comparing the
  the :meth:`~MDAnalysis.core.groups.AtomGroup.center_of_geometry` of
  the leaflet groups, or

* the outer and inner leaflet of a *vesicle* by comparing histograms
  of distances from the centre of geometry (or possibly simply the
  :meth:`~MDAnalysis.core.groups.AtomGroup.radius_of_gyration`).

See example scripts in the MDAnalysisCookbook_ on how to use
:class:`LeafletFinder`. The function :func:`optimize_cutoff` implements a
(slow) heuristic method to find the best cut off for the LeafletFinder
algorithm.

.. _MDAnalysisCookbook: https://github.com/MDAnalysis/MDAnalysisCookbook/tree/master/examples


Algorithm
---------

1. build a graph of all phosphate distances < cutoff
2. identify the largest connected subgraphs
3. analyse first and second largest graph, which correspond to the leaflets

For further details see [Michaud-Agrawal2011]_.


Classes and Functions
---------------------

.. autoclass:: LeafletFinder
   :members:

.. autofunction:: optimize_cutoff

"""
import warnings

import numpy as np

from .. import core
from . import distances
from .. import selections
from ..due import due, Doi


# networkx is an optional import
try:
    import networkx as NX
except ImportError:
    HAS_NX = False
else:
    HAS_NX = True


due.cite(
    Doi("10.1002/jcc.21787"),
    description="LeafletFinder algorithm",
    path="MDAnalysis.analysis.leaflet",
    cite_module=True,
)
del Doi


class LeafletFinder(object):
    """Identify atoms in the same leaflet of a lipid bilayer.

    This class implements the *LeafletFinder* algorithm [Michaud-Agrawal2011]_.

    Parameters
    ----------
    universe : Universe
        :class:`~MDAnalysis.core.universe.Universe` object.
    select : AtomGroup or str
        A AtomGroup instance or a
        :meth:`Universe.select_atoms` selection string
        for atoms that define the lipid head groups, e.g.
        universe.atoms.PO4 or "name PO4" or "name P*"
    cutoff : float (optional)
        head group-defining atoms within a distance of `cutoff`
        Angstroms are deemed to be in the same leaflet [15.0]
    pbc : bool (optional)
        take periodic boundary conditions into account [``False``]
    sparse : bool (optional)
        ``None``: use fastest possible routine; ``True``: use slow
        sparse matrix implementation (for large systems); ``False``:
        use fast :func:`~MDAnalysis.lib.distances.distance_array`
        implementation [``None``].

    Raises
    ------
    ImportError
      This class requires the optional package `networkx`. If not found
      an ImportError is raised.

    Example
    -------
    The components of the graph are stored in the list
    :attr:`LeafletFinder.components`; the atoms in each component are numbered
    consecutively, starting at 0. To obtain the atoms in the input structure
    use :meth:`LeafletFinder.groups`::

       u = mda.Universe(PDB)
       L = LeafletFinder(u, 'name P*')
       leaflet0 = L.groups(0)
       leaflet1 = L.groups(1)

    The residues can be accessed through the standard MDAnalysis mechanism::

       leaflet0.residues

    provides a :class:`~MDAnalysis.core.groups.ResidueGroup`
    instance. Similarly, all atoms in the first leaflet are then ::

       leaflet0.residues.atoms


    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`
    .. versionchanged:: 2.0.0
       The universe keyword no longer accepts non-Universe arguments. Please
       create a :class:`~MDAnalysis.core.universe.Universe` first.
    """

    def __init__(self, universe, select, cutoff=15.0, pbc=False, sparse=None):
        # Raise an error if networkx is not installed
        if not HAS_NX:
            errmsg = (
                "The LeafletFinder class requires an installation of "
                "networkx. Please install networkx "
                "https://networkx.org/documentation/stable/install.html"
            )
            raise ImportError(errmsg)

        self.universe = universe
        self.selectionstring = select
        if isinstance(self.selectionstring, core.groups.AtomGroup):
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
                adj = distances.contact_matrix(
                    coord, cutoff=self.cutoff, returntype="numpy", box=box
                )
            except ValueError:  # pragma: no cover
                warnings.warn(
                    "N x N matrix too big, use sparse=True or sparse=None",
                    category=UserWarning,
                    stacklevel=2,
                )
                raise
        elif self.sparse is True:
            # only try sparse
            adj = distances.contact_matrix(
                coord, cutoff=self.cutoff, returntype="sparse", box=box
            )
        else:
            # use distance_array and fall back to sparse matrix
            try:
                # this works for small-ish systems and depends on system memory
                adj = distances.contact_matrix(
                    coord, cutoff=self.cutoff, returntype="numpy", box=box
                )
            except ValueError:  # pragma: no cover
                # but use a sparse matrix method for larger systems for memory reasons
                warnings.warn(
                    "N x N matrix too big - switching to sparse matrix method (works fine, but is currently rather "
                    "slow)",
                    category=UserWarning,
                    stacklevel=2,
                )
                adj = distances.contact_matrix(
                    coord, cutoff=self.cutoff, returntype="sparse", box=box
                )
        return NX.Graph(adj)

    def _get_components(self):
        """Return connected components (as sorted numpy arrays), sorted by size."""
        return [
            np.sort(list(component))
            for component in NX.connected_components(self.graph)
        ]

    def update(self, cutoff=None):
        """Update components, possibly with a different *cutoff*"""
        if cutoff is None:
            cutoff = self.cutoff
        self._init_graph(cutoff)

    def sizes(self):
        """Dict of component index with size of component."""
        return dict(
            (
                (idx, len(component))
                for idx, component in enumerate(self.components)
            )
        )

    def groups(self, component_index=None):
        """Return a :class:`MDAnalysis.core.groups.AtomGroup` for *component_index*.

        If no argument is supplied, then a list of all leaflet groups is returned.

        See Also
        --------
        :meth:`LeafletFinder.group`
        :meth:`LeafletFinder.groups_iter`
        """
        if component_index is None:
            return list(self.groups_iter())
        else:
            return self.group(component_index)

    def group(self, component_index):
        """Return a :class:`MDAnalysis.core.groups.AtomGroup` for *component_index*."""
        # maybe cache this?
        indices = [i for i in self.components[component_index]]
        return self.selection[indices]

    def groups_iter(self):
        """Iterator over all leaflet :meth:`groups`"""
        for component_index in range(len(self.components)):
            yield self.group(component_index)

    def write_selection(self, filename, **kwargs):
        """Write selections for the leaflets to *filename*.

        The format is typically determined by the extension of *filename*
        (e.g. "vmd", "pml", or "ndx" for VMD, PyMol, or Gromacs).

        See :class:`MDAnalysis.selections.base.SelectionWriter` for all
        options.
        """
        sw = selections.get_writer(filename, kwargs.pop("format", None))
        with sw(
            filename,
            mode=kwargs.pop("mode", "w"),
            preamble="leaflets based on select={selectionstring!r} cutoff={cutoff:f}\n".format(
                **vars(self)
            ),
            **kwargs,
        ) as writer:
            for i, ag in enumerate(self.groups_iter()):
                name = "leaflet_{0:d}".format((i + 1))
                writer.write(ag, name=name)

    def __repr__(self):
        return "<LeafletFinder({0!r}, cutoff={1:.1f} A) with {2:d} atoms in {3:d} groups>".format(
            self.selectionstring,
            self.cutoff,
            self.selection.n_atoms,
            len(self.components),
        )


def optimize_cutoff(
    universe,
    select,
    dmin=10.0,
    dmax=20.0,
    step=0.5,
    max_imbalance=0.2,
    **kwargs,
):
    r"""Find cutoff that minimizes number of disconnected groups.

    Applies heuristics to find best groups:

    1. at least two groups (assumes that there are at least 2 leaflets)
    2. reject any solutions for which:

       .. math::

              \frac{|N_0 - N_1|}{|N_0 + N_1|} > \mathrm{max_imbalance}

       with :math:`N_i` being the number of lipids in group
       :math:`i`. This heuristic picks groups with balanced numbers of
       lipids.

    Parameters
    ----------
    universe : Universe
        :class:`MDAnalysis.Universe` instance
    select : AtomGroup or str
        AtomGroup or selection string as used for :class:`LeafletFinder`
    dmin : float (optional)
    dmax : float (optional)
    step : float (optional)
        scan cutoffs from `dmin` to `dmax` at stepsize `step` (in Angstroms)
    max_imbalance : float (optional)
        tuning parameter for the balancing heuristic [0.2]
    kwargs : other keyword arguments
        other arguments for  :class:`LeafletFinder`

    Returns
    -------
    (cutoff, N)
         optimum cutoff and number of groups found


    .. Note:: This function can die in various ways if really no
              appropriate number of groups can be found; it ought  to be
              made more robust.

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`
    """
    kwargs.pop("cutoff", None)  # not used, so we filter it
    _sizes = []
    for cutoff in np.arange(dmin, dmax, step):
        LF = LeafletFinder(universe, select, cutoff=cutoff, **kwargs)
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
