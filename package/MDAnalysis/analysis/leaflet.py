# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
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
import networkx as NX

from .. import core, selections, _TOPOLOGY_ATTRNAMES, _TOPOLOGY_ATTRS
from . import distances
from .base import AnalysisBase

from ..due import due, Doi

due.cite(Doi("10.1002/jcc.21787"),
         description="LeafletFinder algorithm",
         path="MDAnalysis.analysis.leaflet",
         cite_module=True)
del Doi


class LeafletFinder(object):
    """Identify atoms in the same leaflet of a lipid bilayer.

    This class implements the *LeafletFinder* algorithm [Michaud-Agrawal2011]_.

    Parameters
    ----------
    universe : Universe or str
        :class:`MDAnalysis.Universe` or a file name (e.g., in PDB or
        GRO format)
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

    Example
    -------
    The components of the graph are stored in the list
    :attr:`LeafletFinder.components`; the atoms in each component are numbered
    consecutively, starting at 0. To obtain the atoms in the input structure
    use :meth:`LeafletFinder.groups`::

       L = LeafletFinder(PDB, 'name P*')
       leaflet0 = L.groups(0)
       leaflet1 = L.groups(1)

    The residues can be accessed through the standard MDAnalysis mechanism::

       leaflet0.residues

    provides a :class:`~MDAnalysis.core.groups.ResidueGroup`
    instance. Similarly, all atoms in the first leaflet are then ::

       leaflet0.residues.atoms

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`
    """

    def __init__(self, universe, select, cutoff=15.0, pbc=False, sparse=None):
        universe = core.universe.as_Universe(universe)
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
                adj = distances.contact_matrix(coord, cutoff=self.cutoff, returntype="numpy", box=box)
            except ValueError:      # pragma: no cover
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
            except ValueError:       # pragma: no cover
                # but use a sparse matrix method for larger systems for memory reasons
                warnings.warn(
                    'N x N matrix too big - switching to sparse matrix method (works fine, but is currently rather '
                    'slow)',
                    category=UserWarning, stacklevel=2)
                adj = distances.contact_matrix(coord, cutoff=self.cutoff, returntype="sparse", box=box)
        return NX.Graph(adj)

    def _get_components(self):
        """Return connected components (as sorted numpy arrays), sorted by size."""
        return [np.sort(list(component)) for component in NX.connected_components(self.graph)]

    def update(self, cutoff=None):
        """Update components, possibly with a different *cutoff*"""
        if cutoff is None:
            cutoff = self.cutoff
        self._init_graph(cutoff)

    def sizes(self):
        """Dict of component index with size of component."""
        return dict(((idx, len(component)) for idx, component in enumerate(self.components)))

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
        sw = selections.get_writer(filename, kwargs.pop('format', None))
        with sw(filename, mode=kwargs.pop('mode', 'w'),
                preamble="leaflets based on select={selectionstring!r} cutoff={cutoff:f}\n".format(
                    **vars(self)),
                **kwargs) as writer:
            for i, ag in enumerate(self.groups_iter()):
                name = "leaflet_{0:d}".format((i + 1))
                writer.write(ag, name=name)

    def __repr__(self):
        return "<LeafletFinder({0!r}, cutoff={1:.1f} A) with {2:d} atoms in {3:d} groups>".format(
            self.selectionstring, self.cutoff, self.selection.n_atoms,
            len(self.components))


def optimize_cutoff(universe, select, dmin=10.0, dmax=20.0, step=0.5,
                    max_imbalance=0.2, **kwargs):
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
    kwargs.pop('cutoff', None)  # not used, so we filter it
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

class LipidEnrichment(AnalysisBase):
    """Calculate the lipid depletion-enrichment index around a protein
    by leaflet.

    The depletion-enrichment index (DEI) of a lipid around a protein indicates
    how enriched or depleted that lipid in the lipid annulus around the
    protein, with respect to the density of the lipid in the bulk
    membrane. If more than one leaflet is specified, then the index is
    calculated for each leaflet.

    The DEI for lipid species or group :math:`L` at cutoff :math:`r` is as follows:

    .. math::

        \text{DEI}_{L, r} = \frac{n(x_{(L, r)})}{n(x_r)} \times \frac{n(x)}{n(x_L)}

    where :math:`n(x_{(L, r)})` is the number of lipids :math:`L` within 
    distance :math:`r` of the protein; :math:`n(x_r)` is total number of lipids 
    within distance :math:`r` of the protein; :math:`n(x_L)` is the total number 
    of lipids :math:`L` in that membrane or leaflet; and :math:`n(x)` is the total 
    number of all lipids in that membrane or leaflet.

    The results of this analysis contain these values:

    * 'Near protein': the number of lipids :math:`L` within 
        cutoff :math:`r` of the protein
    * 'Fraction near protein': the fraction of the lipids :math:`L`
        with respect to all lipids within cutoff :math:`r` of the 
        protein: :math:`\frac{n(x_{(L, r)})}{n(x_r)}`
    * 'Enrichment': the depletion-enrichment index.

    This algorithm was obtained from [Corradi2018]_ . Please cite them if you use 
    this analysis in published work.

    .. note::

        This analysis requires ``scikit-learn`` to be installed to 
        determine leaflets.


    Parameters
    ----------

    universe: Universe or AtomGroup
        The atoms to apply this analysis to.
    select_protein: str (optional)
        Selection string for the protein.
    select_residues: str (optional)
        Selection string for the group of residues / lipids to analyse.
    select_headgroup: str (optional)
        Selection string for the lipid headgroups. This is used to determine
        leaflets and compute distances, unless ``compute_headgroup_only`` is
        ``False``.
    n_leaflets: int (optional)
        How many leaflets to split the membrane into.
    count_by_attr: str (optional)
        How to split the lipid species. By default, this is `resnames` so
        each lipid species is computed individually. However, any topology
        attribute can be used. For example, you could add a new attribute
        describing tail saturation.
    enrichment_cutoff: float (optional)
        Cutoff in ångström
    delta: float (optional)
        Leaflets are determined via spectral clustering on a similarity matrix.
        This similarity matrix is created from a pairwise distance matrix by
        applying the gaussian kernel 
        :math:`\exp{\frac{-\text{dist_matrix}^2}{2*\delta^2}}`. If your
        leaflets are incorrect, try finetuning the `delta` value.
    compute_headgroup_only: bool (optional)
        By default this analysis only uses distance from the lipid headgroup
        to determine whether a lipid is within the `enrichment_cutoff`
        distance from the protein. This choice was made to save computation.
        Compute the distance from every atom by toggling this option off.
    **kwargs
        Passed to :class:`~MDAnalysis.analysis.base.AnalysisBase`.

    Attributes
    ----------

    protein: :class:`~MDAnalysis.core.groups.AtomGroup`
        Protein atoms.
    residues: :class:`~MDAnalysis.core.groups.ResidueGroup`
        Lipid residues.
    headgroups: :class:`~MDAnalysis.core.groups.AtomGroup`
        Headgroup atoms.
    leaflet_residues: list of :class:`~MDAnalysis.core.groups.ResidueGroup`
        Residues per leaflet.
    leaflet_headgroups: list of :class:`~MDAnalysis.core.groups.AtomGroup`
        Headgroup atoms per leaflet.
    residue_counts: list of :class:`numpy.ndarray` (n_lipid_groups, n_frames)
        Number of each residue around protein for each frame, per leaflet.
    total_counts: list of :class:`numpy.ndarray` (n_frames,)
        Total number of residues around protein for each frame, per leaflet.
    leaflets: list of dict of dicts
        Counts, fraction, and enrichment index of each residue per frame, per
        leaflet.
    leaflets_summary: list of dict of dicts
        Counts, fraction, and enrichment index as averages and standard 
        deviations, per leaflet.
    box: :class:`numpy.ndarray` or ``None``
        Universe dimensions.
    n_leaflets: int
        Number of leaflets
    delta: float
        delta used Gaussian kernel to transform pairwise distances into a
        similarity matrix for clustering.
    compute_headgroup_only: bool
        whether to compute distances only using the headgroup.
    attrname: str
        The topology attribute used to group lipid species.
    """

    def __init__(self, universe, select_protein: str='protein', 
                 select_residues: str='all', select_headgroup: str='name PO4',
                 n_leaflets: int=2, count_by_attr:str ='resnames', 
                 enrichment_cutoff: float=6, delta: float=0.5,
                 compute_headgroup_only: bool=True, **kwargs):
        super(LipidEnrichment, self).__init__(universe.universe.trajectory, 
                                              **kwargs)
        self.box = universe.dimensions
        self.protein = universe.select_atoms(select_protein)
        self.residues = universe.select_atoms(select_residues).residues
        self.n_residues = len(self.residues)
        self.headgroups = self.residues.atoms.select_atoms(select_headgroup)
        if n_leaflets < 1:
            raise ValueError('Must have at least one leaflet')
        self.n_leaflets = n_leaflets
        self.delta = delta
        self.compute_headgroup_only = compute_headgroup_only
        self.cutoff = enrichment_cutoff
        self.leaflets = []
        self.leaflets_summary = []

        # process this a bit to remove errors like "alt_locs"
        attrname = _TOPOLOGY_ATTRNAMES[count_by_attr.lower().replace('_', '')]
        self.attrname = _TOPOLOGY_ATTRS[attrname].attrname

    def _prepare(self):
        try:
            import sklearn.cluster as skc
        except ImportError:
            raise ImportError('scikit-learn is required to use this analysis '
                              'but is not installed. Install it with `conda '
                              'install scikit-learn` or `pip install '
                              'scikit-learn`.') from None

        if self.n_leaflets > 1:
            # determine leaflets by clustering center-of-geometry
            coms = self.headgroups.center_of_geometry(compound='residues')
            pdist = distances.distance_array(coms, coms, box=self.box)
            psim = np.exp(-pdist**2/(2*(self.delta**2)))
            sc = skc.SpectralClustering(n_clusters=self.n_leaflets, affinity='precomputed_nearest_neighbors')
            clusters = sc.fit_predict(pdist)
            leaflets = []
            for i in range(self.n_leaflets):
                leaflets.append([])
            for res, i in zip(self.residues, clusters):
                leaflets[i].append(res)
            self.leaflet_residues = lres = [np.array(x).sum() 
                                            for x in leaflets]
            self.leaflet_residues.sort(key=lambda x: x.center_of_geometry()[-1],
                                       reverse=True)
            self.leaflet_headgroups = [(self.headgroups & r.atoms) 
                                       for r in lres]
        else:
            self.leaflet_residues = [self.residues]
            self.leaflet_headgroups = [self.headgroups]
        
        if not self.compute_headgroup_only:
            self.leaflet_headgroups = [h.residues.atoms
                                       for h in self.leaflet_headgroups]

        # set up results
        self.residue_counts = [np.zeros((len(r), self.n_frames)) for r in lres]
        self.total_counts = np.zeros((self.n_leaflets, self.n_frames))
        self.leaflet_ids = ids = [getattr(r, self.attrname) for r in lres]
        self.leaflets = [dict.fromkeys(np.unique(r)) for r in ids]
        self.leaflets_summary = [dict.fromkeys(np.unique(r)) for r in ids]

    def _single_frame(self):
        for i, headgroup in enumerate(self.leaflet_headgroups):
            hcom = headgroup.center_of_geometry(compound='residues')
            pairs = distances.capped_distance(self.protein.positions,
                                              hcom, self.cutoff,
                                              box=self.box,
                                              return_distances=False)
            if pairs.size > 0:
                indices = np.sort(pairs[:, 1])
                indices = np.unique(indices)
            else:
                indices = []
            self.residue_counts[i][indices, self._frame_index] = 1
            self.total_counts[i][self._frame_index] = len(indices)
    
    def _conclude(self):
        for i, leaf in enumerate(self.leaflets):
            ids = self.leaflet_ids[i]
            summary = self.leaflets_summary[i]
            res_counts = self.residue_counts[i]
            all_lipids = {}
            all_lipids_sum = {}
            all_lipids['Near protein'] = nc_all = res_counts.sum(axis=0)
            all_lipids_sum['Total near protein'] = nc_all.sum()
            all_lipids_sum['Average near protein'] = avgall = nc_all.mean()
            all_lipids_sum['Total number'] = all_total = len(res_counts)
            all_lipids_sum['SD near protein'] = nc_all.std()
            for resname in leaf.keys():
                mask = ids == resname
                counts = res_counts[mask]
                results = {}
                results_sum = {}
                results['Near protein'] = nc = counts.sum(axis=0)
                results_sum['Average near protein'] = avg = nc.mean()
                results_sum['Total number'] = total = sum(mask)
                results_sum['SD near protein'] = nc.std()
                results['Fraction near protein'] = frac = nc/nc_all
                results_sum['Average fraction near protein'] = avg_frac = avg/avgall
                results_sum['SD fraction near protein'] = frac.std()
                denom = total/all_total
                results['Enrichment'] = frac/denom
                results_sum['Average enrichment'] = avg_frac/denom
                results_sum['SD enrichment'] = (frac/denom).std()
                leaf[resname] = results
                summary[resname] = results_sum
            leaf['all'] = all_lipids
            summary['all'] = all_lipids_sum
    
    def summary_as_dataframe(self):
        """Convert the results summary into a pandas DataFrame.
        
        This requires pandas to be installed.
        """

        if not self.leaflets_summary:
            raise ValueError('Call run() first to get results')
        try:
            import pandas as pd
        except ImportError:
            raise ImportError('pandas is required to use this function '
                              'but is not installed. Please install with '
                              '`conda install pandas` or '
                              '`pip install pandas`.') from None
        
        dfs = [pd.DataFrame.from_dict(d, orient='index') 
               for d in self.leaflets_summary]
        for i, df in enumerate(dfs, 1):
            df['Leaflet'] = i
        df = pd.concat(dfs)
        return df
