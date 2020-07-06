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

.. autoclass:: LipidEnrichment
    :members:
    
"""
import warnings
import functools

import numpy as np
import networkx as NX

from .. import core, selections, _TOPOLOGY_ATTRNAMES, _TOPOLOGY_ATTRS
from . import distances
from .base import AnalysisBase

from ..due import due, Doi

due.cite(Doi("10.1002/jcc.21787"),
         description="LeafletFinder 'graph' algorithm",
         path="MDAnalysis.analysis.leaflet.LeafletFinder",
         cite_module=True)

due.cite(Doi("10.1021/acscentsci.8b00143"),
         description="Lipid enrichment-depletion index formula",
         path="MDAnalysis.analysis.leaflet.LipidEnrichment",
         cite_module=True)

del Doi


class LeafletFinder(object):
    """Identify atoms in the same leaflet of a lipid bilayer.

    This class implements the *LeafletFinder* algorithm [Michaud-Agrawal2011]_.

    Parameters
    ----------
    universe : Universe or AtomGroup
        Atoms to apply the algorithm to
    select : str
        A :meth:`Universe.select_atoms` selection string
        for atoms that define the lipid head groups, e.g.
        universe.atoms.PO4 or "name PO4" or "name P*"
    cutoff : float (optional)
        cutoff distance for computing distances (for the spectral clustering
        method) or determining connectivity in the same leaflet (for the graph
        method). In spectral clustering, it just has to be suitably large to
        cover a significant part of the leaflet, but lower values increase
        computational efficiency. Please see the :func:`optimize_cutoff`
        function for help in with values for the graph method.
    pbc : bool (optional)
        If ``False``, does not follow the minimum image convention when
        computing distances
    method: str, {"graph", "spectralclustering"}
        method to use to assign groups to leaflets. Choose either 
        "graph" for :func:`~distances.group_coordinates_by_graph` or
        "spectralclustering" for 
        :func:`~distances.group_coordinates_by_spectralclustering`
    **kwargs:
        Passed to ``method``

    Attributes
    ----------
    universe: Universe
    select: str
        Selection string
    selection: AtomGroup
        Atoms that the analysis is applied to
    pbc: bool
        Whether to use PBC or not
    box: numpy.ndarray or None
        Cell dimensions to use in calculating distances
    predictor: :class:`networkx.Graph` or :class:`sklearn.cluster.SpectralClustering`
        The object used to predict the leaflets
    components: list of numpy.ndarray
        List of indices of atoms in each leaflet, corresponding to the
        order of `selection`. ``components[i]`` is the array of indices
        for the ``i``-th leaflet. ``k = components[i][j]`` means that the
        ``k``-th atom in `selection` is in the ``i``-th leaflet. 
        The components are sorted by size.
    groups: list of AtomGroups
        List of AtomGroups in each leaflet. ``groups[i]`` is the ``i``-th
        leaflet. The leaflets are sorted by size.
    leaflets: list of AtomGroup
        List of AtomGroups in each leaflet. ``groups[i]`` is the ``i``-th
        leaflet. The leaflets are sorted by z-coordinate so that the
        upper-most leaflet is first.
    sizes: list of ints
        List of the size of each leaflet in ``groups``.


    Example
    -------
    The components of the graph are stored in the list
    :attr:`LeafletFinder.components`; the atoms in each component are numbered
    consecutively, starting at 0. To obtain the atoms in the input structure
    use :attr:`LeafletFinder.groups`::

       L = LeafletFinder(PDB, 'name P*')
       leaflet_1 = L.groups[0]
       leaflet_2 = L.groups[1]

    The residues can be accessed through the standard MDAnalysis mechanism::

       leaflet_1.residues

    provides a :class:`~MDAnalysis.core.groups.ResidueGroup`
    instance. Similarly, all atoms in the first leaflet are then ::

       leaflet_1.residues.atoms


    See also
    --------
    :func:`~MDAnalysis.analysis.distances.group_coordinates_by_graph`
    :func:`~MDAnalysis.analysis.distances.group_coordinates_by_spectralclustering`


    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`
    """

    def __init__(self, universe, select='all', cutoff=15.0, pbc=True,
                 method="graph", **kwargs):
        method = method.lower().replace('_', '')
        if method == "graph":
            self.method = distances.group_coordinates_by_graph
        elif method == "spectralclustering":
            self.method = distances.group_coordinates_by_spectralclustering
        elif method == "centerofgeometry":
            self.method = distances.group_coordinates_by_cog
        else:
            raise ValueError("`method` must be in {'graph', "
                             "'spectralclustering', 'center_of_geometry'")

        self.universe = universe.universe
        self.select = select
        self.selection = universe.select_atoms(select, periodic=pbc)
        self.pbc = pbc
        self.cutoff = cutoff
        self.box = self.universe.dimensions if pbc else None
        results = self.method(self.selection.positions,
                              cutoff=self.cutoff,
                              box=self.box,
                              return_predictor=True,
                              **kwargs)
        self.components, self.predictor = results
        self.groups = [self.selection[x] for x in self.components]
        self.leaflets = sorted(self.groups,
                               key=lambda x: x.center_of_geometry()[-1],
                               reverse=True)
        self.sizes = [len(ag) for ag in self.groups]

    def groups_iter(self):
        """Iterator over all leaflet :meth:`groups`"""
        for group in self.groups:
            yield group

    def write_selection(self, filename, mode="w", format=None, **kwargs):
        """Write selections for the leaflets to *filename*.

        The format is typically determined by the extension of *filename*
        (e.g. "vmd", "pml", or "ndx" for VMD, PyMol, or Gromacs).

        See :class:`MDAnalysis.selections.base.SelectionWriter` for all
        options.
        """
        sw = selections.get_writer(filename, format)
        with sw(filename, mode=mode,
                preamble=f"Leaflets found by {repr(self)}\n",
                **kwargs) as writer:
            for i, ag in enumerate(self.groups, 1):
                writer.write(ag, name=f"leaflet_{i:d}")

    def __repr__(self):
        return (f"LeafletFinder(select='{self.select}', "
                f"cutoff={self.cutoff:.1f} Å, pbc={self.pbc})")


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
        sizes = LF.sizes
        if len(sizes) < 2:
            continue
        n0 = float(sizes[0])  # sizes of two biggest groups ...
        n1 = float(sizes[1])  # ... assumed to be the leaflets
        imbalance = np.abs(n0 - n1) / (n0 + n1)
        # print "sizes: %(sizes)r; imbalance=%(imbalance)f" % vars()
        if imbalance > max_imbalance:
            continue
        _sizes.append((cutoff, len(LF.sizes)))
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

    This algorithm was obtained from [Corradi2018]_. Please cite them if you use 
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

    def __init__(self, universe, select_protein: str = 'protein',
                 select_residues: str = 'all', select_headgroup: str = 'name PO4',
                 n_leaflets: int = 2, count_by_attr: str = 'resnames',
                 enrichment_cutoff: float = 6, pbc: bool = True,
                 compute_headgroup_only: bool = True,
                 update_leaflet_step: int = 1,
                 group_method: str = "spectralclustering",
                 update_method=None, group_kwargs: dict = {},
                 update_kwargs: dict = {}, **kwargs):
        super(LipidEnrichment, self).__init__(universe.universe.trajectory,
                                              **kwargs)
        if n_leaflets < 1:
            raise ValueError('Must have at least one leaflet')
        self.n_leaflets = n_leaflets
        group_method = group_method.lower().replace('_', '')
        if group_method in ("spectralclustering", "graph"):
            self._get_leaflets = functools.partial(self._update_leafletfinder,
                                                   method=group_method,
                                                   kwargs=group_kwargs)
        else:
            raise ValueError("`group_method` should be one of "
                             "{'spectralclustering', 'graph'}")

        if update_method is None:
            self._update_leaflets = self._get_leaflets
        else:
            update_method = update_method.lower().replace('_', '')
            if update_method in ("spectralclustering", "graph"):
                self._update_leaflets = functools.partial(self._update_leafletfinder,
                                                          method=update_method,
                                                          kwargs=update_kwargs)
            elif update_method == "centerofgeometry":
                self._update_leaflets = self._update_cog
            else:
                raise ValueError("`update_method` should be one of "
                                 "{'spectralclustering', 'graph', "
                                 "'center_of_geometry'}")

        self.group_kwargs = group_kwargs
        self.update_leaflet_step = update_leaflet_step

        self.pbc = pbc
        self.box = universe.dimensions if pbc else None
        self.protein = universe.select_atoms(select_protein)
        self.residues = universe.select_atoms(select_residues).residues
        self._resindices = self.residues.resindices
        self.n_residues = len(self.residues)
        self.headgroups = self.residues.atoms.select_atoms(select_headgroup)

        if compute_headgroup_only:
            self._residue_groups = [r.atoms.select_atoms(select_headgroup)
                                    for r in self.residues]
            self._compute_atoms = self.headgroups
        else:
            self._residue_groups = self.residues
            self._compute_atoms = self.residues.atoms
        self.cutoff = enrichment_cutoff
        self.leaflets = []
        self.leaflets_summary = []

        # process this a bit to remove errors like "alt_locs"
        attrname = _TOPOLOGY_ATTRNAMES[count_by_attr.lower().replace('_', '')]
        self.attrname = _TOPOLOGY_ATTRS[attrname].attrname

    def _prepare(self):
        self.ids = np.unique(getattr(self.residues, self.attrname))
        self.near_counts = np.zeros((self.n_leaflets, len(self.ids),
                                     self.n_frames))
        self.residue_counts = np.zeros((self.n_leaflets, len(self.ids),
                                        self.n_frames))
        self.total_counts = np.zeros((self.n_leaflets, self.n_frames))
        self.leaflet_residue_masks = np.zeros((self.n_frames, self.n_leaflets,
                                               self.n_residues), dtype=bool)
        self.leaflet_residues = np.zeros((self.n_frames, self.n_leaflets),
                                         dtype=object)
        self._get_leaflets()
        self._current_ids = [getattr(r, self.attrname)
                             for r in self._current_leaflets]

    def _single_frame(self):
        if not self._frame_index % self.update_leaflet_step:
            self._update_leaflets()
            self._current_ids = [getattr(r, self.attrname)
                                 for r in self._current_leaflets]

        hcom = self._compute_atoms.center_of_geometry(compound="residues")
        pairs = distances.capped_distance(self.protein.positions, hcom,
                                          self.cutoff, box=self.box,
                                          return_distances=False)
        if pairs.size > 0:
            indices = np.sort(pairs[:, 1])
            indices = np.unique(indices)
        else:
            indices = []

        resix = self._resindices[indices]

        for i, leaf in enumerate(self._current_leaflets):
            ids = self._current_ids[i]
            _, ix1, ix2 = np.intersect1d(resix, leaf.resindices,
                                         assume_unique=True,
                                         return_indices=True)
            self.total_counts[i, self._frame_index] = len(ix1)
            subids = ids[ix2]
            for j, x in enumerate(self.ids):
                self.residue_counts[i, j, self._frame_index] = sum(ids == x)
                self.near_counts[i, j, self._frame_index] = sum(subids == x)

    def _conclude(self):
        for i in range(self.n_leaflets):
            leaf = {}
            summary = {}
            res_counts = self.residue_counts[i]
            near_counts = self.near_counts[i]
            all_lipids = {}
            all_lipids_sum = {}
            all_lipids['Near protein'] = nc_all = near_counts.sum(axis=0)
            all_lipids['Total number'] = tc_all = res_counts.sum(axis=0)
            all_lipids_sum['Total near protein'] = nc_all.sum()
            all_lipids_sum['Average near protein'] = avgall = nc_all.mean()
            all_lipids_sum['Total number'] = all_total = res_counts.sum()
            all_lipids_sum['SD near protein'] = nc_all.std()

            for j, resname in enumerate(self.ids):
                results = {}
                results_sum = {}
                results['Near protein'] = nc = near_counts[j]
                results_sum['Average near protein'] = avg = nc.mean()
                results_sum['Total number'] = total = res_counts[j]
                results_sum['SD near protein'] = nc.std()
                results['Fraction near protein'] = frac = nc / nc_all
                avg_frac = avg / avgall
                results_sum['Average fraction near protein'] = avg_frac
                results_sum['SD fraction near protein'] = frac.std()
                results['Enrichment'] = frac / (total / tc_all)
                denom = total.sum() / all_total
                results_sum['Average enrichment'] = avg_frac / denom
                results_sum['SD enrichment'] = (frac / denom).std()
                leaf[resname] = results
                summary[resname] = results_sum
            leaf['all'] = all_lipids
            summary['all'] = all_lipids_sum
            self.leaflets.append(leaf)
            self.leaflets_summary.append(summary)

    def _update_leafletfinder(self, method="spectralclustering", kwargs={}):
        lf = LeafletFinder(self.headgroups, **kwargs,
                           pbc=self.pbc, method=method)
        self._current_leaflets = lf.leaflets

    def _update_cog(self):
        # this one relies on the leaflets not changing _too_ much. Fragile.
        lcogs = [x.center_of_geometry() for x in self._current_leaflets]
        rcogs = self.headgroups.center_of_geometry(compound='residues')
        ix = distances.group_coordinates_by_cog(rcogs, lcogs, box=self.box,
                                                return_predictor=False)
        self._current_leaflets = [self._residue_groups[r] for r in ix]

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
