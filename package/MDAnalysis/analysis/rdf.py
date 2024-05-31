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

r"""Radial Distribution Functions --- :mod:`MDAnalysis.analysis.rdf`
====================================================================

This module contains two classes to calculate radial
`pair distribution functions`_ (`radial distribution functions`_ or "RDF").
The RDF :math:`g_{ab}(r)` between types of particles :math:`a` and :math:`b` is

.. _equation-gab:

.. math::

   g_{ab}(r) = \frac{1}{N_{a}} \frac{1}{N_{b}/V} \sum_{i=1}^{N_a} \sum_{j=1}^{N_b}
               \langle \delta(|\mathbf{r}_i - \mathbf{r}_j| - r) \rangle

which is normalized so that the RDF becomes 1 for large separations in a
homogenous system. The RDF effectively counts the average number of :math:`b`
neighbours in a shell at distance :math:`r` around a :math:`a` particle and
represents it as a density.

The radial cumulative distribution function is

.. math::

   G_{ab}(r) = \int_0^r \!\!dr' 4\pi r'^2 g_{ab}(r')

and the average number of :math:`b` particles within radius :math:`r`

.. _equation-countab:

.. math::

   N_{ab}(r) = \rho G_{ab}(r)

(with the appropriate density :math:`\rho`). The latter function can be used to
compute, for instance, coordination numbers such as the number of neighbors in
the first solvation shell :math:`N(r_1)` where :math:`r_1` is the position of
the first minimum in :math:`g(r)`.

We provide options for calculating the density of particle :math:`b`
in a shell at distance :math:`r` around a :math:`a` particle, which is

.. _equation-nab:

.. math::
   n_{ab}(r) = \rho g_{ab}(r)

.. _`pair distribution functions`:
   https://en.wikipedia.org/wiki/Pair_distribution_function
.. _`radial distribution functions`:
   https://en.wikipedia.org/wiki/Radial_distribution_function

.. Not Implemented yet:
.. - Structure factor?
.. - Coordination number
"""
import warnings
import numpy as np

from ..lib import distances
from .base import AnalysisBase


class InterRDF(AnalysisBase):
    r"""Radial distribution function

    :class:`InterRDF` is a tool to calculate average radial distribution
    functions between two groups of atoms. Suppose we have two AtomGroups ``A``
    and ``B``. ``A`` contains atom ``A1``, ``A2``, and ``B`` contains ``B1``,
    ``B2``. Given ``A`` and ``B`` to :class:`InterRDF`, the output will be the
    average of RDFs between ``A1`` and ``B1``, ``A1`` and ``B2``, ``A2`` and
    ``B1``, ``A2`` and ``B2``. A typical application is to calculate the RDF of
    solvent with itself or with another solute.

    The :ref:`radial distribution function<equation-gab>` is calculated by
    histogramming distances between all particles in `g1` and `g2` while taking
    periodic boundary conditions into account via the minimum image
    convention.

    The `exclusion_block` keyword may be used to exclude a set of distances
    from the calculations.

    Results are available in the attributes :attr:`results.rdf`
    and :attr:`results.count`.

    Parameters
    ----------
    g1 : AtomGroup
        First AtomGroup
    g2 : AtomGroup
        Second AtomGroup
    nbins : int
        Number of bins in the histogram
    range : tuple or list
        The size of the RDF
    norm : str, {'rdf', 'density', 'none'}
          For 'rdf' calculate :math:`g_{ab}(r)`. For
          'density' the :ref:`single particle density<equation-nab>`
          :math:`n_{ab}(r)` is computed. 'none' computes the number of
          particles occurences in each spherical shell.

          .. versionadded:: 2.3.0

    exclusion_block : tuple
        A tuple representing the tile to exclude from the distance array.
    exclude_same : str
        Will exclude pairs of atoms that share the same "residue", "segment", or "chain".
        Those are the only valid values. This is intended to remove atoms that are
        spatially correlated due to direct bonded connections.
    verbose : bool
        Show detailed progress of the calculation if set to `True`

    Attributes
    ----------
    results.bins : numpy.ndarray
       :class:`numpy.ndarray` of the centers of the `nbins` histogram
       bins.

       .. versionadded:: 2.0.0

    bins : numpy.ndarray
       Alias to the :attr:`results.bins` attribute.

       .. deprecated:: 2.0.0
           This attribute will be removed in 3.0.0.
           Use :attr:`results.bins` instead.

    results.edges : numpy.ndarray

      :class:`numpy.ndarray` of the `nbins + 1` edges of the histogram
      bins.

       .. versionadded:: 2.0.0

    edges : numpy.ndarray

       Alias to the :attr:`results.edges` attribute.

       .. deprecated:: 2.0.0
           This attribute will be removed in 3.0.0.
           Use :attr:`results.edges` instead.

    results.rdf : numpy.ndarray
      :class:`numpy.ndarray` of the :ref:`radial distribution
      function<equation-gab>` values for the :attr:`results.bins`.

       .. versionadded:: 2.0.0

    rdf : numpy.ndarray
       Alias to the :attr:`results.rdf` attribute.

       .. deprecated:: 2.0.0
           This attribute will be removed in 3.0.0.
           Use :attr:`results.rdf` instead.

    results.count : numpy.ndarray
      :class:`numpy.ndarray` representing the radial histogram, i.e.,
      the raw counts, for all :attr:`results.bins`.

       .. versionadded:: 2.0.0

    count : numpy.ndarray
       Alias to the :attr:`results.count` attribute.

       .. deprecated:: 2.0.0
           This attribute will be removed in 3.0.0.
           Use :attr:`results.count` instead.

    Example
    -------
    First create the :class:`InterRDF` object, by supplying two
    AtomGroups then use the :meth:`run` method ::

      rdf = InterRDF(ag1, ag2)
      rdf.run()

    Results are available through the :attr:`results.bins` and
    :attr:`results.rdf` attributes::

      plt.plot(rdf.results.bins, rdf.results.rdf)

    The `exclusion_block` keyword allows the masking of pairs from
    within the same molecule. For example, if there are 7 of each
    atom in each molecule, the exclusion mask ``(7, 7)`` can be used.


    .. versionadded:: 0.13.0

    .. versionchanged:: 1.0.0
       Support for the `start`, `stop`, and `step` keywords has been
       removed. These should instead be passed to :meth:`InterRDF.run`.

    .. versionchanged:: 2.0.0
       Store results as attributes `bins`, `edges`, `rdf` and `count`
       of the `results` attribute of
       :class:`~MDAnalysis.analysis.AnalysisBase`.
    """
    def __init__(self,
                 g1,
                 g2,
                 nbins=75,
                 range=(0.0, 15.0),
                 norm="rdf",
                 exclusion_block=None,
                 exclude_same=None,
                 **kwargs):
        super(InterRDF, self).__init__(g1.universe.trajectory, **kwargs)
        self.g1 = g1
        self.g2 = g2
        self.norm = str(norm).lower()

        self.rdf_settings = {'bins': nbins,
                             'range': range}
        self._exclusion_block = exclusion_block
        if exclude_same is not None and exclude_same not in ['residue', 'segment', 'chain']:
            raise ValueError(
                "The exclude_same argument to InterRDF must be None, 'residue', 'segment' "
                "or 'chain'."
            )
        if exclude_same is not None and exclusion_block is not None:
            raise ValueError(
                "The exclude_same argument to InterRDF cannot be used with exclusion_block."
            )
        name_to_attr = {'residue': 'resindices', 'segment': 'segindices', 'chain': 'chainIDs'}
        self.exclude_same = name_to_attr.get(exclude_same)

        if self.norm not in ['rdf', 'density', 'none']:
            raise ValueError(f"'{self.norm}' is an invalid norm. "
                             "Use 'rdf', 'density' or 'none'.")

    def _prepare(self):
        # Empty histogram to store the RDF
        count, edges = np.histogram([-1], **self.rdf_settings)
        count = count.astype(np.float64)
        count *= 0.0
        self.results.count = count
        self.results.edges = edges
        self.results.bins = 0.5 * (edges[:-1] + edges[1:])

        if self.norm == "rdf":
            # Cumulative volume for rdf normalization
            self.volume_cum = 0
        # Set the max range to filter the search radius
        self._maxrange = self.rdf_settings['range'][1]

    def _single_frame(self):
        pairs, dist = distances.capped_distance(self.g1.positions,
                                                self.g2.positions,
                                                self._maxrange,
                                                box=self._ts.dimensions)
        # Maybe exclude same molecule distances
        if self._exclusion_block is not None:
            idxA = pairs[:, 0]//self._exclusion_block[0]
            idxB = pairs[:, 1]//self._exclusion_block[1]
            mask = np.where(idxA != idxB)[0]
            dist = dist[mask]

        if self.exclude_same is not None:
            # Ignore distances between atoms in the same attribute
            attr_ix_a = getattr(self.g1, self.exclude_same)[pairs[:, 0]]
            attr_ix_b = getattr(self.g2, self.exclude_same)[pairs[:, 1]]
            mask = np.where(attr_ix_a != attr_ix_b)[0]
            dist = dist[mask]

        count, _ = np.histogram(dist, **self.rdf_settings)
        self.results.count += count

        if self.norm == "rdf":
            self.volume_cum += self._ts.volume

    def _conclude(self):
        norm = self.n_frames
        if self.norm in ["rdf", "density"]:
            # Volume in each radial shell
            vols = np.power(self.results.edges, 3)
            norm *= 4/3 * np.pi * np.diff(vols)

        if self.norm == "rdf":
            # Number of each selection
            nA = len(self.g1)
            nB = len(self.g2)
            N = nA * nB

            # If we had exclusions, take these into account
            if self._exclusion_block:
                xA, xB = self._exclusion_block
                nblocks = nA / xA
                N -= xA * xB * nblocks

            # Average number density
            box_vol = self.volume_cum / self.n_frames
            norm *= N / box_vol

        self.results.rdf = self.results.count / norm

    @property
    def edges(self):
        wmsg = ("The `edges` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.bins` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.edges

    @property
    def count(self):
        wmsg = ("The `count` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.bins` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.count

    @property
    def bins(self):
        wmsg = ("The `bins` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.bins` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.bins

    @property
    def rdf(self):
        wmsg = ("The `rdf` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.rdf` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.rdf


class InterRDF_s(AnalysisBase):
    r"""Site-specific radial distribution function

    Calculates site-specific radial distribution
    functions. Instead of two groups of atoms it takes as input a list of
    pairs of AtomGroup, ``[[A, B], [C, D], ...]``. Given the same ``A`` and
    ``B`` to
    :class:`InterRDF_s`, the output will be a list of individual RDFs between
    ``A1`` and ``B1``, ``A1`` and ``B2``, ``A2`` and ``B1``, ``A2`` and ``B2``
    (and
    similarly for ``C`` and ``D``). These site-specific radial distribution
    functions are typically calculated if one is interested in the solvation
    shells of a ligand in a binding site or the solvation of specific residues
    in a protein.

    Parameters
    ----------
    u : Universe
        a Universe that contains atoms in `ags`

       .. deprecated:: 2.3.0
           This parameter is superflous and will be removed in
           MDAnalysis 3.0.0.

    ags : list
        a list of pairs of :class:`~MDAnalysis.core.groups.AtomGroup`
        instances
    nbins : int
        Number of bins in the histogram
    range : tuple or list
        The size of the RDF
    norm : str, {'rdf', 'density', 'none'}
        For 'rdf' calculate :math:`g_{ab}(r)`. For
        'density' the :ref:`single particle density<equation-nab>`
        :math:`n_{ab}(r)` is computed. 'none' computes the number of
        particles occurences in each spherical shell.

        .. versionadded:: 2.3.0

    density : bool
        `False`: calculate :math:`g_{ab}(r)`; `True`: calculate
        the true :ref:`single particle density<equation-nab>`
        :math:`n_{ab}(r)`. `density` overwrites the `norm` parameter.

        .. versionadded:: 1.0.1

            This keyword was available since 0.19.0 but was not
            documented. Furthermore, it had the opposite
            meaning. Since 1.0.1 it is officially supported as
            documented.

        .. deprecated:: 2.3.0
            Instead of `density=True` use `norm='density'`

    Attributes
    ----------
    results.bins : numpy.ndarray
        :class:`numpy.ndarray` of the centers of the `nbins` histogram
        bins; all individual site-specific RDFs have the same bins.

        .. versionadded:: 2.0.0

    bins : numpy.ndarray
        Alias to the :attr:`results.bins` attribute.

        .. deprecated:: 2.0.0
            This attribute will be removed in 3.0.0.
            Use :attr:`results.bins` instead.

    results.edges : numpy.ndarray
        array of the ``nbins + 1`` edges of the histogram
        bins; all individual site-specific RDFs have the same bins.

        .. versionadded:: 2.0.0

    edges : numpy.ndarray
        Alias to the :attr:`results.edges` attribute.

        .. deprecated:: 2.0.0
            This attribute will be removed in 3.0.0.
            Use :attr:`results.edges` instead.

    results.rdf : list
        :class:`list` of the site-specific :ref:`radial distribution
        functions<equation-gab>` if `norm='rdf'` or :ref:`density
        functions<equation-nab>` for the :attr:`bins`
        if `norm='density'`. The list contains
        ``len(ags)`` entries. Each entry for the ``i``-th pair `[A, B]
        = ags[i]` in `ags` is a :class:`numpy.ndarray` with shape
        ``(len(A), len(B))``, i.e., a stack of RDFs. For example,
        ``results.rdf[i][0, 2]`` is the RDF between atoms ``A[0]``
        and ``B[2]``.

        .. versionadded:: 2.0.0

    rdf : list
        Alias to the :attr:`results.rdf` attribute.

        .. deprecated:: 2.0.0
            This attribute will be removed in 3.0.0.
            Use :attr:`results.rdf` instead.

    results.count : list
        :class:`list` of the site-specific radial histograms, i.e., the
        raw counts, for all :attr:`results.bins`. The data have the same
        structure as :attr:`results.rdf` except that the arrays contain
        the raw counts.

        .. versionadded:: 2.0.0

    count : list
        Alias to the :attr:`results.count` attribute.

        .. deprecated:: 2.0.0
            This attribute will be removed in 3.0.0.
            Use :attr:`results.count` instead.

    results.cdf : list
        :class:`list` of the site-specific :ref:`cumulative
        counts<equation-countab>`, for all :attr:`results.bins`. The data
        have the same structure as :attr:`results.rdf` except that the arrays
        contain the cumulative counts.

        This attribute only exists after :meth:`get_cdf` has been run.

        .. versionadded:: 2.0.0

    cdf : list
        Alias to the :attr:`results.cdf` attribute.

        .. deprecated:: 2.0.0
            This attribute will be removed in 3.0.0.
            Use :attr:`results.cdf` instead.

    Example
    -------
    First create the :class:`InterRDF_s` object, by supplying one Universe and
    one list of pairs of AtomGroups, then use the :meth:`~InterRDF_s.run`
    method::

      from MDAnalysisTests.datafiles import GRO_MEMPROT, XTC_MEMPROT
      u = mda.Universe(GRO_MEMPROT, XTC_MEMPROT)

      s1 = u.select_atoms('name ZND and resid 289')
      s2 = u.select_atoms('(name OD1 or name OD2) and resid 51 and sphzone 5.0 (resid 289)')
      s3 = u.select_atoms('name ZND and (resid 291 or resid 292)')
      s4 = u.select_atoms('(name OD1 or name OD2) and sphzone 5.0 (resid 291)')
      ags = [[s1, s2], [s3, s4]]

      rdf = InterRDF_s(u, ags)
      rdf.run()

    Results are available through the :attr:`results.bins`
    and :attr:`results.rdf` attributes::

      plt.plot(rdf.results.bins, rdf.results.rdf[0][0, 0])

    (Which plots the rdf between the first atom in ``s1`` and the first atom in
    ``s2``)

    To generate the *cumulative distribution function* (cdf) in the sense of
    "particles within radius :math:`r`", i.e., :math:`N_{ab}(r)`, use the
    :meth:`~InterRDF_s.get_cdf` method ::

      cdf = rdf.get_cdf()

    Results are available through the :attr:`results.cdf` attribute::

      plt.plot(rdf.results.bins, rdf.results.cdf[0][0, 0])

    (Which plots the cdf between the first atom in ``s1`` and the first atom in
    ``s2``)


    .. versionadded:: 0.19.0

    .. versionchanged:: 1.0.0
       Support for the `start`, `stop`, and `step` keywords has been
       removed. These should instead be passed to :meth:`InterRDF_s.run`.

    .. versionchanged:: 2.0.0
       Store results as attributes `bins`, `edges`, `rdf`, `count`
       and `cdf` of the `results` attribute
       of :class:`~MDAnalysis.analysis.AnalysisBase`.

    .. versionchanged:: 2.3.0
       Introduce `norm` and `exclusion_blocks` attributes.
    .. deprecated:: 2.3.0
       Instead of `density=True` use `norm='density'`
    .. deprecated:: 2.3.0
       The `universe` parameter is superflous.
    """
    def __init__(self,
                 u,
                 ags,
                 nbins=75,
                 range=(0.0, 15.0),
                 norm="rdf",
                 density=False,
                 **kwargs):
        super(InterRDF_s, self).__init__(ags[0][0].universe.trajectory,
                                         **kwargs)

        warnings.warn("The `u` attribute is superflous and will be removed "
                      "in MDAnalysis 3.0.0.", DeprecationWarning)

        self.ags = ags
        self.norm = str(norm).lower()
        self.rdf_settings = {'bins': nbins,
                             'range': range}

        if self.norm not in ['rdf', 'density', 'none']:
            raise ValueError(f"'{self.norm}' is an invalid norm. "
                             "Use 'rdf', 'density' or 'none'.")

        if density:
            warnings.warn("The `density` attribute was deprecated in "
                          "MDAnalysis 2.3.0 and will be removed in "
                          "MDAnalysis 3.0.0. Please use `norm=density` "
                          "instead.", DeprecationWarning)
            self.norm = "density"

    def _prepare(self):
        count, edges = np.histogram([-1], **self.rdf_settings)
        self.results.count = [np.zeros((ag1.n_atoms, ag2.n_atoms, len(count)),
                                        dtype=np.float64) for ag1, ag2 in self.ags]
        self.results.edges = edges
        self.results.bins = 0.5 * (edges[:-1] + edges[1:])

        if self.norm == "rdf":
            # Cumulative volume for rdf normalization
            self.volume_cum = 0
        self._maxrange = self.rdf_settings['range'][1]

    def _single_frame(self):
        for i, (ag1, ag2) in enumerate(self.ags):
            pairs, dist = distances.capped_distance(ag1.positions,
                                                    ag2.positions,
                                                    self._maxrange,
                                                    box=self._ts.dimensions)

            for j, (idx1, idx2) in enumerate(pairs):
                count, _ = np.histogram(dist[j], **self.rdf_settings)
                self.results.count[i][idx1, idx2, :] += count

        if self.norm == "rdf":
            self.volume_cum += self._ts.volume

    def _conclude(self):
        norm = self.n_frames
        if self.norm in ["rdf", "density"]:
            # Volume in each radial shell
            vols = np.power(self.results.edges, 3)
            norm *= 4/3 * np.pi * np.diff(vols)

        if self.norm == "rdf":
            # Average number density
            norm *= 1 / (self.volume_cum / self.n_frames)

        # Empty lists to restore indices, RDF
        self.results.indices = []
        self.results.rdf = []

        for i, (ag1, ag2) in enumerate(self.ags):
            # Number of each selection
            self.results.indices.append([ag1.indices, ag2.indices])
            self.results.rdf.append(self.results.count[i] / norm)

    def get_cdf(self):
        r"""Calculate the cumulative counts for all sites.

        This is the :ref:`cumulative count<equation-countab>` within a given
        radius, i.e., :math:`N_{ab}(r)`.

        The result is returned and also stored in the attribute
        :attr:`results.cdf`.

        Returns
        -------
        cdf : list
              list of arrays with the same structure as :attr:`results.rdf`
        """
        self.results.cdf = []

        for count in self.results.count:
            self.results.cdf.append(np.cumsum(count, axis=2) / self.n_frames)

        return self.results.cdf

    @property
    def edges(self):
        wmsg = ("The `edges` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.bins` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.edges

    @property
    def count(self):
        wmsg = ("The `count` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.bins` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.count

    @property
    def bins(self):
        wmsg = ("The `bins` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.bins` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.bins

    @property
    def rdf(self):
        wmsg = ("The `rdf` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.rdf` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.rdf

    @property
    def cdf(self):
        wmsg = ("The `cdf` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.cdf` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.cdf
