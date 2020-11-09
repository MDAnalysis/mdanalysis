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
================================================================

This module contains two classes to calculate radial
`pair distribution functions`_ (`radial distribution functions`_ or "RDF").
The RDF :math:`g_{ab}(r)` between types of particles :math:`a` and :math:`b` is

.. _equation-gab:

.. math::

   g_{ab}(r) = (N_{a} N_{b})^{-1} \sum_{i=1}^{N_a} \sum_{j=1}^{N_b}
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

In :class:`InterRDF_s`, we provide an option `density`. When `density` is
``False``, it will return the RDF :math:`g_{ab}(r)`; when `density` is
``True``, it will return the density of particle :math:`b` in a shell at
distance :math:`r` around a :math:`a` particle, which is

.. _equation-nab:

.. math::
   n_{ab}(r) = \rho g_{ab}(r)

.. _`pair distribution functions`:
   https://en.wikipedia.org/wiki/Pair_distribution_function
.. _`radial distribution functions`:
   https://en.wikipedia.org/wiki/Radial_distribution_function


Average radial distribution function
------------------------------------

:class:`InterRDF` is a tool to calculate average radial distribution functions
between two groups of atoms. Suppose we have two AtomGroups ``A`` and
``B``. ``A`` contains atom ``A1``, ``A2``, and ``B`` contains ``B1``,
``B2``. Given ``A`` and ``B`` to :class:`InterRDF`, the output will be the
average of RDFs between ``A1`` and ``B1``, ``A1`` and ``B2``, ``A2`` and
``B1``, ``A2`` and ``B2``. A typical application is to calculate the RDF of
solvent with itself or with another solute.

.. autoclass:: InterRDF
   :members:
   :inherited-members:

   .. attribute:: bins

      :class:`numpy.ndarray` of the centers of the `nbins` histogram
      bins.

   .. attribute:: edges

      :class:`numpy.ndarray` of the `nbins + 1` edges of the histogram
      bins.

   .. attribute:: rdf

      :class:`numpy.ndarray` of the :ref:`radial distribution
      function<equation-gab>` values for the :attr:`bins`.

   .. attribute:: count

      :class:`numpy.ndarray` representing the radial histogram, i.e.,
      the raw counts, for all :attr:`bins`.


Site-specific radial distribution function
------------------------------------------

:class:`InterRDF_s` calculates site-specific radial distribution
functions. Instead of two groups of atoms it takes as input a list of pairs of
AtomGroup, ``[[A, B], [C, D], ...]``. Given the same ``A`` and ``B`` to
:class:`InterRDF_s`, the output will be a list of individual RDFs between
``A1`` and ``B1``, ``A1`` and ``B2``, ``A2`` and ``B1``, ``A2`` and ``B2`` (and
similarly for ``C`` and ``D``). These site-specific radial distribution
functions are typically calculated if one is interested in the solvation shells
of a ligand in a binding site or the solvation of specific residues in a
protein.

.. autoclass:: InterRDF_s
   :members:
   :inherited-members:

   .. attribute:: bins

      :class:`numpy.ndarray` of the centers of the `nbins` histogram
      bins; all individual site-specific RDFs have the same bins.

   .. attribute:: edges

      :class:`numpy.ndarray` of the `nbins + 1` edges of the histogram
      bins; all individual site-specific RDFs have the same bins.

   .. attribute:: rdf

      :class:`list` of the site-specific :ref:`radial distribution
      functions<equation-gab>` or :ref:`density
      functions<equation-nab>` for the :attr:`bins`. The list contains
      ``len(ags)`` entries. Each entry for the ``i``-th pair ``[A, B]
      = ags[i]`` in `ags` is a :class:`numpy.ndarray` with shape
      ``(len(A), len(B))``, i.e., a stack of RDFs. For example,
      ``rdf[i][0, 2]`` is the RDF between atoms ``A[0]`` and ``B[2]``.

   .. attribute:: count

      :class:`list` of the site-specific radial histograms, i.e., the
      raw counts, for all :attr:`bins`. The data have the same
      structure as :attr:`rdf` except that the arrays contain the raw
      counts.

   .. attribute:: cdf

      :class:`list` of the site-specific :ref:`cumulative
      counts<equation-countab>`, for all :attr:`bins`. The data have the same
      structure as :attr:`rdf` except that the arrays contain the cumulative
      counts.

      This attribute only exists after :meth:`get_cdf` has been run.




.. Not Implemented yet:
.. - Structure factor?
.. - Coordination number

"""
import numpy as np

from ..lib.util import blocks_of
from ..lib import distances
from .base import AnalysisBase


class InterRDF(AnalysisBase):
    r"""Intermolecular pair distribution function

    The :ref:`radial distribution function<equation-gab>` is calculated by
    histogramming distances between all particles in `g1` and `g2` while taking
    periodic boundary conditions into account via the minimum image
    convention.

    The `exclusion_block` keyword may be used to exclude a set of distances
    from the calculations.

    Results are available in the attributes :attr:`rdf` and :attr:`count`.

    Arguments
    ---------
    g1 : AtomGroup
      First AtomGroup
    g2 : AtomGroup
      Second AtomGroup
    nbins : int (optional)
          Number of bins in the histogram
    range : tuple or list (optional)
          The size of the RDF
    exclusion_block : tuple (optional)
          A tuple representing the tile to exclude from the distance
          array.
    verbose : bool (optional)
          Show detailed progress of the calculation if set to ``True``

    Example
    -------
    First create the :class:`InterRDF` object, by supplying two
    AtomGroups then use the :meth:`run` method ::

      rdf = InterRDF(ag1, ag2)
      rdf.run()

    Results are available through the :attr:`bins` and :attr:`rdf`
    attributes::

      plt.plot(rdf.bins, rdf.rdf)

    The `exclusion_block` keyword allows the masking of pairs from
    within the same molecule.  For example, if there are 7 of each
    atom in each molecule, the exclusion mask `(7, 7)` can be used.


    .. versionadded:: 0.13.0

    .. versionchanged:: 1.0.0
       Support for the ``start``, ``stop``, and ``step`` keywords has been
       removed. These should instead be passed to :meth:`InterRDF.run`.

    """
    def __init__(self, g1, g2,
                 nbins=75, range=(0.0, 15.0), exclusion_block=None,
                 **kwargs):
        super(InterRDF, self).__init__(g1.universe.trajectory, **kwargs)
        self.g1 = g1
        self.g2 = g2
        self.u = g1.universe

        self.rdf_settings = {'bins': nbins,
                             'range': range}
        self._exclusion_block = exclusion_block

    def _prepare(self):
        # Empty histogram to store the RDF
        count, edges = np.histogram([-1], **self.rdf_settings)
        count = count.astype(np.float64)
        count *= 0.0
        self.count = count
        self.edges = edges
        self.bins = 0.5 * (edges[:-1] + edges[1:])

        # Need to know average volume
        self.volume = 0.0
        # Set the max range to filter the search radius
        self._maxrange = self.rdf_settings['range'][1]


    def _single_frame(self):
        pairs, dist = distances.capped_distance(self.g1.positions,
                                                self.g2.positions,
                                                self._maxrange,
                                                box=self.u.dimensions)
        # Maybe exclude same molecule distances
        if self._exclusion_block is not None:
            idxA, idxB = pairs[:, 0]//self._exclusion_block[0], pairs[:, 1]//self._exclusion_block[1]
            mask = np.where(idxA != idxB)[0]
            dist = dist[mask]


        count = np.histogram(dist, **self.rdf_settings)[0]
        self.count += count

        self.volume += self._ts.volume

    def _conclude(self):
        # Number of each selection
        nA = len(self.g1)
        nB = len(self.g2)
        N = nA * nB

        # If we had exclusions, take these into account
        if self._exclusion_block:
            xA, xB = self._exclusion_block
            nblocks = nA / xA
            N -= xA * xB * nblocks

        # Volume in each radial shell
        vol = np.power(self.edges[1:], 3) - np.power(self.edges[:-1], 3)
        vol *= 4/3.0 * np.pi

        # Average number density
        box_vol = self.volume / self.n_frames
        density = N / box_vol

        rdf = self.count / (density * vol * self.n_frames)

        self.rdf = rdf


class InterRDF_s(AnalysisBase):
    r"""Site-specific intermolecular pair distribution function

    Arguments
    ---------
    u : Universe
          a Universe that contains atoms in `ags`
    ags : list
          a list of pairs of :class:`~MDAnalysis.core.groups.AtomGroup`
          instances
    nbins : int (optional)
          Number of bins in the histogram
    range : tuple or list (optional)
          The size of the RDF
    density : bool (optional)
          ``False``: calculate :math:`g_{ab}(r)`; ``True``: calculate
          the true :ref:`single particle density<equation-nab>`
          :math:`n_{ab}(r)`.

          .. versionadded:: 1.0.1

             This keyword was available since 0.19.0 but was not
             documented. Furthermore, it had the opposite
             meaning. Since 1.0.1 it is officially supported as
             documented.


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

    Results are available through the :attr:`bins` and :attr:`rdf` attributes::

      plt.plot(rdf.bins, rdf.rdf[0][0, 0])

    (Which plots the rdf between the first atom in ``s1`` and the first atom in
    ``s2``)

    To generate the *cumulative distribution function* (cdf) in the sense of
    "particles within radius :math:`r`", i.e., :math:`N_{ab}(r)`, use the
    :meth:`~InterRDF_s.get_cdf` method ::

      cdf = rdf.get_cdf()

    Results are available through the :attr:`cdf` attribute::

      plt.plot(rdf.bins, rdf.cdf[0][0, 0])

    (Which plots the cdf between the first atom in ``s1`` and the first atom in
    ``s2``)


    .. versionadded:: 0.19.0

    .. versionchanged:: 1.0.0
       Support for the ``start``, ``stop``, and ``step`` keywords has been
       removed. These should instead be passed to :meth:`InterRDF_s.run`.

    """
    def __init__(self, u, ags,
                 nbins=75, range=(0.0, 15.0), density=False, **kwargs):
        super(InterRDF_s, self).__init__(u.universe.trajectory, **kwargs)

        # List of pairs of AtomGroups
        self.ags = ags
        self.u = u
        self._density = density
        self.rdf_settings = {'bins': nbins,
                             'range': range}

    def _prepare(self):
        # Empty list to store the RDF
        count_list = []
        count, edges = np.histogram([-1], **self.rdf_settings)
        count_list = [np.zeros((ag1.n_atoms, ag2.n_atoms, len(count)), dtype=np.float64)
                         for ag1, ag2 in self.ags]

        self.count = count_list
        self.edges = edges
        self.bins = 0.5 * (edges[:-1] + edges[1:])

        # Need to know average volume
        self.volume = 0.0
        self._maxrange = self.rdf_settings['range'][1]


    def _single_frame(self):
        for i, (ag1, ag2) in enumerate(self.ags):
            pairs, dist = distances.capped_distance(ag1.positions,
                                                    ag2.positions,
                                                    self._maxrange,
                                                    box=self.u.dimensions)

            for j, (idx1, idx2) in enumerate(pairs):
                self.count[i][idx1, idx2, :] += np.histogram(dist[j],
                                                             **self.rdf_settings)[0]

        self.volume += self._ts.volume


    def _conclude(self):
        # Volume in each radial shell
        vol = np.power(self.edges[1:], 3) - np.power(self.edges[:-1], 3)
        vol *= 4/3.0 * np.pi

        # Empty lists to restore indices, RDF
        indices = []
        rdf = []

        for i, (ag1, ag2) in enumerate(self.ags):
            # Number of each selection
            indices.append([ag1.indices, ag2.indices])

            # Average number density
            box_vol = self.volume / self.n_frames
            density = 1 / box_vol

            if self._density:
                rdf.append(self.count[i] / (vol * self.n_frames))
            else:
                rdf.append(self.count[i] / (density * vol * self.n_frames))

        self.rdf = rdf
        self.indices = indices

    def get_cdf(self):
        r"""Calculate the cumulative counts for all sites.

        This is the :ref:`cumulative count<equation-countab>` within a given
        radius, i.e., :math:`N_{ab}(r)`.

        The result is returned and also stored in the attribute
        :attr:`cdf`.


        Returns
        -------
        cdf : list
              list of arrays with the same structure as :attr:`rdf`

        """
        # Calculate cumulative distribution function
        # Empty list to restore CDF
        cdf = []

        for count in self.count:
            cdf.append(np.cumsum(count, axis=2) / self.n_frames)

        # Results stored in self.cdf
        # self.cdf is a list of cdf between pairs of AtomGroups in ags
        self.cdf = cdf

        return cdf
