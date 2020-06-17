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

.. math::

   N_{ab}(r) = \rho G_{ab}(r)

(with the appropriate density :math:`\rho`). The latter function can be used to
compute, for instance, coordination numbers such as the number of neighbors in
the first solvation shell :math:`N(r_1)` where :math:`r_1` is the position of
the first minimum in :math:`g(r)`.


.. _`pair distribution functions`:
   https://en.wikipedia.org/wiki/Pair_distribution_function
.. _`radial distribution functions`:
   https://en.wikipedia.org/wiki/Radial_distribution_function


Average radial distribution function
------------------------------------

:class:`InterRDF` is a tool to calculate average radial distribution functions
between two groups of atoms. Suppose we have two AtomGroups ``A`` and
``B``. ``A`` contains atom ``A1``, ``A2``, and ``B`` contains ``B1``,
``B2``. Give ``A`` and ``B`` to class:`InterRDF`, the output will be the
average of RDFs between ``A1`` and ``B1``, ``A1`` and ``B2``, ``A2`` and
``B1``, ``A2`` and ``B2``. A typical application is to calculate the RDF of
solvent with itself or with another solute.

.. autoclass:: InterRDF
   :members:
   :inherited-members:


Site-specific radial distribution function
------------------------------------------

:class:`InterRDF_s` calculates site-specific radial distribution
functions. Instead of two groups of atoms it takes as input a list of pairs of
AtomGroup, ``[[A, B], [C, D], ...]``. Give the same ``A`` and ``B`` to
:class:`InterRDF_s`, the output will be a list of RDFs between ``A1`` and
``B1``, ``A1`` and ``B2``, ``A2`` and ``B1``, ``A2`` and ``B2`` (and similarly
for ``C`` and ``D``). These site-specific radial distribution functions are
typically calculated if one is interested in the solvation shells of a ligand
in a binding site or the solvation of specific residues in a protein. A common
use case is to choose ``A`` and ``C`` to be AtomGroups that only contain a
single atom and ``W`` all solvent molecules: ``InterRDF_s(u, [[A, W], [B,
W]])`` will then produce the RDF of solvent around atom ``A[0]`` and around
atom ``B[0]``.


.. autoclass:: InterRDF_s
   :members:
   :inherited-members:


.. Not Implemented yet:
.. - Structure factor?
.. - Coordination number

"""
import numpy as np

from ..lib.util import blocks_of
from ..lib import distances
from .base import AnalysisBase


class InterRDF(AnalysisBase):
    """Intermolecular pair distribution function

    InterRDF(g1, g2, nbins=75, range=(0.0, 15.0))

    Arguments
    ---------
    g1 : AtomGroup
      First AtomGroup
    g2 : AtomGroup
      Second AtomGroup
    nbins : int (optional)
          Number of bins in the histogram [75]
    range : tuple or list (optional)
          The size of the RDF [0.0, 15.0]
    exclusion_block : tuple (optional)
          A tuple representing the tile to exclude from the distance
          array. [None]
    verbose : bool (optional)
          Show detailed progress of the calculation if set to ``True``; the
          default is ``False``.


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
    """Site-specific intermolecular pair distribution function

    Arguments
    ---------
    u : Universe
          a Universe that contains atoms in `ags`
    ags : list
          a list of pairs of :class:`~MDAnalysis.core.groups.AtomGroup`
          instances
    nbins : int (optional)
          Number of bins in the histogram [75]
    range : tuple or list (optional)
          The size of the RDF [0.0, 15.0]

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

      plt.plot(rdf.bins, rdf.rdf[0][0][0])

    (Which plots the rdf between the first atom in ``s1`` and the first atom in
    ``s2``)

    To generate the *cumulative distribution function* (cdf), use the
    :meth:`~InterRDF_s.get_cdf` method ::

      cdf = rdf.get_cdf()

    Results are available through the :attr:'cdf' attribute::

      plt.plot(rdf.bins, rdf.cdf[0][0][0])

    (Which plots the cdf between the first atom in ``s1`` and the first atom in
    ``s2``)


    .. versionadded:: 0.19.0

    .. versionchanged:: 1.0.0
       Support for the ``start``, ``stop``, and ``step`` keywords has been
       removed. These should instead be passed to :meth:`InterRDF_s.run`.

    """
    def __init__(self, u, ags,
                 nbins=75, range=(0.0, 15.0), density=True, **kwargs):
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
            nA = len(ag1)
            nB = len(ag2)
            N = nA * nB
            indices.append([ag1.indices, ag2.indices])

            # Average number density
            box_vol = self.volume / self.n_frames
            density = N / box_vol

            if self._density:
                rdf.append(self.count[i] / (density * vol * self.n_frames))
            else:
                rdf.append(self.count[i] / (vol * self.n_frames))

        self.rdf = rdf
        self.indices = indices

    def get_cdf(self):
        """Calculate the cumulative distribution functions (CDF) for all sites.
        Note that this is the actual count within a given radius, i.e.,
        :math:`N(r)`.
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
