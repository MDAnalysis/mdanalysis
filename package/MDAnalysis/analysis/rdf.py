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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Radial Distribution Functions --- :mod:`MDAnalysis.analysis.rdf`
================================================================

This module contains two classes to calculate pair distribution functions
 ("radial distribution functions" or "RDF").

:class:`InterRDF` is a tool to calculate average radial distribution functions
between two groups of atoms. Suppose we have two atom group A and B. A contains
atom A1, A2, and B contains B1, B2. Give A and B to class:`InterRDF`, the output
will be the average of RDFs bewteen A1 and B1, A1 and B2, A2 and B1, A2 and B2.

:class:`InterRDF_s` is a tool to calculate site-specific radial distribution
functions. Give the same A and B to class:`InterRDF_s`, the output will be
a list of RDFs between A1 and B1, A1 and B2, A2 and B1, A2 and B2, which are the
site-specific radial distribution functions.


Classes:
-------

.. autoclass:: InterRDF
.. autoclass:: InterRDF_s

-------

.. Not Implemented yet:
.. - Structure factor?
.. - Coordination number

"""
from __future__ import division, absolute_import
import numpy as np

from ..lib.util import blocks_of
from ..lib import distances
from .base import AnalysisBase
from six.moves import zip, range


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
    start : int (optional)
          The frame to start at (default is first)
    stop : int (optional)
          The frame to end at (default is last)
    step : int (optional)
          The step size through the trajectory in frames (default is
          every frame)
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

        # Allocate a results array which we will reuse
        self._result = np.zeros((len(self.g1), len(self.g2)), dtype=np.float64)
        # If provided exclusions, create a mask of _result which
        # lets us take these out
        if self._exclusion_block is not None:
            self._exclusion_mask = blocks_of(self._result,
                                             *self._exclusion_block)
            self._maxrange = self.rdf_settings['range'][1] + 1.0
        else:
            self._exclusion_mask = None

    def _single_frame(self):
        distances.distance_array(self.g1.positions, self.g2.positions,
                                 box=self.u.dimensions, result=self._result)
        # Maybe exclude same molecule distances
        if self._exclusion_mask is not None:
            self._exclusion_mask[:] = self._maxrange

        count = np.histogram(self._result, **self.rdf_settings)[0]
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



"""
Site-specific Radial Distribution Functions --- :mod:`MDAnalysis.analysis.rdf_s`
================================================================

Tools for calculating site-specific pair distribution functions ("radial
distribution functions" or "RDF").

"""

class InterRDF_s(AnalysisBase):
    """Site-specific intermolecular pair distribution function

    Arguments
    ---------
    u : Universe
       A Universe contains atoms in ags
    ags : list
         A list of pairs of AtomGroups
    nbins : int (optional)
          Number of bins in the histogram [75]
    range : tuple or list (optional)
          The size of the RDF [0.0, 15.0]
    start : int (optional)
          The frame to start at (default is first)
    stop : int (optional)
          The frame to end at (default is last)
    step : int (optional)
          The step size through the trajectory in frames (default is
          every frame)

    Example
    -------
    First create the :class:`InterRDF_s` object, by supplying one Universe
    and one list of pairs of AtomGroups then use the :meth:`run` method ::

      from MDAnalysisTests.datafiles import GRO_MEMPROT, XTC_MEMPROT
      u = mda.Universe(GRO_MEMPROT, XTC_MEMPROT)

      s1 = u.select_atoms('name ZND and resid 289')
      s2 = u.select_atoms('(name OD1 or name OD2) and resid 51 and sphzone 5.0 (resid 289)')
      s3 = u.select_atoms('name ZND and (resid 291 or resid 292)')
      s4 = u.select_atoms('(name OD1 or name OD2) and sphzone 5.0 (resid 291)')
      ags = [[s1, s2], [s3, s4]]

      rdf = InterRDF_s(u, ags)
      rdf.run()

    Results are available through the :attr:`bins` and :attr:`rdf`
    attributes::

      plt.plot(rdf.bins, rdf.rdf[0][0][0])

    (Which plots the rdf between the first atom in s1 and the first atom in s2)

    To generate cdf, use the 'cdf' method

      cdf = rdf.get_cdf()

    Results are available through the :attr:'cdf' attributes::

      plt.plot(rdf.bins, rdf.cdf[0][0][0])

    (Which plots the cdf between the first atom in s1 and the first atom in s2)

    .. versionadded:: 0.19.0

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


    def _single_frame(self):
        for i, (ag1, ag2) in enumerate(self.ags):
            result=distances.distance_array(ag1.positions, ag2.positions,
                                            box=self.u.dimensions)
            for j in range(ag1.n_atoms):
                for k in range(ag2.n_atoms):
                    count = np.histogram(result[j, k], **self.rdf_settings)[0]
                    self.count[i][j, k, :] += count

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
