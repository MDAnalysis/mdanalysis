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

Tools for calculating pair distribution functions ("radial
distribution functions" or "RDF").

.. Not Implemented yet:
.. - Structure factor?
.. - Coordination number

"""
from __future__ import division, absolute_import
import numpy as np

from ..lib.util import blocks_of
from ..lib import distances
from .base import AnalysisBase


class InterRDF_s(AnalysisBase):
    """Intermolecular pair distribution function

    InterRDF_s(u, g1, g2, nbins=75, range=(0.0, 15.0))

    Arguments
    ---------
    u : AtomGroup
       An AtomGroup contains atoms in g1 and g2
    g1 : list
        A list of first AtomGroup selection strings
    g2 : list
        A list of second AtomGroup selection strings
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

    Example
    -------
    First create the :class:`InterRDF` object, by supplying two
    AtomGroups then use the :meth:`run` method ::

      rdf = InterRDF(u, ag1, ag2)
      rdf.run()

    Results are available through the :attr:`bins` and :attr:`rdf`
    attributes::

      plt.plot(rdf.bins, rdf.rdf)

    The `exclusion_block` keyword allows the masking of pairs from
    within the same molecule.  For example, if there are 7 of each
    atom in each molecule, the exclusion mask `(7, 7)` can be used.

    .. versionadded:: 0.13.0

    """
    def __init__(self, u, g1, g2,
                 nbins=75, range=(0.0, 15.0), density=True, exclusion_block=None,
                 **kwargs):
        super(InterRDF_s, self).__init__(u.universe.trajectory, **kwargs)
        self.g1 = []
        self.g2 = []
        for sel1 in g1:
            self.g1.append(u.select_atoms(sel1))
        for sel2 in g2:
            self.g2.append(u.select_atoms(sel2))
        self.u = u.universe
        self._density = density

        self.rdf_settings = {'bins': nbins,
                             'range': range}
        self._exclusion_block = exclusion_block

    def _prepare(self):
        # Empty list to store the RDF
        count_list = []
        count, edges = np.histogram([-1], **self.rdf_settings)
        for n  in range(len(self.g1)):
            count_list.append(np.zeros((len(self.g1[n]), len(self.g2[n]),
                              len(count)), dtype=np.float64))
        self.count = count_list
        self.edges = edges
        self.bins = 0.5 * (edges[:-1] + edges[1:])

        # Need to know average volume
        self.volume = 0.0

        # Allocate a results array which we will reuse
        self._result = np.zeros((len(self.g1), len(self.g2)), dtype=np.float64)

        # exclusoin
        # If provided exclusions, create a mask of _result which
        # lets us take these out
        #if self._exclusion_block is not None:
        #    self._exclusion_mask = blocks_of(self._result,
        #                                     *self._exclusion_block)
        #    self._maxrange = self.rdf_settings['range'][1] + 1.0
        #else:
        #    self._exclusion_mask = None

    def _single_frame(self):
        for i in range(len(self.g1)):
            self._result=distances.distance_array(self.g1[i].positions, self.g2[i].positions,
                                                  box=self.u.dimensions)

            # Maybe exclude same molecule distances
            #if self._exclusion_mask is not None:
            #   self._exclusion_mask[:] = self._maxrange


            for j in range(len(self.g1[i])):
                for k in range(len(self.g2[i])):
                    count = np.histogram(self._result[j,k], **self.rdf_settings)[0]
                    self.count[i][j,k,:] += count

        self.volume += self._ts.volume


    def _conclude(self):
        # Volume in each radial shell
        vol = np.power(self.edges[1:], 3) - np.power(self.edges[:-1], 3)
        vol *= 4/3.0 * np.pi

        # Empty lists to restore index, RDF and CDF
        index = []
        rdf_s = []
        cdf_s = []

        for i in range(len(self.g1)):

            cdf_s.append(np.zeros(np.shape(self.count[i])))
            for j in range(len(self.g1[i])):
                for k in range(len(self.g2[i])):
                    cdf_s[i][j,k,:] = np.array([sum(self.count[i][j,k,0:(n+1)])
                                               for n in range(len(self.count[i][j,k,:]))]) / self.n_frames

            # Number of each selection
            nA = len(self.g1[i])
            nB = len(self.g2[i])
            N = nA * nB
            index.append([self.g1[i].indices, self.g2[i].indices])

            # If we had exclusions, take these into account
            #if self._exclusion_block:
            #    xA, xB = self._exclusion_block
            #    nblocks = nA / xA
            #    N -= xA * xB * nblocks

            # Average number density
            box_vol = self.volume / self.n_frames
            density = N / box_vol

            if self._density:
                rdf_s.append(self.count[i] / (density * vol * self.n_frames))
            else:
                rdf_s.append(self.count[i] / (vol * self.n_frames))

        self.cdf_s = cdf_s
        self.rdf_s = rdf_s
        self.index = index
