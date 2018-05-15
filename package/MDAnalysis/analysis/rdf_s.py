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
Site-specific Radial Distribution Functions --- :mod:`MDAnalysis.analysis.rdf_s`
================================================================

Tools for calculating site-specific pair distribution functions ("radial
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
from six.moves import zip, range


class InterRDF_s(AnalysisBase):
    """Intermolecular pair distribution function

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

    Results are available through the :attr:`bins` and :attr:`rdf_s`
    attributes::

      plt.plot(rdf.bins, rdf.rdf_s[0])

    To generate cdf, use the 'cdf' mehthod

      cdf = rdf.get_cdf()

    Results are available through the :attr:'cdf' attributes::

      plt.plot(rdf.bins, rdf.cdf[0])

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
        rdf_s = []

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
                rdf_s.append(self.count[i] / (density * vol * self.n_frames))
            else:
                rdf_s.append(self.count[i] / (vol * self.n_frames))

        self.rdf_s = rdf_s
        self.indices = indices

    def get_cdf(self):
        # Calculate cumulative distribution function
        # Empty list to restore CDF
        cdf = []

        for count in self.count:
            cdf.append(np.cumsum(count, axis=2) / self.n_frames)

        # Results stored in self.cdf
        # self.cdf is a list of cdf between pairs of AtomGroups in ags
        self.cdf = cdf
