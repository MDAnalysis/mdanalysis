# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
Radial Distribution Functions --- :mod:`MDAnalysis.analysis.rdf`
================================================================

Tools for calculating pair distribution functions

# TODO
 - Structure factor?
 - Coordination number
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
    g1
      First AtomGroup
    g2
      Second AtomGroup

    Keywords
    --------
    nbins
          Number of bins in the histogram [75]
    range
          The size of the RDF [0.0, 15.0]
    exclusion_block
          A tuple representing the tile to exclude from the distance
          array. [None]
    start
          The frame to start at [0]
    stop
          The frame to end at [-1]
    step
          The step size through the trajectory in frames [0]

    Example
    -------
    First create the InterRDF object, by supplying two AtomGroups
    then use the `run` method

      rdf = InterRDF(ag1, ag2)
      rdf.run()

    Results are available through the .bins and .rdf attributes

      plt.plot(rdf.bins, rdf.rdf)

    The `exclusion_block` keyword allows the masking of pairs from
    within the same molecule.  For example, if there are 7 of each
    atom in each molecule, the exclusion mask (7, 7) can be used.

    .. versionadded:: 0.13.0
    """
    def __init__(self, g1, g2,
                 nbins=75, range=(0.0, 15.0), exclusion_block=None,
                 start=None, stop=None, step=None):
        self.g1 = g1
        self.g2 = g2
        self.u = g1.universe

        self._setup_frames(self.u.trajectory,
                           start=start,
                           stop=stop,
                           step=step)

        self.rdf_settings = {'bins':nbins,
                             'range':range}

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
        if exclusion_block is not None:
            self._exclusion_block = exclusion_block
            self._exclusion_mask = blocks_of(self._result, *exclusion_block)
            self._maxrange = range[1] + 1.0
        else:
            self._exclusion_block = None
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
        box_vol = self.volume / self.nframes
        density = N / box_vol

        rdf = self.count / (density * vol * self.nframes)

        self.rdf = rdf
