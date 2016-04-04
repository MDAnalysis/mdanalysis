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
Polymer analysis --- :mod:`MDAnalysis.analysis.polymer`
=======================================================


:Author: Richard J. Gowers
:Year: 2015
:Copyright: GNU Public License v3

This module contains various commonly used tools in analysing polymers.

"""
from six.moves import range

import numpy as np
import logging

from .. import NoDataError
from ..lib.util import blocks_of
from ..lib.distances import calc_bonds
from .base import AnalysisBase

logger = logging.getLogger(__name__)


class PersistenceLength(AnalysisBase):
    r"""Calculate the persistence length for polymer chains

    The persistence length is the length at which two points on the polymer
    chain become decorrelated.

    Notes
    -----
    This analysis requires that the trajectory supports indexing

    .. versionadded:: 0.13.0
    """
    def __init__(self, atomgroups,
                 start=None, stop=None, step=None):
        """Calculate the persistence length for polymer chains

        Parameters
        ----------
        atomgroups : list
            List of atomgroups.  Each atomgroup should represent a single
            polymer chain, ordered in the correct order.
        start : int, optional
            First frame of trajectory to analyse, Default: 0
        stop : int, optional
            Last frame of trajectory to analyse, Default: -1
        step : int, optional
            Step between frames to analyse, Default: 1
        """
        self._atomgroups = atomgroups

        # Check that all chains are the same length
        lens = [len(ag) for ag in atomgroups]
        chainlength = len(atomgroups[0])
        if not all( l == chainlength for l in lens):
            raise ValueError("Not all AtomGroups were the same size")

        self._setup_frames(atomgroups[0].universe.trajectory,
                           start, stop, step)

        self._results = np.zeros(chainlength - 1, dtype=np.float32)

    def _single_frame(self):
        # could optimise this by writing a "self dot array"
        # we're only using the upper triangle of np.inner
        # function would accept a bunch of coordinates and spit out the
        # decorrel for that 
        n = len(self._atomgroups[0])

        for chain in self._atomgroups:
            # Vector from each atom to next
            vecs = chain.positions[1:] - chain.positions[:-1]
            # Normalised to unit vectors
            vecs /= np.sqrt((vecs * vecs).sum(axis=1))[:, None]

            inner_pr = np.inner(vecs, vecs)
            for i in range(n-1):
                self._results[:(n-1)-i] += inner_pr[i, i:]

    def _conclude(self):
        n = len(self._atomgroups[0])

        norm = np.linspace(n - 1, 1, n - 1)
        norm *= len(self._atomgroups) * self.nframes

        self.results = self._results / norm
        self._calc_bond_length()

    def _calc_bond_length(self):
        """calculate average bond length"""
        bs = []
        for ag in self._atomgroups:
            pos = ag.positions
            b = calc_bonds(pos[:-1], pos[1:]).mean()
            bs.append(b)
        self.lb = np.mean(bs)

    def perform_fit(self):
        """Fit the results to an exponential decay"""
        from scipy.optimize import curve_fit

        try:
            results = self.results
        except AttributeError:
            raise NoDataError("Use the run method first")
        self.x = np.arange(len(self.results)) * self.lb

        self.lp = fit_exponential_decay(self.x, self.results)

        self.fit = np.exp(-self.x/self.lp)

    def plot(self):
        """Oooh fancy"""
        import matplotlib.pyplot as plt
        plt.ylabel('C(x)')
        plt.xlabel('x')
        plt.xlim([0.0, 40 * self.lb])
        plt.plot(self.x, self.results, 'ro')
        plt.plot(self.x, self.fit)
        plt.show()


def fit_exponential_decay(x, y):
    r"""Fit a function to an exponential decay

    .. math::  y = \exp(-x/a)

    Parameters
    ----------
    x, y : array_like
      The two arrays of data

    Returns
    -------
    a : float
      The coefficient *a* for this decay

    Notes
    -----
    This function assumes that data starts at 1.0 and decays to 0.0

    Requires scipy
    """
    from scipy.optimize import curve_fit

    def expfunc(x, a):
        return np.exp(-x/a)

    a = curve_fit(expfunc, x, y)[0][0]

    return a

    
