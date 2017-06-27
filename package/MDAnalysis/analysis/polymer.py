# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
Polymer analysis --- :mod:`MDAnalysis.analysis.polymer`
=======================================================


:Author: Richard J. Gowers
:Year: 2015
:Copyright: GNU Public License v3

This module contains various commonly used tools in analysing polymers.

"""
from __future__ import division, absolute_import
from six.moves import range

import numpy as np
import scipy.optimize

import logging

from .. import NoDataError
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
    def __init__(self, atomgroups, **kwargs):
        """Calculate the persistence length for polymer chains

        Parameters
        ----------
        atomgroups : list
            List of atomgroups.  Each atomgroup should represent a single
            polymer chain, ordered in the correct order.
        start : int, optional
            First frame of trajectory to analyse, Default: None becomes 0.
        stop : int, optional
            Last frame of trajectory to analyse, Default: None becomes
            n_frames.
        step : int, optional
            Frame index to stop analysis. Default: None becomes
            n_frames. Iteration stops *before* this frame number.
        """
        super(PersistenceLength, self).__init__(
            atomgroups[0].universe.trajectory, **kwargs)
        self._atomgroups = atomgroups

        # Check that all chains are the same length
        lens = [len(ag) for ag in atomgroups]
        chainlength = len(atomgroups[0])
        if not all(l == chainlength for l in lens):
            raise ValueError("Not all AtomGroups were the same size")

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
        norm *= len(self._atomgroups) * self.n_frames

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

        try:
            self.results
        except AttributeError:
            raise NoDataError("Use the run method first")
        self.x = np.arange(len(self.results)) * self.lb

        self.lp = fit_exponential_decay(self.x, self.results)

        self.fit = np.exp(-self.x/self.lp)

    def plot(self, ax=None):
        """Oooh fancy"""
        import matplotlib.pyplot as plt
        if ax is None:
            ax = plt.gca()
        ax.plot(self.x, self.results, 'ro', label='Result')
        ax.plot(self.x, self.fit, label='Fit')
        ax.set(xlabel='x', ylabel='C(x)', xlim=[0.0, 40 * self.lb])
        ax.legend(loc='best')
        return ax


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

    """
    def expfunc(x, a):
        return np.exp(-x/a)

    a = scipy.optimize.curve_fit(expfunc, x, y)[0][0]

    return a
