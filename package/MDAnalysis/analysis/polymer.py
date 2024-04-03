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
Polymer analysis --- :mod:`MDAnalysis.analysis.polymer`
=======================================================


:Author: Richard J. Gowers
:Year: 2015, 2018
:Copyright: GNU Public License v3

This module contains various commonly used tools in analysing polymers.
"""
import numpy as np
import scipy.optimize
import warnings
import logging

from .. import NoDataError
from ..core.groups import requires, AtomGroup
from ..lib.distances import calc_bonds
from .base import AnalysisBase

logger = logging.getLogger(__name__)


@requires('bonds')
def sort_backbone(backbone):
    """Rearrange a linear AtomGroup into backbone order

    Requires that the backbone has bond information,
    and that only backbone atoms are provided (ie no side
    chains or hydrogens).

    Parameters
    ----------
    backbone : AtomGroup
      the backbone atoms, not necessarily in order

    Returns
    -------
    sorted_backbone : AtomGroup
      backbone in order, so `sorted_backbone[i]` is bonded to
      `sorted_backbone[i - 1]` and `sorted_backbone[i + 1]`


    .. versionadded:: 0.20.0
    """
    if not backbone.n_fragments == 1:
        raise ValueError("{} fragments found in backbone.  "
                         "backbone must be a single contiguous AtomGroup"
                         "".format(backbone.n_fragments))

    branches = [at for at in backbone
                if len(at.bonded_atoms & backbone) > 2]
    if branches:
        # find which atom has too many bonds for easier debug
        raise ValueError(
            "Backbone is not linear.  "
            "The following atoms have more than two bonds in backbone: {}."
            "".format(','.join(str(a) for a in branches)))

    caps = [atom for atom in backbone
           if len(atom.bonded_atoms & backbone) == 1]
    if not caps:
        # cyclical structure
        raise ValueError("Could not find starting point of backbone, "
                         "is the backbone cyclical?")

    # arbitrarily choose one of the capping atoms to be the startpoint
    sorted_backbone = AtomGroup([caps[0]])

    # iterate until the sorted chain length matches the backbone size
    while len(sorted_backbone) < len(backbone):
        # current end of the chain
        end_atom = sorted_backbone[-1]

        # look at all bonded atoms which are also part of the backbone
        # and subtract any that have already been added
        next_atom = (end_atom.bonded_atoms & backbone) - sorted_backbone

        # append this to the sorted backbone
        sorted_backbone += next_atom

    return sorted_backbone


class PersistenceLength(AnalysisBase):
    r"""Calculate the persistence length for polymer chains

    The persistence length is the length at which two points on the polymer
    chain become decorrelated.  This is determined by first measuring the
    autocorrelation (:math:`C(n)`) of two bond vectors
    (:math:`\mathbf{a}_i, \mathbf{a}_{i + n}`) separated by :math:`n` bonds

    .. math::

       C(n) = \langle \cos\theta_{i, i+n} \rangle =
               \langle \mathbf{a_i} \cdot \mathbf{a_{i+n}} \rangle

    An exponential decay is then fitted to this, which yields the
    persistence length

    .. math::

       C(n) \approx \exp\left( - \frac{n \bar{l_B}}{l_P} \right)

    where :math:`\bar{l_B}` is the average bond length, and :math:`l_P` is
    the persistence length which is fitted

    Parameters
    ----------
    atomgroups : iterable
       List of AtomGroups. Each should represent a single
       polymer chain, ordered in the correct order.
    verbose : bool, optional
       Show detailed progress of the calculation if set to ``True``.

    Attributes
    ----------
    results.bond_autocorrelation : numpy.ndarray
       the measured bond autocorrelation
    results.lb : float
       the average bond length

       .. versionadded:: 2.0.0

    lb : float

       Alias to the :attr:`results.lb`.

       .. deprecated:: 2.0.0
            Will be removed in MDAnalysis 3.0.0. Please use
            :attr:`results.lb` instead.

    results.x : numpy.ndarray
        length of the decorrelation predicted by *lp*

        .. versionadded:: 2.0.0

    results.lp : float
       calculated persistence length

       .. versionadded:: 2.0.0

    lp : float

       Alias to the :attr:`results.lp`.

       .. deprecated:: 2.0.0
            Will be removed in MDAnalysis 3.0.0. Please use
            :attr:`results.lp` instead.

    results.fit : numpy.ndarray
       the modelled backbone decorrelation predicted by *lp*

       .. versionadded:: 2.0.0

    fit : float

       Alias to the :attr:`results.fit`.

       .. deprecated:: 2.0.0
            Will be removed in MDAnalysis 3.0.0. Please use
            :attr:`results.fit` instead.

    See Also
    --------
    :func:`sort_backbone`
       for producing the sorted AtomGroup required for input.

    Example
    -------
    .. code-block:: python

        from MDAnalysis.tests.datafiles import TRZ_psf, TRZ
        import MDAnalysis as mda
        from MDAnalysis.analysis import polymer
        u = mda.Universe(TRZ_psf, TRZ)

        # this system is a pure polymer melt of polyamide,
        # so we can select the chains by using the .fragments attribute
        chains = u.atoms.fragments

        # select only the backbone atoms for each chain
        backbones = [chain.select_atoms('not name O* H*') for chain in chains]

        # sort the chains, removing any non-backbone atoms
        sorted_backbones = [polymer.sort_backbone(bb) for bb in backbones]
        persistence_length = polymer.PersistenceLength(sorted_backbones)

        # Run the analysis, this will average over all polymer chains
        # and all timesteps in trajectory
        persistence_length = persistence_length.run()
        print(f'The persistence length is: {persistence_length.results.lp}')

        # always check the visualisation of this:
        persistence_length.plot()


    .. versionadded:: 0.13.0
    .. versionchanged:: 0.20.0
       The run method now automatically performs the exponential fit
    .. versionchanged:: 1.0.0
       Deprecated :meth:`PersistenceLength.perform_fit` has now been removed.
    .. versionchanged:: 2.0.0
       Former ``results`` are now stored as ``results.bond_autocorrelation``.
       :attr:`lb`, :attr:`lp`, :attr:`fit` are now stored in a
       :class:`MDAnalysis.analysis.base.Results` instance.
    """
    def __init__(self, atomgroups, **kwargs):
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
            # Normalized to unit vectors
            vecs /= np.sqrt((vecs * vecs).sum(axis=1))[:, None]

            inner_pr = np.inner(vecs, vecs)
            for i in range(n-1):
                self._results[:(n-1)-i] += inner_pr[i, i:]

    @property
    def lb(self):
        wmsg = ("The `lb` attribute was deprecated in "
                "MDAnalysis 2.0.0 and will be removed in MDAnalysis 3.0.0. "
                "Please use `results.variance` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.lb

    @property
    def lp(self):
        wmsg = ("The `lp` attribute was deprecated in "
                "MDAnalysis 2.0.0 and will be removed in MDAnalysis 3.0.0. "
                "Please use `results.variance` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.lp

    @property
    def fit(self):
        wmsg = ("The `fit` attribute was deprecated in "
                "MDAnalysis 2.0.0 and will be removed in MDAnalysis 3.0.0. "
                "Please use `results.variance` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.fit

    def _conclude(self):
        n = len(self._atomgroups[0])

        norm = np.linspace(n - 1, 1, n - 1)
        norm *= len(self._atomgroups) * self.n_frames

        self.results.bond_autocorrelation = self._results / norm
        self._calc_bond_length()

        self._perform_fit()

    def _calc_bond_length(self):
        """calculate average bond length"""
        bs = []
        for ag in self._atomgroups:
            pos = ag.positions
            b = calc_bonds(pos[:-1], pos[1:]).mean()
            bs.append(b)
        self.results.lb = np.mean(bs)

    def _perform_fit(self):
        """Fit the results to an exponential decay"""
        try:
            self.results.bond_autocorrelation
        except AttributeError:
            raise NoDataError("Use the run method first") from None
        self.results.x = self.results.lb *\
                            np.arange(len(self.results.bond_autocorrelation))

        self.results.lp = fit_exponential_decay(self.results.x,
                                                self.results.bond_autocorrelation)

        self.results.fit = np.exp(-self.results.x/self.results.lp)

    def plot(self, ax=None):
        """Visualize the results and fit

        Parameters
        ----------
        ax : matplotlib.Axes, optional
          if provided, the graph is plotted on this axis

        Returns
        -------
        ax : the axis that the graph was plotted on
        """
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(self.results.x,
                self.results.bond_autocorrelation,
                'ro',
                label='Result')
        ax.plot(self.results.x,
                self.results.fit,
                label='Fit')
        ax.set_xlabel(r'x')
        ax.set_ylabel(r'$C(x)$')
        ax.set_xlim(0.0, 40 * self.results.lb)

        ax.legend(loc='best')

        return ax


def fit_exponential_decay(x, y):
    r"""Fit a function to an exponential decay

    .. math::  y = \exp\left(- \frac{x}{a}\right)

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
