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
Principal Component Analysis (PCA) --- :mod:`MDAnalysis.analysis.pca`
=====================================================================

:Authors: John Detlefs
:Year: 2016
:Copyright: GNU Public License v3

This module contains the linear dimension reduction method Principal
Component Analysis. This module constructs a covariance matrix wherein each
element of the matrix is denoted by (i,j) row-column coordinates. The (i,j)
coordinate reflects the influence of the of the ith frame's coordinates on the
jth frame's coordinates of a given trajectory. The eponymous components are the
eigenvectors of this matrix.

For each eigenvector, its eigenvalue reflects the variance that the eigenvector
explains. This value is made into a ratio. Stored in
:attribute:`explained_variance`, this ratio divides the accumulated variance
of the nth eigenvector and the n-1 eigenvectors preceding the eigenvector by
the total variance in the data. For most data, :attribute:`explained_variance`
will be approximately equal to one for some n that is significantly smaller
than the total number of components, these are the components of interest given
by Principal Component Analysis.

From here, we can project a trajectory onto these principal components and
attempt to retrieve some structure from our high dimensional data. We have
provided a [notebook](# TODO edit max's notebook to use the new module)
containing a thorough demonstration of Principal Component Analysis.

For a basic introduction to the module, the :ref:`PCA-tutorial` shows how
to perform Principal Component Analysis.


.. _PCA-tutorial:
The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF` and
:data:`~MDAnalysis.tests.datafiles.DCD`). This tutorial shows how to use the
PCA class.


First load all modules and test data ::
    >>> import MDAnalysis as mda
    >>> import numpy as np
    >>> import MDAnalysis.analysis.pca as pca
    >>> from MDAnalysis.tests.datafiles import PSF, DCD

Given a universe containing trajectory data we can perform PCA using
:class:`PCA`:: and retrieve the principal components.
    >>> u = mda.Universe(PSF,DCD)
    >>> PSF_pca = pca.PCA(u)
    >>> cumulated_variance, principal_components = PSF_pca.fit()

Inspect the components to determine the principal components you would like
to retain. The choice is arbitrary, but I will stop when 95 percent of the
variance is explained by the components.
    >>> n_pcs = next(x[0] for x in enumerate(cumulated_variance) if x[1] > 0.95)
    >>> pca_space = PSF_pca.transform(n_components=n_pcs)

From here, inspection of the pca_space and conclusions to be drawn from the
data are left to the user.
"""
from six.moves import range
import logging

import numpy as np
from MDAnalysis import Universe

from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.pca")

class PCA(AnalysisBase):
    """Principal component analysis on an MD trajectory

    Attributes
    ----------
    p_components: array, (n_components, n_atoms)
        The principal components of the feature space,
        representing the directions of maximum variance in the data.
    variance : array (n_components, )
        The raw variance explained by each eigenvector of the covariance
        matrix.
    explained_variance : array, (n_components, )
        Percentage of variance explained by each of the selected components.
        If a subset of components is not chosen then all components are stored
        and the sum of explained variances is equal to 1.0.

    Methods
    -------
    fit(traj=None, select='All', start=None, stop=None, step=None)
        Find the principal components of selected atoms in a subset of the
        trajectory.
    transform(traj=None, n_components=None)
        Take the original trajectory and project it onto the principal
        components.
    inverse_tranform(pca_space)
        Take a pca_space and map it back onto the trajectory used to create it.
    """

    def __init__(self, u, select='all'):
        """
        Parameters
        ----------
        u : MDAnalysis Universe
            The universe containing the trajectory frames for Principal
            Component Analysis.
        """
        self._u = u
        self._atoms = self._u.select_atoms(select)
        self._n_atoms = self._atoms.n_atoms
        self._calculated = False


    def fit(self, n_components=None, start=None, stop=None,
            step=None):
        """ Use a subset of frames from the trajectory to generate the
            principal components.

        Parameters
        ----------

        n_components : int, optional
            The number of principal components to be saved, default saves
            all principal components, Default: -1
        start : int, optional
            First frame of trajectory to use for generation
            of covariance matrix, Default: 0
        stop : int, optional
            Last frame of trajectory to use for generation
            of covariance matrix, Default: -1
        step : int, optional
            Step between frames of trajectory to use for generation
            of covariance matrix, Default: 1

        Return
        ------
        cumulated_variance: array, (n_components, )
            The amount of variance explained by the nth component and the n-1
            components preceding it.
        p_components: array, (n_components, n_atoms * 3)

        """
        self._setup_frames(self._u.trajectory, start, stop, step)
        self.n_components = n_components
        self.run()
        self._calculated = True

        if n_components is None:
            return self.cumulated_variance, self.p_components
        else:
            return (self.cumulated_variance[:n_components],
                    self.p_components[:n_components])

    def _prepare(self):
        self._xyz = np.zeros((self.n_frames, self._n_atoms, 3))


    def _single_frame(self):
        self._xyz[self._frame_index] = self._atoms.positions.copy()


    def _conclude(self):
        self._xyz = self._xyz.reshape(self.n_frames, self._n_atoms * 3,
                                      order='F')
        x = self._xyz - self._xyz.mean(0)
        cov = np.cov(x, rowvar = 0)
        e_vals, e_vects = np.linalg.eig(cov)
        sort_idx = np.argsort(e_vals)[::-1]
        self.variance = e_vals[sort_idx]
        self.p_components = e_vects[sort_idx]
        self.cumulated_variance = (np.cumsum(self.variance) /
                                   np.sum(self.variance))


    def transform(self, n_components=None):
        """Apply the dimensionality reduction on a trajectory

        Parameters
        ----------
        n_components: int, optional
            The number of components to be projected onto, default maps onto
            all components.
        Returns
        -------
        pca_space : array, shape (number of frames, number of components)
        """

        xyz = np.array([self._atoms.positions.copy() for ts in self._u.trajectory])
        xyz = xyz.reshape(self._u.trajectory.n_frames,
                          self._n_atoms*3, order='F')
        if n_components is None:
            return np.dot(xyz, self.p_components.T)
        else:
            return np.dot(xyz, self.p_components[:n_components].T)


    def inverse_tranform(self, pca_space):
        """ Transform PCA-transformed data back to original configuration space.

        Parameters
        ----------
        pca_space : array
            The space corresponding to the trajectory coordinates after the PCA
            transformation. Assumes current principal components were the
            components projected onto to create the pca_space.

        Returns
        -------
        original_traj : array, shape (number of frames, number of atoms*3)
        """
        inv = np.linalg.inv(self.p_components)
        return np.dot(pca_space, inv)
