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

r"""Principal Component Analysis (PCA) --- :mod:`MDAnalysis.analysis.pca`
=====================================================================

:Authors: John Detlefs
:Year: 2016
:Copyright: GNU Public License v3

.. versionadded:: 0.16.0

This module contains the linear dimensions reduction method Principal Component
Analysis (PCA). PCA sorts a simulation into 3N directions of descending
variance, with N being the number of atoms. These directions are called
the principal components. The dimensions to be analyzed are reduced by only
looking at a few projections of the first principal components. To learn how to
run a Principal Component Analysis, please refer to the :ref:`PCA-tutorial`.

The PCA problem is solved by solving the eigenvalue problem of the covariance
matrix, a :math:`3N \times 3N` matrix where the element :math:`(i, j)` is the
covariance between coordinates :math:`i` and :math:`j`. The principal
components are the eigenvectors of this matrix.

For each eigenvector, its eigenvalue is the variance that the eigenvector
explains. Stored in :attr:`PCA.cumulated_variance`, a ratio for each number of
eigenvectors up to index :math:`i` is provided to quickly find out how many
principal components are needed to explain the amount of variance reflected by
those :math:`i` eigenvectors. For most data, :attr:`PCA.cumulated_variance`
will be approximately equal to one for some :math:`n` that is significantly
smaller than the total number of components, these are the components of
interest given by Principal Component Analysis.

From here, we can project a trajectory onto these principal components and
attempt to retrieve some structure from our high dimensional data.

For a basic introduction to the module, the :ref:`PCA-tutorial` shows how
to perform Principal Component Analysis.

.. _PCA-tutorial:

PCA Tutorial
------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF` and
:data:`~MDAnalysis.tests.datafiles.DCD`). This tutorial shows how to use the
PCA class.

First load all modules and test data

    >>> import MDAnalysis as mda
    >>> import MDAnalysis.analysis.pca as pca
    >>> from MDAnalysis.tests.datafiles import PSF, DCD

Given a universe containing trajectory data we can perform Principal Component
Analyis by using the class :class:`PCA` and retrieving the principal
components.

    >>> u = mda.Universe(PSF, DCD)
    >>> PSF_pca = pca.PCA(u, select='backbone')
    >>> PSF_pca.run()

Inspect the components to determine the principal components you would like
to retain. The choice is arbitrary, but I will stop when 95 percent of the
variance is explained by the components. This cumulated variance by the
components is conveniently stored in the one-dimensional array attribute
``cumulated_variance``. The value at the ith index of `cumulated_variance`
is the sum of the variances from 0 to i.

    >>> n_pcs = np.where(PSF_pca.cumulated_variance > 0.95)[0][0]
    >>> atomgroup = u.select_atoms('backbone')
    >>> pca_space = PSF_pca.transform(atomgroup, n_components=n_pcs)

From here, inspection of the ``pca_space`` and conclusions to be drawn from the
data are left to the user.

Classes and Functions
---------------------

.. autoclass:: PCA
.. autofunction:: cosine_content

"""
from __future__ import division, absolute_import
from six.moves import range
import warnings

import numpy as np
import scipy.integrate

from MDAnalysis import Universe
from MDAnalysis.analysis.align import _fit_to
from MDAnalysis.lib.log import ProgressMeter

from .base import AnalysisBase


class PCA(AnalysisBase):
    """Principal component analysis on an MD trajectory.

    After initializing and calling method with a universe or an atom group,
    principal components ordering the atom coordinate data by decreasing
    variance will be available for analysis. As an example:

        >>> pca = PCA(universe, select='backbone').run()
        >>> pca_space =  pca.transform(universe.select_atoms('backbone'), 3)

    generates the principal components of the backbone of the atomgroup and
    then transforms those atomgroup coordinates by the direction of those
    variances. Please refer to the :ref:`PCA-tutorial` for more detailed
    instructions.

    Attributes
    ----------
    p_components: array, (n_components, n_atoms * 3)
        The principal components of the feature space,
        representing the directions of maximum variance in the data.
    variance : array (n_components, )
        The raw variance explained by each eigenvector of the covariance
        matrix.
    cumulated_variance : array, (n_components, )
        Percentage of variance explained by the selected components and the sum
        of the components preceding it. If a subset of components is not chosen
        then all components are stored and the cumulated variance will converge
        to 1.
    pca_space : array (n_frames, n_components)
        After running :meth:`pca.transform` the projection of the
        positions onto the principal components will exist here.
    mean_atoms: MDAnalyis atomgroup
        After running :meth:`PCA.run`, the mean position of all the atoms
        used for the creation of the covariance matrix will exist here.

    Methods
    -------
    transform(atomgroup, n_components=None)
        Take an atomgroup or universe with the same number of atoms as was
        used for the calculation in :meth:`PCA.run` and project it onto the
        principal components.

    Notes
    -----
    Computation can be speed up by supplying a precalculated mean structure
    """

    def __init__(self, universe, select='all', align=False, mean=None,
                 n_components=None, **kwargs):
        """
        Parameters
        ----------
        universe: Universe
            Universe
        select: string, optional
            A valid selection statement for choosing a subset of atoms from
            the atomgroup.
        align: boolean, optional
            If True, the trajectory will be aligned to a reference
            structure.
        mean: MDAnalysis atomgroup, optional
            An optional reference structure to be used as the mean of the
            covariance matrix.
        n_components : int, optional
            The number of principal components to be saved, default saves
            all principal components, Default: None
        start : int, optional
            First frame of trajectory to use for generation
            of covariance matrix, Default: None
        stop : int, optional
            Last frame of trajectory to use for generation
            of covariance matrix, Default: None
        step : int, optional
            Step between frames of trajectory to use for generation
            of covariance matrix, Default: None
        """
        super(PCA, self).__init__(universe.trajectory,
                                  **kwargs)
        if self.n_frames == 1:
            raise ValueError('No covariance information can be gathered from a'
                             'single trajectory frame.\n')

        self._u = universe

        # for transform function
        self.align = align
        # access 0th index
        self._u.trajectory[0]
        # reference will be 0th index
        self._reference = self._u.select_atoms(select)
        self._atoms = self._u.select_atoms(select)
        self.n_components = n_components
        self._n_atoms = self._atoms.n_atoms
        self._calculated = False

        if mean is None:
            self.mean = np.zeros(self._n_atoms*3)
            self._calc_mean = True
        else:
            self.mean = mean.positions
            self._calc_mean = False

    def _prepare(self):
        n_dim = self._n_atoms * 3
        self.cov = np.zeros((n_dim, n_dim))
        self._ref_atom_positions = self._reference.positions
        self._ref_cog = self._reference.center_of_geometry()
        self._ref_atom_positions -= self._ref_cog

        if self._calc_mean:
            interval = int(self.n_frames // 100)
            interval = interval if interval > 0 else 1
            format = ("Mean Calculation Step"
                      "%(step)5d/%(numsteps)d [%(percentage)5.1f%%]\r")
            mean_pm = ProgressMeter(self.n_frames if self.n_frames else 1,
                                    interval=interval, verbose=self._verbose,
                                    format=format)
            for i, ts in enumerate(self._u.trajectory[self.start:self.stop:
                                                      self.step]):
                if self.align:
                    mobile_cog = self._atoms.center_of_geometry()
                    mobile_atoms, old_rmsd = _fit_to(self._atoms.positions,
                                                     self._ref_atom_positions,
                                                     self._atoms,
                                                     mobile_com=mobile_cog,
                                                     ref_com=self._ref_cog)
                else:
                    self.mean += self._atoms.positions.ravel()
                mean_pm.echo(i)
            self.mean /= self.n_frames

        self.mean_atoms = self._atoms
        self.mean_atoms.positions = self._atoms.positions

    def _single_frame(self):
        if self.align:
            mobile_cog = self._atoms.center_of_geometry()
            mobile_atoms, old_rmsd = _fit_to(self._atoms.positions,
                                             self._ref_atom_positions,
                                             self._atoms,
                                             mobile_com=mobile_cog,
                                             ref_com=self._ref_cog)
            # now all structures are aligned to reference
            x = mobile_atoms.positions.ravel()
        else:
            x = self._atoms.positions.ravel()
        x -= self.mean
        self.cov += np.dot(x[:, np.newaxis], x[:, np.newaxis].T)

    def _conclude(self):
        self.cov /= self.n_frames - 1
        e_vals, e_vects = np.linalg.eig(self.cov)
        sort_idx = np.argsort(e_vals)[::-1]
        self.variance = e_vals[sort_idx]
        self.variance = self.variance[:self.n_components]
        self.p_components = e_vects[:self.n_components, sort_idx]
        self.cumulated_variance = (np.cumsum(self.variance) /
                                   np.sum(self.variance))
        self._calculated = True

    def transform(self, atomgroup, n_components=None, start=None, stop=None,
                  step=None):
        """Apply the dimensionality reduction on a trajectory

        Parameters
        ----------
        atomgroup: MDAnalysis atomgroup/ Universe
            The atomgroup or universe containing atoms to be PCA transformed.
        n_components: int, optional
            The number of components to be projected onto, Default none: maps
            onto all components.
        start: int, optional
            The frame to start on for the PCA transform. Default: None becomes
            0, the first frame index.
        stop: int, optional
            Frame index to stop PCA transform. Default: None becomes n_frames.
            Iteration stops *before* this frame number, which means that the
            trajectory would be read until the end.
        step: int, optional
            Number of frames to skip over for PCA transform. Default: None
            becomes 1.

        Returns
        -------
        pca_space : array, shape (number of frames, number of components)
        """
        if not self._calculated:
            self.run()

        if isinstance(atomgroup, Universe):
            atomgroup = atomgroup.atoms

        if(self._n_atoms != atomgroup.n_atoms):
            raise ValueError('PCA has been fit for'
                             '{} atoms. Your atomgroup'
                             'has {} atoms'.format(self._n_atoms,
                                                   atomgroup.n_atoms))
        if not (self._atoms.types == atomgroup.types).all():
            warnings.warn('Atom types do not match with types used to fit PCA')

        traj = atomgroup.universe.trajectory
        start, stop, step = traj.check_slice_indices(start, stop, step)
        n_frames = len(range(start, stop, step))

        dim = (n_components if n_components is not None else
               self.p_components.shape[1])

        dot = np.zeros((n_frames, dim))

        for i, ts in enumerate(traj[start:stop:step]):
            xyz = atomgroup.positions.ravel() - self.mean
            dot[i] = np.dot(xyz, self.p_components[:, :n_components])

        return dot


def cosine_content(pca_space, i):
    """Measure the cosine content of the PCA projection.

    The cosine content of pca projections can be used as an indicator if a
    simulation is converged. Values close to 1 are an indicator that the
    simulation isn't converged. For values below 0.7 no statement can be made.
    If you use this function please cite [BerkHess1]_.


    Parameters
    ----------
    pca_space: array, shape (number of frames, number of components)
        The PCA space to be analyzed.
    i: int
        The index of the pca_component projectection to be analyzed.

    Returns
    -------
    A float reflecting the cosine content of the ith projection in the PCA
    space. The output is bounded by 0 and 1, with 1 reflecting an agreement
    with cosine while 0 reflects complete disagreement.

    References
    ----------
    .. [BerkHess1] Berk Hess. Convergence of sampling in protein simulations.
                   Phys. Rev. E 65, 031910 (2002).
    """

    t = np.arange(len(pca_space))
    T = len(pca_space)
    cos = np.cos(np.pi * t * (i + 1) / T)
    return ((2.0 / T) * (scipy.integrate.simps(cos*pca_space[:, i])) ** 2 /
            scipy.integrate.simps(pca_space[:, i] ** 2))
