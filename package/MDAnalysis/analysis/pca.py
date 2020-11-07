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
import warnings

import numpy as np
import scipy.integrate

from MDAnalysis import Universe
from MDAnalysis.analysis.align import _fit_to
from MDAnalysis.lib.log import ProgressBar

from ..lib import util
from ..due import due, Doi
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
    p_components: array, (n_atoms * 3, n_components)
        The principal components of the feature space,
        representing the directions of maximum variance in the data.
        The column vector p_components[:, i] is the eigenvector
        corresponding to the variance[i].
    variance : array (n_components, )
        The raw variance explained by each eigenvector of the covariance
        matrix.
    cumulated_variance : array, (n_components, )
        Percentage of variance explained by the selected components and the sum
        of the components preceding it. If a subset of components is not chosen
        then all components are stored and the cumulated variance will converge
        to 1.
    mean_atoms: MDAnalyis atomgroup
        After running :meth:`PCA.run`, the mean position of all the atoms
        used for the creation of the covariance matrix will exist here.

    Methods
    -------
    transform(atomgroup, n_components=None)
        Take an atomgroup or universe with the same number of atoms as was
        used for the calculation in :meth:`PCA.run`, and project it onto the
        principal components.

    Notes
    -----
    Computation can be sped up by supplying a precalculated mean structure.

    .. versionchanged:: 1.0.0
       ``n_components`` now limits the correct axis of ``p_components``.
       ``cumulated_variance`` now accurately represents the contribution of
       each principal component and does not change when ``n_components`` is
       given. If ``n_components`` is not None or is less than the number of
       ``p_components``, ``cumulated_variance`` will not sum to 1.
       ``align=True`` now correctly aligns the trajectory and computes the
       correct means and covariance matrix.

    .. versionchanged:: 0.19.0
       The start frame is used when performing selections and calculating
       mean positions.  Previously the 0th frame was always used.
    """

    def __init__(self, universe, select='all', align=False, mean=None,
                 n_components=None, **kwargs):
        """
        Parameters
        ----------
        universe : Universe
            Universe
        select : string, optional
            A valid selection statement for choosing a subset of atoms from
            the atomgroup.
        align : boolean, optional
            If True, the trajectory will be aligned to a reference
            structure.
        mean : MDAnalysis atomgroup, optional
            An optional reference structure to be used as the mean of the
            covariance matrix.
        n_components : int, optional
            The number of principal components to be saved, default saves
            all principal components
        verbose : bool (optional)
             Show detailed progress of the calculation if set to ``True``.
        """
        super(PCA, self).__init__(universe.trajectory, **kwargs)
        self._u = universe

        # for transform function
        self.align = align

        self._calculated = False
        self._n_components = n_components
        self._select = select
        self._mean = mean

    def _prepare(self):
        # access start index
        self._u.trajectory[self.start]
        # reference will be start index
        self._reference = self._u.select_atoms(self._select)
        self._atoms = self._u.select_atoms(self._select)
        self._n_atoms = self._atoms.n_atoms

        if self._mean is None:
            self.mean = np.zeros(self._n_atoms*3)
            self._calc_mean = True
        else:
            self.mean = self._mean.positions
            self._calc_mean = False

        if self.n_frames == 1:
            raise ValueError('No covariance information can be gathered from a'
                             'single trajectory frame.\n')
        n_dim = self._n_atoms * 3
        self.cov = np.zeros((n_dim, n_dim))
        self._ref_atom_positions = self._reference.positions
        self._ref_cog = self._reference.center_of_geometry()
        self._ref_atom_positions -= self._ref_cog

        if self._calc_mean:
            for ts in ProgressBar(self._u.trajectory[self.start:self.stop:self.step],
                                  verbose=self._verbose, desc="Mean Calculation"):
                if self.align:
                    mobile_cog = self._atoms.center_of_geometry()
                    mobile_atoms, old_rmsd = _fit_to(self._atoms.positions - mobile_cog,
                                                     self._ref_atom_positions,
                                                     self._atoms,
                                                     mobile_com=mobile_cog,
                                                     ref_com=self._ref_cog)

                self.mean += self._atoms.positions.ravel()
            self.mean /= self.n_frames

        self.mean_atoms = self._atoms
        self.mean_atoms.positions = self._atoms.positions

    def _single_frame(self):
        if self.align:
            mobile_cog = self._atoms.center_of_geometry()
            mobile_atoms, old_rmsd = _fit_to(self._atoms.positions - mobile_cog,
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
        self._variance = e_vals[sort_idx]
        self._p_components = e_vects[:, sort_idx]
        self._calculated = True
        self.n_components = self._n_components

    @property
    def n_components(self):
        return self._n_components

    @n_components.setter
    def n_components(self, n):
        if self._calculated:
            if n is None:
                n = len(self._variance)
            self.variance = self._variance[:n]
            self.cumulated_variance = (np.cumsum(self._variance) /
                                       np.sum(self._variance))[:n]
            self.p_components = self._p_components[:, :n]
        self._n_components = n

    def transform(self, atomgroup, n_components=None, start=None, stop=None,
                  step=None):
        """Apply the dimensionality reduction on a trajectory

        Parameters
        ----------
        atomgroup : MDAnalysis atomgroup/ Universe
            The atomgroup or universe containing atoms to be PCA transformed.
        n_components : int, optional
            The number of components to be projected onto, The default
            ``None``maps onto all components.
        start : int, optional
            The frame to start on for the PCA transform. The default
            ``None`` becomes 0, the first frame index.
        stop : int, optional
            Frame index to stop PCA transform. The default ``None`` becomes
            the total number of frames in the trajectory.
            Iteration stops *before* this frame number, which means that the
            trajectory would be read until the end.
        step : int, optional
            Number of frames to skip over for PCA transform. If set to ``None``
            (the default) then every frame is analyzed (i.e., same as
            ``step=1``).

        Returns
        -------
        pca_space : array, shape (n_frames, n_components)

        .. versionchanged:: 0.19.0
           Transform now requires that :meth:`run` has been called before,
           otherwise a :exc:`ValueError` is raised.
        """
        if not self._calculated:
            raise ValueError('Call run() on the PCA before using transform')

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
            dot[i] = np.dot(xyz, self._p_components[:, :dim])

        return dot

    @due.dcite(
        Doi('10.1002/(SICI)1097-0134(19990901)36:4<419::AID-PROT5>3.0.CO;2-U'),
        Doi('10.1529/biophysj.104.052449'),
        description="RMSIP",
        path='MDAnalysis.analysis.pca',
    )
    def rmsip(self, other, n_components=None):
        """Compute the root mean square inner product between subspaces.

        This is only symmetric if the number of components is the same for
        both instances. The RMSIP effectively measures how
        correlated the vectors of this instance are to those of ``other``.

        Please cite [Amadei1999]_ and [Leo-Macias2004]_ if you use this function.

        Parameters
        ----------
        other : :class:`~MDAnalysis.analysis.pca.PCA`
            Another PCA class. This must have already been run.
        n_components : int or tuple of ints, optional
            number of components to compute for the inner products.
            ``None`` computes all of them.

        Returns
        -------
        float:
            Root mean square inner product of the selected subspaces.
            0 indicates that they are mutually orthogonal, whereas 1 indicates
            that they are identical.

        See also
        --------
        rmsip

        .. versionadded:: 1.0.0
        """
        try:
            a = self.p_components
        except AttributeError:
            raise ValueError('Call run() on the PCA before using rmsip')

        try:
            b = other.p_components
        except AttributeError:
            if isinstance(other, type(self)):
                raise ValueError(
                    'Call run() on the other PCA before using rmsip')
            else:
                raise ValueError('other must be another PCA class')

        return rmsip(a.T, b.T, n_components=n_components)

    @due.dcite(
        Doi('10.1016/j.str.2007.12.011'),
        description="Cumulative overlap",
        path='MDAnalysis.analysis.pca',
    )
    def cumulative_overlap(self, other, i=0, n_components=None):
        """Compute the cumulative overlap of a vector in a subspace.

        This is not symmetric. The cumulative overlap measures the overlap of
        the chosen vector in this instance, in the ``other`` subspace.

        Please cite [Yang2008]_ if you use this function.

        Parameters
        ----------
        other : :class:`~MDAnalysis.analysis.pca.PCA`
            Another PCA class. This must have already been run.
        i : int, optional
            The index of eigenvector to be analysed.
        n_components : int, optional
            number of components in ``other`` to compute for the cumulative overlap.
            ``None`` computes all of them.

        Returns
        -------
        float:
            Cumulative overlap of the chosen vector in this instance to
            the ``other`` subspace. 0 indicates that they are mutually
            orthogonal, whereas 1 indicates that they are identical.

        See also
        --------
        cumulative_overlap

        .. versionadded:: 1.0.0

        """

        try:
            a = self.p_components
        except AttributeError:
            raise ValueError(
                'Call run() on the PCA before using cumulative_overlap')

        try:
            b = other.p_components
        except AttributeError:
            if isinstance(other, type(self)):
                raise ValueError(
                    'Call run() on the other PCA before using cumulative_overlap')
            else:
                raise ValueError('other must be another PCA class')

        return cumulative_overlap(a.T, b.T, i=i, n_components=n_components)


def cosine_content(pca_space, i):
    """Measure the cosine content of the PCA projection.

    The cosine content of pca projections can be used as an indicator if a
    simulation is converged. Values close to 1 are an indicator that the
    simulation isn't converged. For values below 0.7 no statement can be made.
    If you use this function please cite [BerkHess1]_.


    Parameters
    ----------
    pca_space : array, shape (number of frames, number of components)
        The PCA space to be analyzed.
    i : int
        The index of the pca_component projection to be analyzed.

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


@due.dcite(
    Doi('10.1002/(SICI)1097-0134(19990901)36:4<419::AID-PROT5>3.0.CO;2-U'),
    Doi('10.1529/biophysj.104.052449'),
    description="RMSIP",
    path='MDAnalysis.analysis.pca',
)
def rmsip(a, b, n_components=None):
    """Compute the root mean square inner product between subspaces.

    This is only symmetric if the number of components is the same for
    ``a`` and ``b``. The RMSIP effectively measures how
    correlated the vectors of ``a`` are to those of ``b``.

    Please cite [Amadei1999]_ and [Leo-Macias2004]_ if you use this function.

    Parameters
    ----------
    a : array, shape (n_components, n_features)
        The first subspace. Must have the same number of features as ``b``.
    b : array, shape (n_components, n_features)
        The second subspace. Must have the same number of features as ``a``.
    n_components : int or tuple of ints, optional
        number of components to compute for the inner products.
        ``None`` computes all of them.

    Returns
    -------
    float:
        Root mean square inner product of the selected subspaces.
        0 indicates that they are mutually orthogonal, whereas 1 indicates
        that they are identical.

    .. versionadded:: 1.0.0
    """
    n_components = util.asiterable(n_components)
    if len(n_components) == 1:
        n_a = n_b = n_components[0]
    elif len(n_components) == 2:
        n_a, n_b = n_components
    else:
        raise ValueError('Too many values provided for n_components')

    if n_a is None:
        n_a = len(a)
    if n_b is None:
        n_b = len(b)

    sip = np.matmul(a[:n_a], b[:n_b].T) ** 2
    msip = sip.sum()/n_a
    return msip**0.5


@due.dcite(
    Doi('10.1016/j.str.2007.12.011'),
    description="Cumulative overlap",
    path='MDAnalysis.analysis.pca',
)
def cumulative_overlap(a, b, i=0, n_components=None):
    """Compute the cumulative overlap of a vector in a subspace.

    This is not symmetric. The cumulative overlap measures the overlap of
    the chosen vector in ``a``, in the ``b`` subspace.

    Please cite [Yang2008]_ if you use this function.

    Parameters
    ----------
    a : array, shape (n_components, n_features) or vector, length n_features
        The first subspace containing the vector of interest. Alternatively,
        the actual vector. Must have the same number of features as ``b``.
    b : array, shape (n_components, n_features)
        The second subspace. Must have the same number of features as ``a``.
    i : int, optional
        The index of eigenvector to be analysed.
    n_components : int, optional
        number of components in ``b`` to compute for the cumulative overlap.
        ``None`` computes all of them.

    Returns
    -------
    float:
        Cumulative overlap of the chosen vector in ``a`` to the ``b`` subspace.
        0 indicates that they are mutually orthogonal, whereas 1 indicates
        that they are identical.

    .. versionadded:: 1.0.0
    """

    if len(a.shape) < len(b.shape):
        a = a[np.newaxis, :]

    vec = a[i][np.newaxis, :]
    vec_norm = (vec**2).sum() ** 0.5

    if n_components is None:
        n_components = len(b)

    b = b[:n_components]
    b_norms = (b**2).sum(axis=1) ** 0.5

    o = np.abs(np.matmul(vec, b.T)) / (b_norms*vec_norm)
    return (o**2).sum() ** 0.5
