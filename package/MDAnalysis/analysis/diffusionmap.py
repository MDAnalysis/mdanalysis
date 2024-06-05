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

"""Diffusion map --- :mod:`MDAnalysis.analysis.diffusionmap`
=====================================================================

:Authors: Eugen Hruska, John Detlefs
:Year: 2016
:Copyright: GNU Public License v2

This module contains the non-linear dimension reduction method diffusion map.
The eigenvectors of a diffusion matrix represent the 'collective coordinates'
of a molecule; the largest eigenvalues are the more dominant collective
coordinates. Assigning physical meaning to the 'collective coordinates' is a
fundamentally difficult problem. The time complexity of the diffusion map is
:math:`O(N^3)`, where N is the number of frames in the trajectory, and the in-memory
storage complexity is :math:`O(N^2)`. Instead of a single trajectory a sample of
protein structures can be used. The sample should be equilibrated, at least
locally. The order of the sampled structures in the trajectory is irrelevant.

The :ref:`Diffusion-Map-tutorial` shows how to use diffusion map for dimension
reduction.

More details about diffusion maps are in
:footcite:p:`deLaPorte2008,Coifman-Lafon,Ferguson2011,Clementi2011`.

.. _Diffusion-Map-tutorial:

Diffusion Map tutorial
----------------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF` and
:data:`~MDAnalysis.tests.datafiles.DCD`). This tutorial shows how to use the
Diffusion Map class.

First load all modules and test data

.. code-block:: python

   import MDAnalysis as mda
   import MDAnalysis.analysis.diffusionmap as diffusionmap
   from MDAnalysis.tests.datafiles import PSF, DCD

Given a universe or atom group, we can create and eigenvalue decompose
the Diffusion Matrix from that trajectory using :class:`DiffusionMap` and get
the corresponding eigenvalues and eigenvectors.

.. code-block:: python

   u = mda.Universe(PSF,DCD)

We leave determination of the appropriate scale parameter epsilon to the user,
:footcite:p:`Clementi2011` uses a complex method involving the k-nearest-neighbors
of a trajectory frame, whereas others simple use a trial-and-error approach
with a constant epsilon. Currently, the constant epsilon method is implemented
by MDAnalysis.

.. code-block:: python

   dmap = diffusionmap.DiffusionMap(u, select='backbone', epsilon=2)
   dmap.run()

From here we can perform an embedding onto the k dominant eigenvectors. The
non-linearity of the map means there is no explicit relationship between the
lower dimensional space and our original trajectory. However, this is an
isometry (distance preserving map), which means that points close in the lower
dimensional space are close in the higher-dimensional space and vice versa.
In order to embed into the most relevant low-dimensional space, there should
exist some number of dominant eigenvectors, whose corresponding eigenvalues
diminish at a constant rate until falling off, this is referred to as a
spectral gap and should be somewhat apparent for a system at equilibrium with a
high number of frames.

.. code-block:: python

   import matplotlib.pyplot as plt
   f, ax = plt.subplots()
   upper_limit = # some reasonably high number less than the n_eigenvectors
   ax.plot(dmap.eigenvalues[:upper_limit])
   ax.set(xlabel ='eigenvalue index', ylabel='eigenvalue')
   plt.tight_layout()

From here we can transform into the diffusion space

.. code-block:: python

   num_eigenvectors = # some number less than the number of frames after
   # inspecting for the spectral gap
   fit = dmap.transform(num_eigenvectors, time=1)

It can be difficult to interpret the data, and is left as a task
for the user. The `diffusion distance` between frames i and j is best
approximated by the euclidean distance  between rows i and j of
self.diffusion_space.

Classes
-------

.. autoclass:: DiffusionMap
.. autoclass:: DistanceMatrix

References
----------

If you use this Dimension Reduction method in a publication, please
cite :footcite:p:`Coifman-Lafon`.

If you choose the default metric, this module uses the fast QCP algorithm
:footcite:p:`Theobald2005` to calculate the root mean square distance (RMSD)
between two coordinate sets (as implemented in
:func:`MDAnalysis.lib.qcprot.CalcRMSDRotationalMatrix`).  When using this
module in published work please :footcite:p:`Theobald2005`.


.. footbibliography::
"""
import logging
import warnings

import numpy as np

from MDAnalysis.core.universe import Universe
from MDAnalysis.core.groups import AtomGroup, UpdatingAtomGroup
from .rms import rmsd
from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.diffusionmap")


class DistanceMatrix(AnalysisBase):
    """Calculate the pairwise distance between each frame in a trajectory
    using a given metric

    A distance matrix can be initialized on its own and used as an
    initialization argument in :class:`DiffusionMap`.

    Parameters
    ----------
    universe : `~MDAnalysis.core.universe.Universe` or `~MDAnalysis.core.groups.AtomGroup`
        The MD Trajectory for dimension reduction, remember that
        computational cost of eigenvalue decomposition
        scales at O(N^3) where N is the number of frames.
        Cost can be reduced by increasing step interval or specifying a
        start and stop value when calling :meth:`DistanceMatrix.run`.
    select : str, optional
        Any valid selection string for
        :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms`
        This selection of atoms is used to calculate the RMSD between
        different frames. Water should be excluded.
    metric : function, optional
        Maps two numpy arrays to a float, is positive definite and
        symmetric. The API for a metric requires that the arrays must have
        equal length, and that the function should have weights as an
        optional argument. Weights give each index value its own weight for
        the metric calculation over the entire arrays. Default: metric is
        set to rms.rmsd().
    cutoff : float, optional
        Specify a given cutoff for metric values to be considered equal,
        Default: 1EO-5
    weights : array, optional
        Weights to be given to coordinates for metric calculation
    verbose : bool, optional
            Show detailed progress of the calculation if set to ``True``; the
            default is ``False``.

    Attributes
    ----------
    atoms : `~MDAnalysis.core.groups.AtomGroup`
        Selected atoms in trajectory subject to dimension reduction
    results.dist_matrix : numpy.ndarray, (n_frames, n_frames)
        Array of all possible ij metric distances between frames in trajectory.
        This matrix is symmetric with zeros on the diagonal.

        .. versionadded:: 2.0.0

    dist_matrix : numpy.ndarray, (n_frames, n_frames)

        .. deprecated:: 2.0.0
                 Will be removed in MDAnalysis 3.0.0. Please use
                 :attr:`results.dist_matrix` instead.

    Example
    -------
    Often, a custom distance matrix could be useful for local
    epsilon determination or other manipulations on the diffusion
    map method. The :class:`DistanceMatrix` exists in
    :mod:`~MDAnalysis.analysis.diffusionmap` and can be passed
    as an initialization argument for :class:`DiffusionMap`.

    .. code-block:: python

        import MDAnalysis as mda
        import MDAnalysis.analysis.diffusionmap as diffusionmap
        from MDAnalysis.tests.datafiles import PSF, DCD

    Now create the distance matrix and pass it as an argument to
    :class:`DiffusionMap`.

        u = mda.Universe(PSF,DCD)
        dist_matrix = diffusionmap.DistanceMatrix(u, select='all')
        dist_matrix.run()
        dmap = diffusionmap.DiffusionMap(dist_matrix)
        dmap.run()

    .. versionchanged:: 1.0.0
       ``save()`` method has been removed. You can use ``np.save()`` on
       :attr:`DistanceMatrix.results.dist_matrix` instead.
    .. versionchanged:: 2.0.0
         :attr:`dist_matrix` is now stored in a
         :class:`MDAnalysis.analysis.base.Results` instance.
    .. versionchanged:: 2.2.0
         :class:`DistanceMatrix` now also accepts `AtomGroup`.
    .. versionchanged:: 2.8.0
         :class:`DistanceMatrix` is now correctly works with `frames=...`
         parameter (#4432) by iterating over `self._sliced_trajectory`
    """
    def __init__(self, universe, select='all', metric=rmsd, cutoff=1E0-5,
                 weights=None, **kwargs):
        # remember that this must be called before referencing self.n_frames
        super(DistanceMatrix, self).__init__(universe.universe.trajectory,
                                             **kwargs)

        if isinstance(universe, UpdatingAtomGroup):
            wmsg = ("U must be a static AtomGroup. Parsing an updating AtomGroup "
                    "will result in a static AtomGroup with the current frame "
                    "atom selection.")
            warnings.warn(wmsg)

        self.atoms = universe.select_atoms(select)
        self._metric = metric
        self._cutoff = cutoff
        self._weights = weights
        self._calculated = False

    def _prepare(self):
        self.results.dist_matrix = np.zeros((self.n_frames, self.n_frames))

    def _single_frame(self):
        iframe = self._frame_index
        i_ref = self.atoms.positions
        # diagonal entries need not be calculated due to metric(x,x) == 0 in
        # theory, _ts not updated properly. Possible savings by setting a
        # cutoff for significant decimal places to sparsify matrix
        for j, ts in enumerate(self._sliced_trajectory[iframe:]):
            self._ts = ts
            j_ref = self.atoms.positions
            dist = self._metric(i_ref, j_ref, weights=self._weights)
            self.results.dist_matrix[self._frame_index,
                                     j+self._frame_index] = (
                                            dist if dist > self._cutoff else 0)
            self.results.dist_matrix[j+self._frame_index,
                                     self._frame_index] = (
                                self.results.dist_matrix[self._frame_index,
                                                         j+self._frame_index])
        self._ts = self._sliced_trajectory[iframe]

    @property
    def dist_matrix(self):
        wmsg = ("The `dist_matrix` attribute was deprecated in "
                "MDAnalysis 2.0.0 and will be removed in MDAnalysis 3.0.0. "
                "Please use `results.dist_matrix` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.dist_matrix

    def _conclude(self):
        self._calculated = True


class DiffusionMap(object):
    """Non-linear dimension reduction method

    Dimension reduction with diffusion mapping of selected structures in a
    trajectory.

    Attributes
    ----------
    eigenvalues: array (n_frames,)
        Eigenvalues of the diffusion map

    Methods
    -------
    run()
        Constructs an anisotropic diffusion kernel and performs eigenvalue
        decomposition on it.
    transform(n_eigenvectors, time)
        Perform an embedding of a frame into the eigenvectors representing
        the collective coordinates.


    .. versionchanged:: 2.2.0
         :class:`DiffusionMap` now also accepts `AtomGroup`.
    """

    def __init__(self, u, epsilon=1, **kwargs):
        """
        Parameters
        -------------
        u : MDAnalysis Universe or AtomGroup or DistanceMatrix object.
            Can be a Universe or AtomGroup, in which case one must supply kwargs for the
            initialization of a DistanceMatrix. Otherwise, this can be a
            DistanceMatrix already initialized. Either way, this will be made
            into a diffusion kernel.
        epsilon : Float
            Specifies the method used for the choice of scale parameter in the
            diffusion map. More information in
            :footcite:p:`Coifman-Lafon,Ferguson2011,Clementi2011`, Default: 1.
        **kwargs
            Parameters to be passed for the initialization of a
            :class:`DistanceMatrix`.
        """
        if isinstance(u, AtomGroup) or isinstance(u, Universe):
            self._dist_matrix = DistanceMatrix(u, **kwargs)
        elif isinstance(u, DistanceMatrix):
            self._dist_matrix = u
        else:
            raise ValueError("U is not a Universe or AtomGroup or DistanceMatrix and"
                             " so the DiffusionMap has no data to work with.")
        self._epsilon = epsilon

    def run(self, start=None, stop=None, step=None):
        """ Create and decompose the diffusion matrix in preparation
        for a diffusion map.

        Parameters
        ----------
        start : int, optional
            start frame of analysis
        stop : int, optional
            stop frame of analysis
        step : int, optional
            number of frames to skip between each analysed frame

        .. versionchanged:: 0.19.0
           Added start/stop/step kwargs
        """
        # run only if distance matrix not already calculated
        if not self._dist_matrix._calculated:
            self._dist_matrix.run(start=start, stop=stop, step=step)
        # important for transform function and length of .run() method
        self._n_frames = self._dist_matrix.n_frames
        if self._n_frames > 5000:
            warnings.warn("The distance matrix is very large, and can "
                          "be very slow to compute. Consider picking a larger "
                          "step size in distance matrix initialization.")
        self._scaled_matrix = (self._dist_matrix.results.dist_matrix ** 2 /
                               self._epsilon)
        # take negative exponent of scaled matrix to create Isotropic kernel
        self._kernel = np.exp(-self._scaled_matrix)
        D_inv = np.diag(1 / self._kernel.sum(1))
        self._diff = np.dot(D_inv, self._kernel)
        self._eigenvals, self._eigenvectors = np.linalg.eig(self._diff)
        sort_idx = np.argsort(self._eigenvals)[::-1]
        self.eigenvalues = self._eigenvals[sort_idx]
        self._eigenvectors = self._eigenvectors[sort_idx]
        self._calculated = True
        return self

    def transform(self, n_eigenvectors, time):
        """ Embeds a trajectory via the diffusion map

        Parameters
        ---------
        n_eigenvectors : int
            The number of dominant eigenvectors to be used for
            diffusion mapping
        time : float
            Exponent that eigenvalues are raised to for embedding, for large
            values, more dominant eigenvectors determine diffusion distance.
        Return
        ------
        diffusion_space : array (n_frames, n_eigenvectors)
            The diffusion map embedding as defined by :footcite:p:`Ferguson2011`.
        """
        return (self._eigenvectors[1:n_eigenvectors+1,].T *
                (self.eigenvalues[1:n_eigenvectors+1]**time))
