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

"""Diffusion map --- :mod:`MDAnalysis.analysis.diffusionmap`
=====================================================================

:Authors: Eugen Hruska, John Detlefs
:Year: 2016
:Copyright: GNU Public License v2

This module contains the non-linear dimension reduction method diffusion map.
The eigenvectors of a diffusion matrix represent the 'collective coordinates'
of a molecule; the largest eigenvalues are the more dominant collective
coordinates. Assigning phyiscal meaning to the 'collective coordinates' is a
fundamentally difficult problem. The time complexity of the diffusion map is
:math:`O(N^3)`, where N is the number of frames in the trajectory, and the in-memory
storage complexity is :math:`O(N^2)`. Instead of a single trajectory a sample of
protein structures can be used. The sample should be equiblibrated, at least
locally. The order of the sampled structures in the trajectory is irrelevant.

The :ref:`Diffusion-Map-tutorial` shows how to use diffusion map for dimension
reduction.

More details about diffusion maps are in [deLaPorte1]_, [Lafon1]_ ,
[Ferguson1]_, and [Clementi1]_.

.. _Diffusion-Map-tutorial:

Diffusion Map tutorial
----------------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF` and
:data:`~MDAnalysis.tests.datafiles.DCD`). This tutorial shows how to use the
Diffusion Map class.

First load all modules and test data ::

   >>> import MDAnalysis as mda
   >>> import MDAnalysis.analysis.diffusionmap as diffusionmap
   >>> from MDAnalysis.tests.datafiles import PSF, DCD

Given a universe or atom group, we can create and eigenvalue decompose
the Diffusion Matrix from that trajectory using :class:`DiffusionMap`:: and get
the corresponding eigenvalues and eigenvectors.

   >>> u = mda.Universe(PSF,DCD)

We leave determination of the appropriate scale parameter epsilon to the user,
[Clementi1]_ uses a complex method involving the k-nearest-neighbors of a
trajectory frame, whereas others simple use a trial-and-error approach with
a constant epsilon. Currently, the constant epsilon method is implemented
by MDAnalysis.

   >>> dmap = diffusionmap.DiffusionMap(u, select='backbone', epsilon=2)
   >>> dmap.run()

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

   >>>  # first cell of  a jupyter notebook should contain: %matplotlib inline
   >>>  import matplotlib.pyplot as plt
   >>>  f, ax = plt.subplots()
   >>>  upper_limit = # some reasonably high number less than the n_eigenvectors
   >>>  ax.plot(dmap.eigenvalues[:upper_limit])
   >>>  ax.set(xlabel ='eigenvalue index', ylabel='eigenvalue')
   >>>  plt.tight_layout()

From here we can transform into the diffusion space

   >>> num_eigenvectors = # some number less than the number of frames after
   >>> # inspecting for the spectral gap
   >>> fit = dmap.transform(num_eigenvectors, time=1)

It can be difficult to interpret the data, and is left as a task
for the user. The `diffusion distance` between frames i and j is best
approximated by the euclidean distance  between rows i and j of
self.diffusion_space.


.. _Distance-Matrix-tutorial:

Distance Matrix tutorial
------------------------

Often a, a custom distance matrix could be useful for local
epsilon determination or other manipulations on the diffusion
map method. The :class:`DistanceMatrix` exists in
:mod:`~MDAnalysis.analysis.diffusionmap` and can be passed
as an initialization argument for :class:`DiffusionMap`.

    >>> import MDAnalysis as mda
    >>> import MDAnalysis.analysis.diffusionmap as diffusionmap
    >>> from MDAnalysis.tests.datafiles import PSF, DCD

Now create the distance matrix and pass it as an argument to
:class:`DiffusionMap`.

    >>> u = mda.Universe(PSF,DCD)
    >>> dist_matrix = diffusionmap.DistanceMatrix(u, select='all')
    >>> dist_matrix.run()
    >>> dmap = diffusionmap.DiffusionMap(dist_matrix)
    >>> dmap.run()

Classes
-------

.. autoclass:: DiffusionMap
.. autoclass:: DistanceMatrix

References
----------

If you use this Dimension Reduction method in a publication, please
cite [Lafon1]_.

If you choose the default metric, this module uses the fast QCP algorithm
[Theobald2005]_ to calculate the root mean square distance (RMSD) between two
coordinate sets (as implemented in
:func:`MDAnalysis.lib.qcprot.CalcRMSDRotationalMatrix`).  When using this
module in published work please cite [Theobald2005]_.


.. [Lafon1] Coifman, Ronald R., Lafon, Stephane. Diffusion
            maps. Appl. Comput. Harmon.  Anal. 21, 5–30 (2006).

.. [deLaPorte1] J. de la Porte, B. M. Herbst, W. Hereman, S. J. van der Walt.
             An Introduction to Diffusion Maps. In: The 19th Symposium of the
             Pattern Recognition Association of South Africa (2008).

.. [Clementi1] Rohrdanz, M. A, Zheng, W, Maggioni, M, & Clementi, C.
             Determination of reaction coordinates via locally scaled diffusion
             map. J. Chem. Phys. 134, 124116 (2011).

.. [Ferguson1] Ferguson, A. L.; Panagiotopoulos, A. Z.; Kevrekidis, I. G.
             Debenedetti, P. G. Nonlinear dimensionality reduction in molecular
             simulation: The diffusion map approach Chem. Phys. Lett. 509, 1−11
             (2011)

"""
from __future__ import division, absolute_import
from six.moves import range
import logging
import warnings

import numpy as np

from MDAnalysis.core.universe import Universe
from .rms import rmsd
from .base import AnalysisBase
from MDAnalysis.lib.util import deprecate

logger = logging.getLogger("MDAnalysis.analysis.diffusionmap")


class DistanceMatrix(AnalysisBase):
    """Calculate the pairwise distance between each frame in a trajectory
    using a given metric

    A distance matrix can be initialized on its own and used as an
    initialization argument in :class:`DiffusionMap`. Refer to the
    :ref:`Distance-Matrix-tutorial` for a demonstration.

    Attributes
    ----------
    atoms : AtomGroup
        Selected atoms in trajectory subject to dimension reduction
    dist_matrix : array, (n_frames, n_frames)
        Array of all possible ij metric distances between frames in trajectory.
        This matrix is symmetric with zeros on the diagonal.

    """
    def __init__(self, u, select='all', metric=rmsd, cutoff=1E0-5,
                 weights=None, start=None, stop=None, step=None,
                 **kwargs):
        """
        Parameters
        ----------
        u : universe `~MDAnalysis.core.universe.Universe`
            The MD Trajectory for dimension reduction, remember that
            computational cost of eigenvalue decomposition
            scales at O(N^3) where N is the number of frames.
            Cost can be reduced by increasing step interval or specifying a
            start and stop.
        select: str, optional
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
        start : int, optional
            First frame of trajectory to analyse, Default: None becomes 0.
        stop : int, optional
            Frame index to stop analysis. Default: None becomes
            n_frames. Iteration stops *before* this frame number,
            which means that the trajectory would be read until the end.
        step : int, optional
            Step between frames to analyse, Default: None becomes 1.
        verbose : bool (optional)
             Show detailed progress of the calculation if set to ``True``; the
             default is ``False``.
        """
        self._u = u
        traj = self._u.trajectory

        # remember that this must be called before referencing self.n_frames
        super(DistanceMatrix, self).__init__(self._u.trajectory,
                                           start=start, stop=stop, step=step, **kwargs)

        self.atoms = self._u.select_atoms(select)
        self._metric = metric
        self._cutoff = cutoff
        self._weights = weights
        self._calculated = False

    def _prepare(self):
        self.dist_matrix = np.zeros((self.n_frames, self.n_frames))

    def _single_frame(self):
        iframe = self._ts.frame
        i_ref = self.atoms.positions
        # diagonal entries need not be calculated due to metric(x,x) == 0 in
        # theory, _ts not updated properly. Possible savings by setting a
        # cutoff for significant decimal places to sparsify matrix
        for j, ts in enumerate(self._u.trajectory[iframe:self.stop:self.step]):
            self._ts = ts
            j_ref = self.atoms.positions
            dist = self._metric(i_ref, j_ref, weights=self._weights)
            self.dist_matrix[self._frame_index, j+self._frame_index] = (
                dist if dist > self._cutoff else 0)
            self.dist_matrix[j+self._frame_index, self._frame_index] = (
                self.dist_matrix[self._frame_index, j+self._frame_index])
        self._ts = self._u.trajectory[iframe]

    def _conclude(self):
        self._calculated = True

    @deprecate(release="0.19.0", remove="1.0.0",
               message="Use ``np.save(filename, DistanceMatrix.dist_matrix)`` instead.")
    def save(self, filename):
        """save squared distance matrix

        Parameters
        ----------
        outfile : str
            file to save distance matrix

        """
        np.save(filename, self.dist_matrix)
        logger.info("Wrote the distance-squared matrix to file %r", filename)


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
    """

    def __init__(self, u, epsilon=1, **kwargs):
        """
        Parameters
        -------------
        u : MDAnalysis Universe or DistanceMatrix object
            Can be a Universe, in which case one must supply kwargs for the
            initialization of a DistanceMatrix. Otherwise, this can be a
            DistanceMatrix already initialized. Either way, this will be made
            into a diffusion kernel.
        epsilon : Float
            Specifies the method used for the choice of scale parameter in the
            diffusion map. More information in [Lafon1]_, [Ferguson1]_ and
            [Clementi1]_, Default: 1.
        **kwargs
            Parameters to be passed for the initialization of a
            :class:`DistanceMatrix`.
            """
        if isinstance(u, Universe):
            self._dist_matrix = DistanceMatrix(u, **kwargs)
        elif isinstance(u, DistanceMatrix):
            self._dist_matrix = u
        else:
            raise ValueError("U is not a Universe or DistanceMatrix and"
                             " so the DiffusionMap has no data to work with.")
        self._epsilon = epsilon
        # important for transform function and length of .run() method
        self._n_frames = self._dist_matrix.n_frames
        if self._n_frames > 5000:
            warnings.warn("The distance matrix is very large, and can "
                          "be very slow to compute. Consider picking a larger "
                          "step size in distance matrix initialization.")


    def run(self):
        """ Create and decompose the diffusion matrix in preparation
        for a diffusion map."""
        # run only if distance matrix not already calculated
        if not self._dist_matrix._calculated:
            self._dist_matrix.run()
        self._scaled_matrix = (self._dist_matrix.dist_matrix ** 2 /
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
            The diffusion map embedding as defined by [Ferguson1]_.
        """
        return (self._eigenvectors[1:n_eigenvectors+1,].T *
                (self.eigenvalues[1:n_eigenvectors+1]**time))
