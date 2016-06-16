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
Diffusion map --- :mod:`MDAnalysis.analysis.diffusionmap`
=====================================================================

:Authors: Eugen Hruska, John Detlefs
:Year: 2016
:Copyright: GNU Public License v3

The module contains the non-linear dimension reduction method diffusion map.
The diffusion map provides an estimate of the slowest collective
coordinates for a trajectory. This non-linear dimension reduction method
assumes that the trajectory is long enough to represent a probability
distribution of as protein close to the equilibrium. Furthermore, the diffusion
map assumes that the diffusion coefficients are constant. The eigenvectors with
the largest eigenvalues are the more dominant collective coordinates. Assigning
phyiscal meaning to the reaction coordinates is a fundamentally difficult
problem. The time complexity of the diffusion map is O(N^3), where N is the
number of frames in the trajectory, and the storage complexity is O(N^3).
Instead of a single trajectory a sample of protein structures
can be used. The sample should be equiblibrated, at least locally. Different
weights can be used to determine the anisotropy of the diffusion kernel.
The motivation for the creation of an anisotropic kernel is given on page 14 of
[Lafon1]_. The order of the sampled structures in the trajectory is irrelevant.

The :ref:`Diffusion-Map-tutorial` shows how to use diffusion map for dimension
reduction.

More details about diffusion maps are in [Lafon1]_ , [Ferguson1], and [Clementi1]_.

.. _Diffusion-Map-tutorial:

Diffusion Map tutorial
--------------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF` and
:data:`~MDAnalysis.tests.datafiles.DCD`). Notice that this trajectory
does not represent a sample from equilibrium. This violates a basic
assumption of the anisotropic kernel. The anisotropic kernel is created to act
as a discrete approximation for the Fokker-Plank equation from statistical
mechanics, results from a non-equilibrium trajectory shouldn't be interpreted
for this reason. This tutorial shows how to use the diffusionmap function.
First load all modules and test data ::

   >>> import MDAnalysis
   >>> import numpy as np
   >>> import MDAnalysis.analysis.diffusionmap as diffusionmap
   >>> from MDAnalysis.tests.datafiles import PSF,DCD

Given a universe or atom group, we can calculate the diffusion map from
that trajectory using :class:`DiffusionMap`:: and get the corresponding
eigenvalues and eigenvectors. Fr


   >>> u = MDAnalysis.Universe(PSF,DCD)
   >>> d_matrix = DistMatrix(u)
   >>> dmap = diffusionmap.DiffusionMap(DistMatrix)
   >>> dmap.run()
   >>> eigenvalues = dmap.eigenvalues
   >>> eigenvectors = dmap.eigenvectors
From here we can
   >>> ts_one = u.trajectory[0]
   >>> params_on_frame = dmap.embedding(ts_one)

Classes
-------

.. autoclass:: DiffusionMap
   :members:
.. autoclass:: DistMatrix
.. autoclass:: Epsilon

References
---------

If you use this Dimension Reduction method in a publication, please
reference:
..[Ferguson1] Ferguson, A. L.; Panagiotopoulos, A. Z.; Kevrekidis, I. G.;
Debenedetti, P. G. Chem. Phys. Lett. 2011, 509, 1−11
..[Lafon1] Coifman, Ronald R., Lafon, Stephane (2006) Diffusion maps.
Appl. Comput. Harmon. Anal. 21, 5–30.
..[Clementi1] Rohrdanz, M. A, Zheng, W, Maggioni, M, & Clementi, C. (2013)
Determination of reaction coordinates via locally scaled
diffusion map. Journal of Chemical Physics.

.. If you choose the default metric, this module uses the fast QCP algorithm
[Theobald2005]_ to calculate the root mean square distance (RMSD) between
two coordinate sets (as implemented in :func:`MDAnalysis.lib.qcprot.CalcRMSDRotationalMatrix`).
When using this module in published work please cite [Theobald2005]_.


"""
from six.moves import range
import logging

import numpy as np

from MDAnalysis.analysis import rms

from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.diffusionmap")


class DistMatrix(AnalysisBase):
    """ Calculate the pairwise distance between each frame in a trajectory using
        a given metric

    Attributes
    ----------
    atoms : AtomGroup
        Selected atoms in trajectory subject to dimension reduction
    dist_matrix : array
        Array of all possible ij metric distances between frames in trajectory.
        This matrix is symmetric with zeros on the diagonal.

    Methods
    -------
    save(filename)
        Save the `dist_matrix` to a given filename
    """
    def __init__(self, u, select='all', metric=rms.rmsd, cutoff=1E0-5,
                 weights=None, start=None, stop=None, step=None):
        """
        Parameters
        ----------
        u : trajectory `~MDAnalysis.core.AtomGroup.Universe`
            The MD Trajectory for dimension reduction, remember that
            computational cost of eigenvalue decomposition
            scales at O(N^3) where N is the number of frames.
            Cost can be reduced by increasing step interval or specifying a
            start and stop.
        select: str, optional
            Any valid selection string for
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms`
            This selection of atoms is used to calculate the RMSD between
            different frames. Water should be excluded.
        metric : function, optional
            Maps two numpy arrays to a float, is positive definite and
            symmetric, Default: metric is set to rms.rmsd().
        cutoff : float, optional
            Specify a given cutoff for metric values to be considered equal,
            Default: 1EO-5
        weights : array, optional
            Weights to be given to coordinates for metric calculation
        start : int, optional
            First frame of trajectory to analyse, Default: 0
        stop : int, optional
            Last frame of trajectory to analyse, Default: -1
        step : int, optional
            Step between frames to analyse, Default: 1
        """

        self._u = u
        traj = self._u.trajectory
        self.atoms = self._u.select_atoms(select)
        self._metric = metric
        self._cutoff = cutoff
        self._weights = weights
        # remember that this must be called before referencing self.nframes
        self._setup_frames(traj, start, stop, step)

    def _prepare(self):
        logger.info('numframes {0}'.format(self.nframes))
        self.dist_matrix = np.zeros((self.nframes, self.nframes))

    def _single_frame(self):
        traj_index = self._ts.frame
        i_ref = self.atoms.positions - self.atoms.center_of_mass()
        # diagonal entries need not be calculated due to metric(x,x) == 0 in
        # theory, _ts not updated properly. Possible savings by setting a
        # cutoff for significant decimal places to sparsify matrix
        for j in range(self.stop-1, self._ts.frame-1, -self.step):
            self._ts = self._u.trajectory[j]
            j_ref = self.atoms.positions-self.atoms.center_of_mass()
            dist = self._metric(i_ref, j_ref, weights=self._weights)
            # distance squared in preparation for kernel calculation
            # don't think this will come up in other areas, but might be
            # a fix we should make later
            self.dist_matrix[traj_index, j] = dist**2 if dist > self._cutoff
            else 0
            self.dist_matrix[j, traj_index] = self.dist_matrix[traj_index, j]

    def save(self, filename):
        np.savetxt(filename, self.dist_matrix)
        logger.info("Wrote the distance-squared matrix to file %r", filename)


class Epsilon(object):
    """Manipulates a distance matrix by local scale parameters

    Attributes
    ----------
    scaledMatrix : DistMatrix object
        A matrix with each term divided by a local scale parameter

    Methods
    -------
    determineEpsilon()
        Determine local scale parameters using a chosen algorithm
    """
    def __init__(self, DistMatrix, **kwargs):
        self._distMatrix = DistMatrix

    def determineEpsilon(self):
        pass


class DiffusionMap(DistMatrix):
    """Non-linear dimension reduction method

    Dimension reduction with diffusion map of the structures in the universe.

    Attributes
    ----------
    eigenvalues: array
        Eigenvalues of the diffusion map
    eigenvectors: array
        Eigenvectors of the diffusion map

    Methods
    -------
    spectral_gap()
        Retrieve a guess for the set of eigenvectors reflecting the intrinsic
        dimensionality of the molecular system.

    embedding(timestep)
        Perform an embedding of a frame into the eigenvectors representing
        the collective coordinates.
    """

    def __init__(self, DistMatrix, epsilon=EpsilonConstant, weights=None):
        """
        Parameters
        -------------
        DistMatrix : DistMatrix object
            Distance matrix to be made into a diffusion kernel and perform
            an eigenvector decomposition on.
        epsilon : function, optional
            Specifies the method used for the choice of scale parameter in the
            diffusion map. More information in [1] and [2]. With 'average'
            the average of the RMSD to the k-th greatest value is used.
        weights: list, optional
            The list has to have the same length as the trajectory.
            With 'None' the weight of each frame of the trajectory will be the
            same.
        """
        self._epsilon = epsilon
        # this should be a reference to the same object as
        # epsilon.scaledMatrix
        self.dist_matrix = DistMatrix.dist_matrix
        if weights_ker is None:
            # weights do not apply to metric
            self._weights_ker = np.ones((self.nframes,))
        else:
            if weights.shape[0] != self.nframes:
                raise ValueError("The weight should have the same length as "
                                 'the trajectory')
            else:
                # weights are constructed as relative to the mean
                self._weights_ker = (np.asarray(weights, dtype=np.float64) /
                                     np.mean(weights))

    def _conclude(self):
        self._epsilon.determineEpsilon()
        self._kernel = self._epsilon.scaledMatrix
        # take negative exponent of scaled matrix to create Isotropic kernel
        self._kernel = np.exp(-self._kernel)

        # define an anistropic diffusion term q
        q_vector = np.zeros((self.nframes, ))
        # weights should reflect the density of the points on the manifold
        for i in range(self.nframes):
            q_vector[i] = np.dot(self._kernel2[i, :], self._weights_ker)

        # Form a new kernel from the anisotropic diffusion term q
        self._kernel /= np.sqrt(q_vector[:, np.newaxis].dot(q_vector[np.newaxis]))

        # Weighted Graph Laplacian normalization on this new graph
        d_vector = np.zeros((self.nframes, ))
        for i in range(self.nframes):
            d_vector[i] = np.dot(self._kernel2[i, :], self._weights_ker)

        for i in range(self.nframes):
            self._kernel[i, :] = self._kernel2[i, :] * self._weights_ker

        # Define anisotropic transition kernel from this term
        self._kernel /= np.sqrt(d_vector[:, np.newaxis].dot(d_vector[np.newaxis]))

        eigenvals, eigenvectors = np.linalg.eig(self._kernel)

        eg_arg = np.argsort(eigenvals)
        self.eigenvalues = eigenvals[eg_arg[::-1]]
        self.eigenvectors = eigenvectors[eg_arg[::-1], :]

    def spectral_gap(self):
        # TODO
        pass

    def embedding(self, timestep):
        # TODO
        pass


class EpsilonConstant(Epsilon):
    """Premade function for the determination of epsilon based on providing
    an arbitrary constant"""

    def __init__(self, DistMatrix, epsilon):
        """
        Parameters
        ----------
        epsilon : int
            The value of epsilon to be used as a local scale parameter
        """
        self._epsilon = epsilon
        epsilon = np.full((nframes, ), epsilon)
