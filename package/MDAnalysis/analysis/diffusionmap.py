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
distribution of a protein close to the equilibrium. Furthermore, the diffusion
map assumes that the diffusion coefficients associated with the dynamical
motion of molecules in the system are constant. The eigenvectors with
the largest eigenvalues are the more dominant collective coordinates. Assigning
phyiscal meaning of the 'collective coordinates' is a fundamentally difficult
problem. The time complexity of the diffusion map is O(N^3), where N is the
number of frames in the trajectory, and the in-memory storage complexity is
O(N^2). Instead of a single trajectory a sample of protein structures
can be used. The sample should be equiblibrated, at least locally. Different
weights can be used to determine the anisotropy of the diffusion kernel.
The motivation for the creation of an anisotropic kernel is given on page 14 of
[Lafon1]_. The order of the sampled structures in the trajectory is irrelevant.

The :ref:`Diffusion-Map-tutorial` shows how to use diffusion map for dimension
reduction.

More details about diffusion maps are in [Lafon1]_ , [Ferguson1]_, and
[Clementi1]_.

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
eigenvalues and eigenvectors.

   >>> u = MDAnalysis.Universe(PSF,DCD)
   >>> dist_matrix = DistanceMatrix(u)

We leave determination of the appropriate scale parameter epsilon to the user,
[Clementi1]_ uses a complex method involving the k-nearest-neighbors of a
trajectory frame, whereas others simple use a trial-and-error approach with
a constant epsilon. For those users, a
`~MDAnalysis.analysis.diffusionmap.EpsilonConstant` class has been provided.
Any user choosing to write their own Epsilon class must write a class that
inherits from `~MDAnalysis.analysis.diffusionmap.Epsilon`, and call the
`determine_epsilon` function required by the API before running the diffusion
map.

   >>> epsilon_matrix = EpsilonConstant(DistanceMatrix, 500)
   >>> epsilon_matrix.determine_epsilon()
   >>> dmap = diffusionmap.DiffusionMap(dist_matrix, epsilon_matrix)
   >>> dmap.decompose_kernel()
   >>> eigenvalues = dmap.eigenvalues
   >>> eigenvectors = dmap.eigenvectors

From here we can perform an embedding onto the k dominant eigenvectors. This
is similar to the idea of a transform in Principal Component Analysis, but the
non-linearity of the map means there is no explicit relationship between the
lower dimensional space and our original trajectory. However, this is an
isometry (distance preserving map), which means that points close in the lower
dimensional space are close in the higher-dimensional space and vice versa.
In order to embed into the most relevant low-dimensional space, there should
exist some number of dominant eigenvectors, whose corresponding eigenvalues
diminish at a constant rate until falling off, this is referred to as a
spectral gap and should be somewhat apparent for a system at equilibrium with a
high number of frames.

   >>> num_dominant_eigenvectors = # some number less than the number of frames
   >>> fit = dmap.transform(num_dominant_eigenvectors)

From here it can be difficult to interpret the data, and is left as a task
for the user. The `diffusion distance` between frames i and j is best
approximated by the euclidean distance  between rows i and j of self.fit.

Classes
-------

.. autoclass:: DiffusionMap
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
two coordinate sets (as implemented
in :func:`MDAnalysis.lib.qcprot.CalcRMSDRotationalMatrix`).
When using this module in published work please cite [Theobald2005]_.


"""
from six.moves import range
import logging

import numpy as np
from MDAnalysis.analysis import rms

from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.diffusionmap")


class DistanceMatrix(AnalysisBase):
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
        self.dist_matrix = np.zeros((self.nframes, self.nframes))
        self._i = -1

    def _single_frame(self):
        self._i = self._i + 1
        iframe = self._ts.frame
        i_ref = self.atoms.positions - self.atoms.center_of_mass()
        # diagonal entries need not be calculated due to metric(x,x) == 0 in
        # theory, _ts not updated properly. Possible savings by setting a
        # cutoff for significant decimal places to sparsify matrix
        for j, ts in enumerate(self._u.trajectory[self.start:self.stop:self.step]):
            self._ts = ts
            if self._ts.frame >= iframe:
                if self._i == j:
                    self.dist_matrix[self._i, self._i] = 0
                else:
                    j_ref = self.atoms.positions-self.atoms.center_of_mass()
                    dist = self._metric(i_ref, j_ref, weights=self._weights)
                    self.dist_matrix[self._i, j] = (dist if dist > self._cutoff
                                                    else 0)
                    self.dist_matrix[j, self._i] = self.dist_matrix[self._i, j]
        self._ts = self._u.trajectory[iframe]

    def save(self, filename):
        np.savetxt(filename, self.dist_matrix)
        logger.info("Wrote the distance-squared matrix to file %r", filename)


class Epsilon(object):
    """Manipulates a distance matrix by local scale parameters

    Attributes
    ----------
    scaledMatrix : DistanceMatrix object
        A matrix with each term divided by a local scale parameter

    Methods
    -------
    determineEpsilon()
        Determine local scale parameters using a chosen algorithm
    """
    def __init__(self, DistanceMatrix, **kwargs):
        self._DistanceMatrix = DistMatrix

    def determine_epsilon(self):
        pass


class EpsilonConstant(Epsilon):
    """Premade function for the determination of epsilon based on providing
    an arbitrary constant"""

    def __init__(self, DistanceMatrix, epsilon):
        """
        Parameters
        ----------
        epsilon : int
            The value of epsilon to be used as a local scale parameter
        """
        self.scaledMatrix = DistanceMatrix.dist_matrix
        self._epsilon = epsilon

    def determine_epsilon(self):
        self.scaledMatrix /= self._epsilon
        return


class DiffusionMap(AnalysisBase):
    """Non-linear dimension reduction method

    Dimension reduction with diffusion mapping of selected structures in a
    trajectory.

    Attributes
    ----------
    eigenvalues: array
        Eigenvalues of the diffusion map
    eigenvectors: array
        Eigenvectors of the diffusion map
    fit : array
        After calling `transform(num_eigenvectors)` the diffusion map embedding
        into the lower dimensional diffusion space will exist here.

    Methods
    -------
    decompose_kernel()
        Constructs an anisotropic diffusion kernel and performs eigenvalue
        decomposition on it.

    spectral_gap()
        Retrieve a guess for the set of eigenvectors reflecting the intrinsic
        dimensionality of the molecular system.

    transform(num_eigenvectors)
        Perform an embedding of a frame into the eigenvectors representing
        the collective coordinates.
    """

    def __init__(self, DistanceMatrix, epsilon, weights=None, timescale=1):
        """
        Parameters
        -------------
        DistanceMatrix : DistanceMatrix object
            Distance matrix to be made into a diffusion kernel and perform
            an eigenvector decomposition on.
        epsilon : Epsilon object
            Specifies the method used for the choice of scale parameter in the
            diffusion map. More information in [1], [2] and [3].
        weights: list, optional
            The list has to have the same length as the trajectory.
            With 'None' the weight of each frame of the trajectory will be the
            same.
        timescale: int, optional
            The number of steps in the random walk, large t reflects global
            structure whereas small t indicates local structure.
        """
        self.DistanceMatrix = DistanceMatrix
        self._nframes = DistanceMatrix.nframes
        self._t = timescale
        self._epsilon = epsilon

        if weights is None:
            # weights do not apply to metric but density of data
            self._weights_ker = np.ones((self._nframes,))
        else:
            if weights.shape[0] != self._nframes:
                raise ValueError("The weight should have the same length as "
                                 'the trajectory')
            else:
                # weights are constructed as relative to the mean
                self._weights_ker = (np.asarray(weights, dtype=np.float64) /
                                     np.mean(weights))

    def decompose_kernel(self):
        # this should be a reference to the same object as
        # self.DistanceMatrix.dist_matrix
        self._kernel = self._epsilon.scaledMatrix

        # take negative exponent of scaled matrix to create Isotropic kernel
        self._kernel = np.exp(-self._kernel)

        # define an anistropic diffusion term q
        q_vector = np.zeros((self._nframes, ))
        # weights should reflect the density of the points on the manifold
        for i in range(self._nframes):
            q_vector[i] = np.dot(self._kernel[i, :], self._weights_ker)

        # Form a new kernel from the anisotropic diffusion term q
        self._kernel /= np.sqrt(q_vector[:, np.newaxis].dot(q_vector[np.newaxis]))

        # Weighted Graph Laplacian normalization on this new graph
        d_vector = np.zeros((self._nframes, ))
        for i in range(self._nframes):
            d_vector[i] = np.dot(self._kernel[i, :], self._weights_ker)

        for i in range(self._nframes):
            self._kernel[i, :] = self._kernel[i, :] * self._weights_ker

        # Define anisotropic transition by dividing kernel by this term
        self._kernel /= np.sqrt(d_vector[:, np.newaxis].dot(d_vector[np.newaxis]))

        # Apply timescaling
        for i in range(self._t):
            if i > 1:
                self._kernel.__matmul__(self._kernel)

        eigenvals, eigenvectors = np.linalg.eig(self._kernel)

        eg_arg = np.argsort(eigenvals)
        self.eigenvalues = eigenvals[eg_arg[::-1]]
        self.eigenvectors = eigenvectors[eg_arg[::-1], :]

    def spectral_gap(self):
        """An ad hoc method to determine the spectral
            gap in a set of eigenvalues

        """
        pass

    def _prepare(self):
        self.fit = np.zeros((self.eigenvectors.shape[0],
                                  self.num_eigenvectors))

    def _single_frame(self):
        # The diffusion map embedding takes the ith sample in the
        # data matrix and maps it to each of the ith coordinates
        # in the set of k-dominant eigenvectors
        for k in range(self.num_eigenvectors):
            self.fit[self._ts.frame][k] = self.eigenvectors[k][self._ts.frame]

    def transform(self, num_eigenvectors):
        """ Embeds a trajectory via the diffusion map

        Parameter
        ---------
        num_eigenvectors : int
            The number of dominant eigenvectors to be used for
            diffusion mapping

        Return
        ------
        fit : array
            The diffusion map embedding as defined by [Ferguson1]_.
            This isn't a linear transformation, but an isometry
            between the higher dimensional space and the space spanned by
            the eigenvectors.
        """
        self.num_eigenvectors = num_eigenvectors
        self._setup_frames(self.DistanceMatrix._trajectory,
                           self.DistanceMatrix.start, self.DistanceMatrix.stop,
                           self.DistanceMatrix.step)
        self.run()
        return self.fit
