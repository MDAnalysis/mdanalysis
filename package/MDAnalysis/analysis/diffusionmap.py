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
the largest eigenvalues are the slowest collective coordinates. The complexity of
the diffusion map is O(N^3), where N is the number of frames in the trajectory.
Instead of a single trajectory a sample of protein structures can be used.
The sample should be equiblibrated, at least locally.  Different weights can be
used. The order of the sampled structures in the trajectory is irrelevant.

The :ref:`Diffusion-Map-tutorial` shows how to use diffusion map for dimension
reduction.

More details about diffusion maps are in [Lafon1]_  and [Clementi1]_.

.. _Diffusion-Map-tutorial:

Diffusion Map tutorial
--------------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF` and
:data:`~MDAnalysis.tests.datafiles.DCD`). Notice that this trajectory
does not represent a sample from the equilibrium. This violates a basic
assumption of the diffusion map and the results shouldn't be interpreted for
this reason/ This tutorial shows how to use the diffusionmap function.
First load all modules and test data ::

   >>> import MDAnalysis
   >>> import numpy as np
   >>> import MDAnalysis.analysis.diffusionmap as diffusionmap
   >>> from MDAnalysis.tests.datafiles import PSF,DCD

Given an universe or atom group, we can calculate the diffusion map from
that trajectory using :class:`DiffusionMap`:: and get the corresponding eigenvalues
and eigenvectors

   >>> u = MDAnalysis.Universe(PSF,DCD)
   >>> dmap = diffusionmap.DiffusionMap(u)
   >>> dmap.run()
   >>> eigenvalues = dmap.eigenvalues
   >>> eigenvectors = dmap.eigenvectors

Classes
-------

.. autoclass:: DiffusionMap
   :members:

References
---------

If you use this Dimension Reduction method in a publication, please
reference:
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
    d_matrix : array
        Array of all possible ij metric distances between frames in trajectory.
        This matrix is symmetric with zeros on the diagonal.
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
        self.d_matrix = np.zeros((self.nframes, self.nframes))

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
            self.d_matrix[traj_index, j] = dist if dist > self._cutoff else 0
            self.d_matrix[j, traj_index] = self.d_matrix[traj_index, j]

    def store(self, filename):
        pass

    def load(self, filename):
        pass


class DiffusionMap(DistMatrix):
    """Non-linear dimension reduction method

    Dimension reduction with diffusion map of the structures in the universe.

    Attributes
    ----------
    atoms : AtomGroup
        Selected atoms in trajectory subject to dimension reduction
    d_matrix : array
        Array of all possible ij metric distances between frames in trajectory
    eigenvalues: array
        Eigenvalues of the diffusion map
    eigenvectors: array
        Eigenvectors of the diffusion map
        The second and higher eigenvectors ev[i+1,:] represent the i-th slowest
        collective coordinates.
    """

    def __init__(self, DistMatrix, epsilon=epsilon_constant, weights=None,
                 **kwargs):
        """
        Parameters
        -------------
        DistMatrix : DistMatrix object
            Distance matrix to be made into a diffusion kernel and perform
            an Eigenvector decomposition on.
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
        self.d_matrix = DistMatrix.d_matrix
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
        # TODO: fix this, simply a proof of concept, need to add kwargs for
        # epsilon calculation
        self._epsilon = epsilon_constant(10, self.nframes)
        self._kernel2 = np.zeros((self.nframes, self.nframes))
        # possibly mappable
        for i in range(self.nframes):
            self._kernel2[i, :] = (np.exp((-self.d_matrix[i, :] ** 2) /
                                   (self._epsilon[i]*self._epsilon[:])))

        p_vector = np.zeros((self.nframes, ))
        d_vector = np.zeros((self.nframes, ))

        for i in range(self.nframes):
            p_vector[i] = np.dot(self._kernel2[i, :], self._weights_ker)

        self._kernel2 /= np.sqrt(p_vector[:, np.newaxis].dot(p_vector[np.newaxis]))

        for i in range(self.nframes):
            d_vector[i] = np.dot(self._kernel2[i, :], self._weights_ker)

        for i in range(self.nframes):
            self._kernel2[i, :] = self._kernel2[i, :] * self._weights_ker

        self._kernel2 /= np.sqrt(d_vector[:, np.newaxis].dot(d_vector[np.newaxis]))
        # eigenvalues and eigenvector are the collective coordinates
        eigenvals, eigenvectors = np.linalg.eig(self._kernel2)

        eg_arg = np.argsort(eigenvals)
        self.eigenvalues = eigenvals[eg_arg[::-1]]
        self.eigenvectors = eigenvectors[eg_arg[::-1], :]




def epsilon_average(dist_matrix, k, nframes):
    """Premade function for the determination of epsilon based on averaging after
    picking a neighbor based on connectivity.

    Parameters
    ----------
    dist_matrix : matrix
        The nframes-by-nframes matrix generated by DistMatrix
    k : int
        The index of the diffusion value in i_th matrix column used for epsilon
        value picking.
    nframes : int
        The number of frames used in the dimension reduction.
    """
    epsilon = np.zeros((nframes, ), )
    # modulus to prevent array index out of bounds exception
    k = k % self.nframes
    for i in range(nframes):
        # np.argsort(diffusion_matrix[i,:])#[10]]
        # picks k largest rmsd value in column for choice of epsilon
        epsilon[i] = d_matrix[i, np.argsort(d_matrix[i, :])[k]]

    epsilon = np.full((nframes, ), epsilon.mean())
    return epsilon


def epsilon_constant(epsilon, nframes):
    """Premade function for the determination of epsilon based on providing
    an arbitrary constant

    Parameters
    ----------
    epsilon : int
        The value of epsilon to be used as a local scale parameter
    nframes : int
        The number of frames used in the dimension reduction.
    """
    epsilon = np.full((nframes, ), epsilon)
    return epsilon
