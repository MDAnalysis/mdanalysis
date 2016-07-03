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

We leave determination of the appropriate scale parameter epsilon to the user,
[Clementi1]_ uses a complex method involving the k-nearest-neighbors of a
trajectory frame, whereas others simple use a trial-and-error approach with
a constant epsilon. Users wishing to use a complex method to determine epsilon
will have to initialize a distance matrix and scale it appropriately themselves
and then pass the new scaled distance matrix as a parameter.

   >>> dmap = diffusionmap.DiffusionMap(dist_matrix, epsilo =2)
   >>> dmap.run()
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

   >>> num_eigenvectors = # some number less than the number of frames
   >>> fit = dmap.transform(num_eigenvectors)

From here it can be difficult to interpret the data, and is left as a task
for the user. The `diffusion distance` between frames i and j is best
approximated by the euclidean distance  between rows i and j of
self.diffusion_space. A Jupyter notebook providing an analysis of protein
opening and closing is provided [here](#TODO url to notebook)

Classes
-------

.. autoclass:: DiffusionMap
.. autoclass:: DistMatrix

References
---------

If you use this Dimension Reduction method in a publication, please
reference:
..[Ferguson1] Ferguson, A. L.; Panagiotopoulos, A. Z.; Kevrekidis, I. G.
Debenedetti,  P. G. Nonlinear dimensionality reduction in molecular
simulation: The diffusion map approach  Chem. Phys. Lett. 509, 1−11 (2011)
..[Lafon1] Coifman, Ronald R., Lafon, Stephane Diffusion maps.
Appl. Comput. Harmon. Anal. 21, 5–30  (2006).
..[Lafon2] Boaz Nadler, Stéphane Lafon, Ronald R. Coifman, Ioannis G. Kevrekidis.
Diffusion maps, spectral clustering and reaction coordinates
of dynamical systems. Appl. Comput. Harmon. Anal. 21 (2006) 113–127
..[Clementi1] Rohrdanz, M. A, Zheng, W, Maggioni, M, & Clementi, C.
Determination of reaction coordinates via locally scaled
diffusion map. J. Chem. Phys. 134, 124116 (2011).

.. If you choose the default metric, this module uses the fast QCP algorithm
[Theobald2005]_ to calculate the root mean square distance (RMSD) between
two coordinate sets (as implemented
in :func:`MDAnalysis.lib.qcprot.CalcRMSDRotationalMatrix`).
When using this module in published work please cite [Theobald2005]_.


"""
from six.moves import range
import logging
import warnings

import MDAnalysis as mda
import numpy as np

from .rms import rmsd
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
    def __init__(self, u, select='all', metric=rmsd, cutoff=1E0-5,
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
        self._calculated = False
        # remember that this must be called before referencing self.nframes
        self._setup_frames(traj, start, stop, step)

    def _prepare(self):
        self.dist_matrix = np.zeros((self.nframes, self.nframes))

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

    def save(self, filename):
        np.save(filename, self.dist_matrix)
        logger.info("Wrote the distance-squared matrix to file %r", filename)


class DiffusionMap(object):
    """Non-linear dimension reduction method

    Dimension reduction with diffusion mapping of selected structures in a
    trajectory.

    Attributes
    ----------
    eigenvalues: array
        Eigenvalues of the diffusion map
    eigenvectors: array
        Eigenvectors of the diffusion map
    diffusion_space : array
        After calling `transform(n_eigenvectors)` the diffusion map embedding
        into the lower dimensional diffusion space will exist here.

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
            diffusion map. More information in [1], [2] and [3], Default: 1.
        timescale: int, optional
            The number of steps in the random walk, large t reflects global
            structure whereas small t indicates local structure, Default: 1
            """
        if isinstance(u, mda.Universe):
            self._dist_matrix = DistanceMatrix(u, **kwargs)
        elif isinstance(u, DistanceMatrix):
            self._dist_matrix = u
        else:
            raise ValueError("U is not a Universe or DistanceMatrix and"
                             " so the DiffusionMap has no data to work with.")
        self._epsilon = epsilon
        # important for transform function and length of .run() method
        self._nframes = self._dist_matrix.nframes
        if self._nframes > 2000:
            warnings.warn("The distance matrix is very large, and can "
                          "be very slow to compute. Consider picking a larger "
                          "step size in distance matrix initialization.")


    def run(self):
        # run only if distance matrix not already calculated
        if not self._dist_matrix._calculated:
            self._dist_matrix.run()
        self._scaled_matrix = (self._dist_matrix.dist_matrix ** 2 /
                               self._epsilon)
        # take negative exponent of scaled matrix to create Isotropic kernel
        self._kernel = np.exp(-self._scaled_matrix)
        D_inv = np.diag(1 / self._kernel.sum(1))
        self._diff = np.dot(D_inv, self._kernel)
        eigenvals, eigenvectors = np.linalg.eig(self._diff)
        sort_idx = np.argsort(eigenvals)[::-1]
        self.eigenvalues = eigenvals[sort_idx]
        self.eigenvectors = eigenvectors[sort_idx]
        self._calculated = True

    def transform(self, n_eigenvectors, time):
        """ Embeds a trajectory via the diffusion map

        Parameter
        ---------
        n_eigenvectors : int
            The number of dominant eigenvectors to be used for
            diffusion mapping
        time : float
            Exponent that eigenvalues are raised to for embedding, for large
            values, more dominant eigenvectors determine diffusion distance.
        Return
        ------
        diffusion_space : array
            The diffusion map embedding as defined by [Ferguson1]_.
        """
        return (self.eigenvectors[1:n_eigenvectors+1,].T *
                (self.eigenvalues[1:n_eigenvectors+1]**time))
