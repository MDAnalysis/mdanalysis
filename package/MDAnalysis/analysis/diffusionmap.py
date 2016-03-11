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

:Author: Eugen Hruska
:Year: 2016
:Copyright: GNU Public License v3

The module contains the non-linear dimension reduction method diffusion map.
The diffusion map allows to get a quick estimate of the slowest collective 
coordinates for a trajectory. This non-linear dimension reduction method 
assumes that the trajectory is long enough to represents a probability 
distribution of as protein close to the equilibrium. Further the diffusion map 
assumes that the diffusion coefficients are constant. The eigenvectors with the 
largest eigenvalues are the slowest collective coordinates. The complexity of 
the diffusion map is O(N^2), where N is the number of frames in the trajectory. 
Instead of a single trajectory a sample of protein structures can be used. 
The sample should be equiblibrated, at least locally.  Different weights can be used. 
The order of the sampled structures in the trajectory is irrelevant.

The :ref:`Diffusion-Map-tutorial` shows how to use diffusion map for dimension reduction.

More details about diffusion maps are in [Lafon1]_  and [Clementi1]_.

.. _Diffusion-Map-tutorial:

Diffusion Map tutorial
--------------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF` and
:data:`~MDAnalysis.tests.datafiles.DCD`). Notice that these test data 
aren't representing a sample from the equilibrium. This violates a basic 
assumption of the diffusion map and the results shouldn't be interpreted for this reason. 
This tutorial shows how to use the diffusionmap function.
First load all modules and test data ::

   >>> import MDAnalysis
   >>> import numpy as np
   >>> import MDAnalysis.analysis.diffusionmap as diffusionmap
   >>> from MDAnalysis.tests.datafiles import PSF,DCD

In the simplest case, we can simply calculate the diffusion map from 
one trajectory :func:`diffusionmap`::

   >>> u = MDAnalysis.Universe(PSF,DCD)
   >>> eg,ev=diffusionmap.diffusionmap(u)

To see how the two slowest collective coordinates how the Other stuff in paper

   >>> import matplotlib.pyplot as plt
   >>> plt.scatter(ev[:,1],ev[:,2])
   >>> plt.title('diffusion map')
   >>> plt.xlabel("slowest collective coordinate")
   >>> plt.ylabel("second slowest collective coordinate")
   >>> plt.show()


Functions
---------

.. autofunction:: diffusionmap

References
---------

If you use this QCP rotation calculation method in a publication, please
reference:
..[Lafon1] Coifman, Ronald R., Lafon, Stephane (2006) Diffusion maps. 
Appl. Comput. Harmon. Anal. 21, 5–30.
..[Clementi1] Rohrdanz, M. A, Zheng, W, Maggioni, M, & Clementi, C. (2013) 
Determination of reaction coordinates via locally scaled 
diffusion map. Journal of Chemical Physics.



"""

import logging
import numpy as np
import MDAnalysis.lib.qcprot as qcp
from six.moves import range


def diffusionmap(u, select='all', epsilon='average', k=10, weight=None):
    """Non-linear dimension reduction method diffusion map

    diffusionmap(u, select, epsilon, k, weight)

    Dimension reduction with diffusion map of the structures in the universe.

    Parameters
    -------------
      *u*
         trajectory :class:`~MDAnalysis.core.AtomGroup.Universe`
         The trajectory can be a long trajectory 
      select: str, optional 
         1. any valid selection string for
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` 
         This selection of atoms is used to calculate the RMSD between different frames. Water should be excluded.
      epsilon : float, optional 
          Specifies the epsilon used for the diffusion map. More information in [1] and [2] 
          With 'average' the average of the RMSD to the k-nearest-neighbor will be used.  
      k : int, optional
          specifies the k for the k-nearest-neighbor is average epsilon is used.
      weight: numpy array, optional
         The numpy array has to have the same length as the trajectory.
         With 'None' the weight of each frame of the trajectory will be the same.
         If order of the weights has to be the same as the order of framesin the trajectory.


    Returns
    ----------------
    eigenvalues: numpy array
         Eigenvalues of the diffusion map
    eigenvectors: numpy array
         Eigenvectors of the diffusion map
         The second and higher eigenvectors ev[i+1,:] represent the i-th slowest collective coordinates.

    Notes
    ---------------
    
    The dimension reduction works in the following way:

    1. A RMSD between each every pair of frames is calculated.
    2. The normalized kernel is obtain as in [1] and [2].
    3. The eigenvalues and eigenvectors of the normalized kernel are the output.

    """
    logger = logging.getLogger('MDAnalysis.analysis.diffusionmap')
    frames = u.trajectory
    ref_atoms = u.select_atoms(select)
    nframes = len(frames)
    natoms = ref_atoms.n_atoms

    rmsd_matrix = np.zeros((nframes, nframes))
    kernel2 = np.zeros((nframes, nframes))

    if epsilon == 'average':
        epsilon = np.zeros((nframes, ), )
        type_epsilon = 'average'
    else:
        value_epsilon = epsilon
        epsilon = np.full((nframes, ), value_epsilon)
        type_epsilon = 'constant'


    rot = np.zeros(9)
    if weight is None:
        weights = np.full((nframes, ), 1)
    else:
        if weight.shape[0] != nframes:
            raise ValueError("The weight should have the same length as the trajectroy")
        else:
            weights = weight

    for i in range(nframes):
        logger.info("calculating rmsd from structure {0} to all".format(i))
        i_ref = np.copy(u.trajectory[i].positions-ref_atoms.center_of_mass())
        for j in range(i, nframes):
            j_ref = np.copy(u.trajectory[j].positions-ref_atoms.center_of_mass())
            rmsd_matrix[i, j] = qcp.CalcRMSDRotationalMatrix(i_ref.T.astype(np.float64), \
              j_ref.T.astype(np.float64), natoms, rot, weight)

    #fill in symmetric values
    rmsd_matrix = rmsd_matrix + rmsd_matrix.T - np.diag(rmsd_matrix.diagonal())

    #calculate epsilons 
    if type_epsilon == 'average':
        for i in range(nframes):
           #np.argsort(rmsd_matrix[i,:])#[10]]
            epsilon[i] = rmsd_matrix[i, np.argsort(rmsd_matrix[i, :])[k]]
        epsilon = np.full((nframes, ), epsilon.mean())


    logger.info('epsilon: {0}'.format(epsilon))

    #calculate normalized kernel
    for i in range(nframes):
        kernel2[i, :] = np.exp(-rmsd_matrix[i, :]**2/(epsilon[i]*epsilon[:]))

    p_vector = np.zeros((nframes, ))
    d_vector = np.zeros((nframes, ))
    for i in range(nframes):
        p_vector[i] = np.dot(kernel2[i, :], weights)

    kernel2 /= np.sqrt(p_vector[:, np.newaxis].dot(p_vector[np.newaxis]))

    for i in range(nframes):
        d_vector[i] = np.dot(kernel2[i, :], weights)

    for i in range(nframes):
        kernel2[i, :] = kernel2[i, :]*weights

    kernel2 /= np.sqrt(d_vector[:, np.newaxis].dot(d_vector[np.newaxis]))

    #eigenvalues and eigenvector are the collective coordinates
    eg, ev = np.linalg.eig(kernel2)

    eg_arg = np.argsort(eg)
    eg = eg[eg_arg[::-1]]
    ev = ev[eg_arg[::-1],:]

    return eg, ev
