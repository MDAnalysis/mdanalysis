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

More details how the correst slowest collective coordinates can be calculated.in:
[1] Nadler, B, Lafon, S, Coifman, R. R, & Kevrekidis, I. G. (2013) Diffusion maps, 
spectral clustering and reaction coordinates of dynamical systems. Applied and 
Computational Harmonic Analysis 21, 113–127.
[2] Rohrdanz, M. A, Zheng, W, Maggioni, M, & Clementi, C. (2013) Determi- nation 
of reaction coordinates via locally scaled diffusion map. Journal of Chemical Physics.

The :ref:`Diffusion-Map-tutorial` shows how to use diffusion map for dimension reduction.


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
   >>> import MDAnalysis.lib.qcprot as qcp
   >>> import MDAnalysis.analysis.diffusionmap as diffusionmap
   >>> from MDAnalysis.tests.datafiles import PSF,DCD

In the simplest case, we can simply calculate the diffusion map from one trajectory :func:`diffusionmap`::

   >>> u = MDAnalysis.Universe(PSF,DCD)
   >>> eg,ev=diffusionmap.diffusionmap(u)

To see how the two slowest collective coordinates how the Other stuff in paper

   >>> import matplotlib.pyplot as plt
   >>> import matplotlib.pyplot as plt
   >>> plt.scatter(ev[:,1],ev[:,2])
   >>> plt.title('diffusion map')
   >>> plt.xlabel("slowest collective coordinate")
   >>> plt.ylabel("second slowest collective coordinate")
   >>> plt.show()


Common usage
------------

To **fit a single structure** with :func:`alignto`::

   >>> ref = Universe(PSF, PDB_small)
   >>> mobile = Universe(PSF, DCD)     # we use the first frame
   >>> alignto(mobile, ref, select="protein and name CA", mass_weighted=True)

This will change *all* coordinates in *mobile* so that the protein
C-alpha atoms are optimally superimposed (translation and rotation).


Functions
---------

.. autofunction:: diffusionmap


"""

import os.path
from six.moves import range, zip, zip_longest
import numpy as np
import warnings
import logging


import MDAnalysis.lib.qcprot as qcp
from MDAnalysis.exceptions import SelectionError, SelectionWarning
from MDAnalysis.lib.log import ProgressMeter
import MDAnalysis.analysis.rms as rms


logger = logging.getLogger('MDAnalysis.analysis.diffusionmap')


def diffusionmap(u, select='all', epsilon='average', k=10, weight=None):
    """ Non-linear dimension reduction method diffusion map
    
    The dimension reduction is done in the following way:

    1. A RMSD between each every pair of frames is calculated.
    2. The normalized kernel is obtain as in [1] and [2].
    3. The eigenvalues and eigenvectors of the normalized kernel are the output.

    :Arguments:
      *u*
         trajectory :class:`~MDAnalysis.core.AtomGroup.Universe`
         The trajectory can be a long trajectory 
      *select*
         1. any valid selection string for
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` 
         This selection of atoms is used to calculate the RMSD between different frames. Water should be excluded.
      *epsilon* : float 
          Specifies the epsilon used for the diffusion map. More information in [1] and [2] 
          With 'average' the average of the RMSD to the k-nearest-neighbor will be used.  
      *k* : specifies the k for the k-nearest-neighbor is average epsilon is used.
      *weight* : None or numpy array of the same length as the trajectory
         With 'None' the weight of each frame of the trajectory will be the same.
         If order of the weights has to be the same as the order of framesin the trajectory.


    :Returns: eigenvalues, eigenvectors 
    the second and higher eigenvectors ev[i+1,:] represent the i-th slowest collective coordinates.

    """

    frames=u.trajectory
    ref_atoms = u.select_atoms(select)
    nframes=len(frames)
    natoms = ref_atoms.n_atoms

    rmsd_matrix=np.zeros((nframes,nframes), dtype=np.float64)
    epsilon=np.zeros((nframes,), dtype=np.float64)
    kernel2=np.zeros((nframes,nframes), dtype=np.float64)

    rmsd = np.zeros((nframes,))
    rot= np.zeros(9, dtype=np.float64)
    weight=None
    type_epsilon='average'
    if weight is None:
      weights=np.full((nframes,),1)
    else:
      if weight.shape[0]!=nframes:
        raise ValueError("The weight should have the same length as the trajectroy")
      else:
        weights=weight

    for i in range(nframes):
       logger.info("calculating rmsd from structure "+str(i)+" to all")
       i_ref=np.copy(u.trajectory[i].positions-ref_atoms.center_of_mass())
       for j in range(nframes):
         j_ref=np.copy(u.trajectory[j].positions-ref_atoms.center_of_mass())
         rmsd_matrix[i,j] = qcp.CalcRMSDRotationalMatrix(i_ref.T.astype(np.float64), j_ref.T.astype(np.float64), natoms, rot, weight)

    #calculate epsilons 
    for i in range(nframes):
      #np.argsort(rmsd_matrix[i,:])#[10]]
      epsilon[i]=rmsd_matrix[i,np.argsort(rmsd_matrix[i,:])[k]]

    if epsilon=='average':
      epsilon=np.full((nframes,),epsilon.mean())
    else:
      epsilon=np.full((nframes,),epsilon)

    logger.info("epsilon: "+str(epsilon))

    #calculate normalized kernel
    for i in range(nframes):
      kernel2[i,:]=np.exp(-rmsd_matrix[i,:]**2/(epsilon[i]*epsilon[:]))

    p_vector = np.zeros((nframes,))
    d_vector = np.zeros((nframes,))
    for i in range(nframes):
      p_vector[i]=np.dot(kernel2[i,:],weights)

    kernel2 /= np.sqrt(p_vector[:,np.newaxis].dot(p_vector[np.newaxis]))

    for i in range(nframes):
      d_vector[i]=np.dot(kernel2[i,:],weights)

    for i in range(nframes):
      kernel2[i,:]= kernel2[i,:]*weights

    kernel2 /= np.sqrt(d_vector[:,np.newaxis].dot(d_vector[np.newaxis]))

    #eigenvalues and eigenvector are the collective coordinates
    eg, ev=np.linalg.eig(kernel2)

    if np.abs(eg[0]-1)>1.0e-14:
      raise ValueError("Lowest eigenvalue should be 1 up to numeric precision")

    return eg, ev
