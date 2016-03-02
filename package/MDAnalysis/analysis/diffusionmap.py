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

BETTER
The module contains functions to fit a target structure to a reference
structure. They use the fast QCP algorithm to calculate the root mean
square distance (RMSD) between two coordinate sets [Theobald2005]_ and
the rotation matrix *R* that minimizes the RMSD [Liu2010]_. (Please
cite these references when using this module.).

Typically, one selects a group of atoms (such as the C-alphas),
calculates the RMSD and transformation matrix, and applys the
transformation to the current frame of a trajectory to obtain the
rotated structure. The :func:`alignto` and :func:`rms_fit_trj`
functions can be used to do this for individual frames and
trajectories respectively.

The :ref:`Diffusion-Map-tutorial` shows how to use diffusion map for dimension reduction.


.. _Diffusion-Map-tutorial:

Diffusion Map tutorial
--------------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF` and
:data:`~MDAnalysis.tests.datafiles.DCD`). First load all modules and test data ::

   >>> import MDAnalysis
   >>> import numpy as np
   >>> import MDAnalysis.lib.qcprot as qcp
   >>> import MDAnalysis.analysis.diffusionmap as diffusionmap
   >>> from MDAnalysis.tests.datafiles import PSF,DCD


In the simplest case, we can simply calculate the diffusion map from one trajectory :func:`diffusionmap`::

   >>> u = MDAnalysis.Universe(PSF,DCD)
   >>> eg,ev=diffusionmap.diffusionmap(u)

Since o(N^2) long trajectories take too long, increase stride

weight

plot

intro

what eg, ev do
eg should after several steps fall off to 0
ev[i,:] is the ith slowest eigenvector, this can be used to dimension reduction

Other stuff in paper


   >>> import matplotlib.pyplot as plt
   >>> plt.plot(eg[:10])
   >>> plt.title('diffusion map eigenvalues')
   >>> plt.xlabel("eigenvalue")
   >>> plt.ylabel("index of eigenvalue")
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


def diffusionmap(u, distance='rmsd', type_epsilon='average', stride=1):#, select="all", mass_weighted=False,subselection=None, tol_mass=0.1, strict=False):
    """Spatially align *mobile* to *reference* by doing a RMSD fit on *select* atoms.

    The superposition is done in the following way:

    1. A rotation matrix is computed that minimizes the RMSD between
       the coordinates of `mobile.select_atoms(sel1)` and
       `reference.select_atoms(sel2)`; before the rotation, *mobile* is
       translated so that its center of geometry (or center of mass)
       coincides with the one of *reference*. (See below for explanation of
       how *sel1* and *sel2* are derived from *select*.)

    2. All atoms in :class:`~MDAnalysis.core.AtomGroup.Universe` that
       contains *mobile* are shifted and rotated. (See below for how
       to change this behavior through the *subselection* keyword.)

    The *mobile* and *reference* atom groups can be constructed so that they
    already match atom by atom. In this case, *select* should be set to "all"
    (or ``None``) so that no further selections are applied to *mobile* and
    *reference*, therefore preserving the exact atom ordering (see
    :ref:`ordered-selections-label`).

    .. Warning:: The atom order for *mobile* and *reference* is *only*
       preserved when *select* is either "all" or ``None``. In any other case,
       a new selection will be made that will sort the resulting AtomGroup by
       index and therefore destroy the correspondence between the two groups. **It
       is safest not to mix ordered AtomGroups with selection strings.**

    :Arguments:
      *mobile*
         structure to be aligned, a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
         or a whole :class:`~MDAnalysis.core.AtomGroup.Universe`
      *reference*
         reference structure, a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
         or a whole :class:`~MDAnalysis.core.AtomGroup.Universe`
      *select*
         1. any valid selection string for
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` that produces identical
            selections in *mobile* and *reference*; or
         2. dictionary ``{'mobile':sel1, 'reference':sel2}``.
            (the :func:`fasta2select` function returns such a
            dictionary based on a ClustalW_ or STAMP_ sequence alignment); or
         3.  tuple ``(sel1, sel2)``

         When using 2. or 3. with *sel1* and *sel2* then these selections can also each be
         a list of selection strings (to generate a AtomGroup with defined atom order as
         described under :ref:`ordered-selections-label`).
      *mass_weighted* : boolean
         ``True`` uses the masses :meth:`reference.masses` as weights for the
         RMSD fit.
      *tol_mass*
         Reject match if the atomic masses for matched atoms differ by more than
         *tol_mass* [0.1]
      *strict*
         ``True``
             Will raise :exc:`SelectioError` if a single atom does not
             match between the two selections.
         ``False`` [default]
             Will try to prepare a matching selection by dropping
             residues with non-matching atoms. See :func:`get_matching_atoms`
             for details.
      *subselection*
         Apply the transformation only to this selection.

         ``None`` [default]
             Apply to `mobile.universe.atoms` (i.e. all atoms in the
             context of the selection from *mobile* such as the rest of a
             protein, ligands and the surrounding water)
         *selection-string*
             Apply to `mobile.select_atoms(selection-string)`
         :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
             Apply to the arbitrary group of atoms

    :Returns: eigenvalues, eigenvectors 

    """

    frames=u.trajectory
    ref_atoms = u.select_atoms('all')
    nframes=len(frames)
    natoms = ref_atoms.n_atoms

        #define data

    #mobile_atoms = u.atoms

    rmsd_matrix=np.zeros((nframes,nframes), dtype=np.float64)
    epsilon=np.zeros((nframes,), dtype=np.float64)
    kernel2=np.zeros((nframes,nframes), dtype=np.float64)

    rmsd = np.zeros((nframes,))
    rot= np.zeros(9, dtype=np.float64)
    weight=None
    type_epsilon='average'
    if weight==None:
      weights=np.full((nframes,),1)


    for i in range(nframes):
       logger.info("calculating rmsd from structure "+str(i)+" to all")
       i_ref=np.copy(u.trajectory[i].positions-ref_atoms.center_of_mass())
       #i_ref, ref_atoms.center_of_mass()#=i_ref-ref_atoms.center_of_mass()
       #ref_atoms.center_of_mass()
       for j in range(nframes):
         j_ref=np.copy(u.trajectory[j].positions-ref_atoms.center_of_mass())
         rmsd_matrix[i,j] = qcp.CalcRMSDRotationalMatrix(i_ref.T.astype(np.float64), j_ref.T.astype(np.float64), natoms, rot, weight)
         #print i,j,rmsd_matrix[i,j]

    #print rmsd_matrix[:10,:10]  
    #print rmsd_matrix.max() 
       #all_traj=u.trajectory[0]
       #current_ref=u.trajectory[-1]
       #rmsd_matrix[i,:] = rms.rmsd(u2.atoms.coordinates(),u.atoms.coordinates())

    #calculate epsilons (average)
    #rmsd_matrix[0,np.argsort(rmsd_matrix[0,:])[10]]
    for i in range(nframes):
      #np.argsort(rmsd_matrix[i,:])#[10]]
      epsilon[i]=rmsd_matrix[i,np.argsort(rmsd_matrix[i,:])[10]]

    if type_epsilon=='average':
      epsilon=np.full((nframes,),epsilon.mean())

    logger.info("epsilon: "+str(epsilon))

    #calculate exponetial to get 
    for i in range(nframes):
      kernel2[i,:]=np.exp(-rmsd_matrix[i,:]**2/(epsilon[i]*epsilon[:]))

    #normalisation
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

    eg, ev=np.linalg.eig(kernel2)

    if np.abs(eg[0]-1)>1.0e-14:
      raise ValueError("Lowest eigenvalue should be 1")

    return eg, ev
