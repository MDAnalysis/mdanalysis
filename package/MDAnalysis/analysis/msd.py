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

r"""
Mean Squared Displacement --- :mod:`MDAnalysis.analysis.msd`
==============================================================

This module implements the calculation of Mean Squared Displacmements (MSDs).
MSDs can be used to characterise the speed at which particles move and has its roots
in the study of Brownian motion. For a full explanation of the best practices for the computation of MSDs and the subsequent calculation of self diffusivities the reader is directed to [Maginn2019]_.
MSDs are computed from the following expression:

.. math::

   MSD(r_{d}) = \bigg{\langle} \frac{1}{N} \sum_{i=1}^{N} |r_{d}  - r_{d}(t_0)|^2 \bigg{\rangle}_{t_{0}}

Where :math:`N` is the number of equivalent particles the MSD is calculated over, :math:`r` are their coordinates and :math:`d` the desired
dimensionality of the MSD. Note that while the definition of the MSD is universal, there are many practical considerations to computing the MSD
that vary between implementations. In this module, we compute a "windowed" MSD, where the MSD is averaged over all possible lag times :math:`t \le t_{max}`,
where :math:`t_{max}` is the length of the trajectory, thereby maximising the number of samples.

The computation of the MSD in this way can be computationally intensive due to it's :math:`N^2` scaling with respect to :math:`t_{max}`. 
An algorithm to compute the MSD with :math:`N log(N)` scaling based on a Fast Fourier Transform is known and can be accessed by setting fft=True [Calandri2011]_.
The python implementation was originally presented here [SO2015]_. 

Computing an MSD
----------------
This example computes a 2D MSD for the movement of phospholipid headgroups in the xy plane.
This is a common starting point for computing self diffusivity of phospholipids in a bilayer.
Files provided as part of the MDAnalysis test suite are used
(in the variables :data:`~MDAnalysis.tests.datafiles.GRO_MEMPROT` and
:data:`~MDAnalysis.tests.datafiles.XTC_MEMPROT`)

First load all modules and test data

    >>> import MDAnalysis as mda
    >>> import MDAnalysis.analysis.msd as msd
    >>> from MDAnalysis.tests.datafiles import GRO_MEMPROT, XTC_MEMPROT

Given a universe containing trajectory data we can extract the MSD
Analyis by using the class :class:`MeanSquaredDisplacement` 

    >>> u = mda.Universe(GRO_MEMPROT, XTC_MEMPROT)
    >>> MSD = msd.MeanSquaredDisplacement(u, 'name PO4', msd_type='xy', fft=True)
    >>> MSD.run()

The MSD can then be accessed as

    >>> msd =  MSD.timeseries

Visual inspection of the MSD is important, so lets take a look at it with a simple plot.

    >>> import matplotlib.pyplot as plt
    >>> nframes = len(u.trajectory)
    >>> timestep = 2 # this needs to be the time between frames in you trajectory
    >>> lagtimes = np.arange(nframes)*timestep # make the lag time axis
    >>> plt.plot(msd, lagtimes)
    >>> plt.show()

We can see that the MSD is roughly linear between a lag-time of 500 and 1000. Linearity of a segment of the MSD is required to accurately determine self diffusivity.
This linear segment represents the so called "middle" of the MSD plot, where ballistic trajectories at short time-lags are excluded along with poorly averaged data at long time-lags.
This can be confirmed with a log-log plot as is often reccomended [Maginn2019]_ where the "middle" segment can be identified as having a slope of 1.
Now that we have identified what segment of our MSD to analyse, lets compute a self diffusivity.

Computing Self Diffusivity
--------------------------------
Self diffusivity is closely related to the MSD.

.. math::

   D_d = \frac{1}{2d} \lim_{t \to \infty} \frac{d}{dt} MSD(r_{d}) 

From the MSD, self diffusivities :math:`D` with the desired dimensionality :math:`d` can be computed by fitting the MSD with respect to the lag time to a linear model. 
An example of this is shown below.

    >>> import scipy.stats.linregress as lr
    >>> start_time = 500
    >>> start_index = start_time/timestep
    >>> end_time = 1000
    >>> end_index = end_time/timestep
    >>> linear_model = lr(lagtimes[start_index:end_index], msd[start_index:end_index])
    >>> slope = linear_model.slope
    >>> error = linear_model.r-value
    >>> D = slope * 1/(2*MSD._dim_fac) #dim_fac is 2 as we computed a 2D msd ('xy')

We have now computed a self diffusivity!


Notes
_____
There are several factors that must be taken into account when setting up and processing trajectories for computation of self diffusivities.
These include specific instructions around simulation settings, using unwrapped trajectories and maintaining relativley small elapsed time between saved frames. 
Additionally corrections for finite size effects are sometimes employed along with varied means of estimating errors.
The reader is directed to the following review, which describes many of the common pitfalls [Maginn2019]_. There are other ways to compute self diffusivity including from a Green-Kubo integral. At this point in time these methods are beyond the scope of this module



References
----------

.. [Maginn2019] Maginn, E. J.; Messerly, R. A.; Carlson, D. J.; Roe, D. R.; Elliott, J. R. Best Practices for Computing Transport Properties 1. Self-Diffusivity and Viscosity from Equilibrium Molecular Dynamics [Article v1.0]. Living J. Comput. Mol. Sci. 2019, 1 (1).
.. [Calandri2011] Calandrini, V.; Pellegrini, E.; Calligari, P.; Hinsen, K.; Kneller, G. R. NMoldyn-Interfacing Spectroscopic Experiments, Molecular Dynamics Simulations and Models for Time Correlation Functions. Collect. SFN 2011, 12, 201–232.
.. [SO2015] https://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft

Classes and Functions
---------------------

.. autoclass:: MeanSquaredDisplacement
    :members:

"""

from __future__ import division, absolute_import
from six.moves import zip

import os
import errno
import warnings
import bz2
import functools

import numpy as np
from numpy.fft import fft,ifft
import logging
import MDAnalysis
import tidynamics



class MeanSquaredDisplacement(object):
    r"""Class representing Mean Squared Displacement

    Parameters
    ----------
    u : Universe
        An MDAnalysis :class:`Universe`
    selection : str
        An MDAnalysis selection string
    msd_type : str
        The dimensions to use for the msd calculation:
        one of ('xyz', 'xy', 'yz', 'xz', 'x', 'y', 'z')
    fft : bool
        Use a fast FFT based algorithm for computation of the MSD

    Returns
    -------
    timeseries : :class:`np.ndarray`
        the MSD as a function of lag time
    
    """
    # Attributes
    # ----------
    # _dim_fac : int
    #     dimensionality of the MSD
    # _dim : arraylike
    #     array used to slice the trajectory along the xyz axis to acheive the correct dimensionality
    # _position_array : :class:`np.ndarray`
    #     positions used to calculate the MSD
    # N_particles : int
    #     number of particles 


    def __init__(self, u, selection, msd_type='xyz', fft=True):

        #args
        self.u = u
        self.selection = selection
        self.msd_type = msd_type
        self.fft = fft
        
        #local
        self.dim_fac = 0
        self._dim = None
        self.atoms = None
        self.n_frames = len(self.u.trajectory)
        self._position_array = None

        #result
        self.timeseries = None
        #prep
        self._prepare()

            
    def _prepare(self):
        self.parse_msd_type()
        self.select_reference_positions()
        
    def parse_msd_type(self):
        
        if self.msd_type == 'xyz': # full 3d
            self._dim = [0,1,2]
            self.dim_fac = 3.0
            
        elif self.msd_type == 'xy': # xy
            self._dim = [0,1]
            self.dim_fac = 2.0

        elif self.msd_type == 'xz': # xz
            self._dim = [0,2]
            self.dim_fac = 2.0

        elif self.msd_type == 'yz': # yz
            self._dim = [1,2]
            self.dim_fac = 2.0

        elif self.msd_type == 'x': # x
            self._dim = [0]
            self.dim_fac = 1.0

        elif self.msd_type == 'y': # y
            self._dim = [1]
            self.dim_fac = 1.0

        elif self.msd_type == 'z': # z
            self._dim = [2]
            self.dim_fac = 1.0

        else:
            raise ValueError('invalid msd_type specified')

    def select_reference_positions(self):
        self._position_array = self.u.trajectory.timeseries(self.u.select_atoms(self.selection),order='fac') 
        self.N_particles = self._position_array.shape[1] 
    
    def run(self):
        if self.fft == True:
            self.timeseries = self._run_fft()
        else:
            self.timeseries = self._run_simple()

    def _run_simple(self): # naieve algorithm without FFT
        msds_byparticle = np.zeros([self.n_frames, self.N_particles])
        lagtimes = np.arange(1,self.n_frames)
        msds_byparticle[0,:] = np.zeros(self.N_particles) # preset the zero lagtime so we dont have to iterate through
        for n in range(self.N_particles):
            for lag in lagtimes:
                disp = self._position_array[:-lag,n,self._dim if lag else None] - self._position_array[lag:,n,self._dim]
                sqdist = np.square(disp, dtype=np.float64).sum(axis=1, dtype=np.float64) #accumulation in anything other than f64 is innacurate
                msds_byparticle[lag,n] = np.mean(sqdist, dtype=np.float64)
        msds = msds_byparticle.mean(axis=1, dtype=np.float64)
        return msds

    def _run_fft(self): #with FFT, np.f64 bit prescision required.
        msds_byparticle = []
        reshape_positions = self._position_array[:,:,self._dim].astype(np.float64)
        N=reshape_positions.shape[0]
        D=np.square(reshape_positions).sum(axis=2, dtype=np.float64) 
        D=np.append(D,np.zeros(reshape_positions.shape[:2]), axis=0) 
        Q=2*D.sum(axis=0, dtype=np.float64)
        S1=np.zeros(reshape_positions.shape[:2],dtype=np.float64)
        for m in range(N):
            Q=Q-D[m-1,:]-D[N-m,:]
            S1[m,:]=Q/(N-m)
        S2accumulate = []
        for i in range(reshape_positions.shape[2]):
            S2accumulate.append(self.autocorrFFT(reshape_positions[:,:,i]))
        S2= np.sum(S2accumulate,axis=0,dtype=np.float64)

        msds_byparticle.append(S1-2*S2)
        msds = np.concatenate(msds_byparticle,axis=1).mean(axis=-1, dtype=np.float64)
        return msds

    @staticmethod
    def autocorrFFT(x):
        r""" Calculates an autocorrelation function of the input signal via an FFT

        Parameters
        ----------
        x : :class:`np.ndarray`
            array to compute the autocorrelation for

        Returns
        -------
        autocorr : :class:`np.ndarray`
            the autocorrelation 
        """
        N=(x.shape[0])

        #find closest power of 2 (code adapted from tidynamics)
        current_exp = int(np.ceil(np.log2(N+1)))
        if N == 2**current_exp:
            n_fft = N
        if N < 2**current_exp:
            n_fft = 2**current_exp
        elif N > 2**current_exp:
            n_fft = 2**(current_exp+1)

        F = fft(x, n=2*n_fft, axis=0)  #zero pad to get non-cyclic autocorrelation
        PowerSpectralDensity = F * F.conjugate()
        inverse = ifft(PowerSpectralDensity,axis=0)
        autocorr = (inverse[:N]).real   #autocorr convention B
        n = np.arange(1, N+1)[::-1] 
        return autocorr/n[:, np.newaxis] #autocorr convention A


