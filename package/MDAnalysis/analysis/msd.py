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

"""
Mean Squared Displacement --- :mod:`MDAnalysis.analysis.msd`
==============================================================

This module implements the calculation of Mean Squared Displacmements (MSD).
MSDs can be used to characterise the speed at which particles move and has its roots
in the study of Brownian motion.  For a thorough review see XXX et al. MSDs are computed from
the following expression

Where XX represents an ensemble average over 

The computation of the MSD in this way can be computationally intensive due to it's N^2 scaling. 
An algorithm to compute the MSD with Nlogn(N) scaling based on a Fast Fourier Transform is known and can be accessed by setting FFT=True.




Classes and Functions
---------------------

.. autoclass:: MeanSquaredDisplacement


"""

from __future__ import division, absolute_import
from six.moves import zip

import os
import errno
import warnings
import bz2
import functools

import numpy as np
from scipy import fft,ifft
import logging
import MDAnalysis



class MeanSquaredDisplacement(object):
    r"""Class representing a density on a regular cartesian grid.

    Parameters
    ----------∏
    u : 
        An MDAnalysis Universe :class:`Universe`
    selection : 
        An MDAnalysis selection string
    
    

    Attributes
    ----------




    Notes
    -----
    Notes


    See Also
    --------
   

    Examples
    --------
    Typical use:

    

    """

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
        # _position_array is shape time, nparticles, 3
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

    def _run_fft(self):  #with FFT
        # _position_array is shape time, nparticles, 3
        msds_byparticle = []
        reshape_positions = self._position_array[:,:,self._dim]
        N=reshape_positions.shape[0]
        D=np.square(reshape_positions).sum(axis=2, dtype=np.float64) 
        D=np.append(D,np.zeros(reshape_positions.shape[:2]), axis=0) 
        Q=2*D.sum(axis=0)
        S1=np.zeros(reshape_positions.shape[:2])
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
        N=(x.shape[0])
        F = fft(x, n=2*N, axis=0)  #zero pad to get non-cyclic autocorrelation
        PowerSpectralDensity = F * F.conjugate()
        inverse = ifft(PowerSpectralDensity,axis=0)
        autocorr = (inverse[:N]).real   #autocorr convention B
        n = np.arange(1, N+1)[::-1] 
        return autocorr/n[:, np.newaxis] #autocorr convention A


