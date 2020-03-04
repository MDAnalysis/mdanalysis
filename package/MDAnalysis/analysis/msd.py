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

from __future__ import division, absolute_import
from six.moves import zip

import os
import errno
import warnings
import bz2
import functools

import numpy as np
from scipy import fft, ifft


import logging

import MDAnalysis
import MDAnalysis.lib.distances
from MDAnalysis.lib.util import openany
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.core.groups import AtomGroup
from .base import AnalysisBase


class MeanSquaredDisplacement(object):

    def __init__(self, u, selection, msd_type='xyz', position_treatment='atom', mass_weighted=False, kwargs=None, **basekwargs):
        
        self.u = u
        self.selection = selection
        self.msd_type = msd_type
        self.position_treatment = position_treatment
        self.mass_weighted = mass_weighted

        #local
        self.fft = False
        self.dim_fac = 0
        self._dim = None
        self.atoms = None
        self.n_frames = len(self.u.trajectory)
        self._position_array = None

        #result
        self.timeseries = None

        self.check_input()
        self._prepare()

    def check_input(self):
        self.check_masses()
        self.parse_msd_type()
        
    def _prepare(self):
        self.select_reference_positions()
        #self.construct_arrays()
        
    def check_masses(self):
        self.atoms = self.u.select_atoms(self.selection)
        if (self.mass_weighted or self.position_treatment == 'com'):
            masses = self.atoms.masses
            if masses.any == None:
                raise ValueError ('cannot have no mass for mass_weighted=True or position_treatment=com')
            else:
                pass
        else:
            pass
  
        
    def parse_msd_type(self):

        if self.msd_type == 'xyz': # full tensor
            self._dim = [0,1,2]
            self.dim_fac = 3.0
            
        elif self.msd_type == 'xy': # xy
            self._dim = [0,1]
            self.dim_fac = 2.0

        elif self.md_type == 'xz': # xz
            self._dime = [0,2]
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
        if self.position_treatment == 'atom':
            self._position_array = self.u.trajectory.timeseries(self.u.select_atoms(self.selection),order='fac') 
            self.N_particles = self._position_array.shape[1]
        elif self.position_treatment == 'com':
            raise NotImplementedError
            #TODO work out timeseries for com
        else:
            raise ValueError('invalid position_treatment specified')
    
    def run(self):
        
        if self.fft == True:
            self.timeseries = self._run_fft()
        else:
            self.timeseries = self._run_naieve()

    def _run_naieve(self): # naieve algorithm pre vectorisation / without FFT
        # _position_array is shape time, nparticles, 3
        msds_byparticle = np.zeros([self.n_frames, self.N_particles])
        lagtimes = np.arange(1,self.n_frames)
        msds_byparticle[0,:] = np.zeros(self.N_particles) # preset the zero lagtime so we dont have to iterate through
        for n in range(self.N_particles):
            for lag in lagtimes:
                disp = self._position_array[:-lag,n,self._dim if lag else None] - self._position_array[lag:,n,self._dim]
                sqdist = np.square(disp).sum(axis=1)
                msds_byparticle[lag,n] = sqdist.mean()
        msds = msds_byparticle.mean(axis=1)
        return msds

    def _run_fft(self):  #with FFT
        # _position_array is shape time, nparticles, 3
        N=r.shape(0)
        D=np.square(r).sum(axis=1) 
        D=np.append(D,0) 
        S2=sum([self.autocorrFFT(r[:, i]) for i in range(r.shape[1])])
        Q=2*D.sum()
        S1=np.zeros(N)
        for m in range(N):
            Q=Q-D[m-1]-D[N-m]
            S1[m]=Q/(N-m)
        return S1-2*S2

    @classmethod   
    def autocorrFFT(x):
        N=len(x)
        F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
        PowerSpectralDensity = F * F.conjugate()
        res = np.fft.ifft(PowerSpectralDensity)
        res= (res[:N]).real   #now we have the autocorrelation in convention B
        n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
        return res/n #this is the autocorrelation in convention A


