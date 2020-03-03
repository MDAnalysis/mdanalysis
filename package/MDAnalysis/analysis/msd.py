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

import logging

import MDAnalysis
import MDAnalysis.lib.distances
from MDAnalysis.lib.util import openany
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.core.groups import AtomGroup
from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.contacts")



class MeanSquaredDisplacement(AnalysisBase):


    def __init__(self, u, selection, msd_type='xyz', position_treatment='atom', mass_weighted=False)
                 kwargs=None, **basekwargs):
        self.u = u
        super(MeanSquaredDisplacement, self).__init__(self.u.trajectory, **basekwargs)

        self.selection = selection
        self.msd_type = msd_type
        self.position_treatment = position_treatment
        self.mass_weighted = mass_weighted

        #local
        self.dim_fac = 0
        self._dim = None
        self.atoms = None
        self.n_frames = self.u.n_frames
        self._position_array = None
        self.

        #result
        self.timeseries = None

        self.check_input()

    def check_input(self):
        check_masses()
        parse_msd_type()
        
    def _prepare(self):
        select_reference_positions()
        construct_arrays()
        
    def check_masses(self):
        self.atoms = u.select_atoms(self.selection)
        if (self.mass_weighted or self.position_treatment == 'com'):
            masses = self.atoms.masses
            if masses.any == None:
               raise ValueError ('cannot have no mass for mass_weighted=True or position_treatment=com')
            raise NotImplementedError
        else:
            pass
  
        
    def parse_msd_type(self):

        if self.msd_type = 'xyz':
            self._dim = [0,1,2]
            self.dim_fac = 3.0
            
        elif self.msd_type = 'xy':
            self._dim = [0,1]
            self.dim_fac = 2.0

        elif self.md_type = 'xz':
            self._dime = [0,2]
            self.dim_fac = 2.0

        elif self.msd_type = 'yz':
            self._dim = [1,2]
            self.dim_fac = 2.0

        elif self.msd_type = 'x':
            self._dim = [0]
            self.dim_fac = 1.0

        elif self.msd_type = 'y':
            self._dim = [1]
            self.dim_fac = 1.0

        elif self.msd_type = 'z':
            self._dim = [2]
            self.dim_fac = 1.0

        else:
            raise ValueError('invalid msd_type specified')

    def select_reference_positions(self):
        if self.position_treatment == 'atom'
            self._position_array = self.u.trajectory.timeseries(self.u.select_atoms(self.selection),order='fac') 
            self.N_particles = self._position_array.shape[0]
        elif self.position_treatment == 'com':
            raise NotImplementedError
            #TODO work out timeseries for com
        else:
            raise ValueError('invalid position_treatment specified')
        
    def _run(): # naieve algorithm pre vectorisation / without FFT
        # r is shape time, nparticles, 3
        msds_byparticle = np.zeros([self.n_frames, self.N_particles])
        lagtimes = np.arange(self.n_frames)
        for n in range(self.N_particles):
            for i,lag in enumerate(lag_times): 
                disp = r[:-lag,n,self._dim if lag else None] - r[lag:,n,self._dim] )
                sqdist = np.square(disp).sum(axis=1)
                msds[i,n] = sqdist.mean()
        msds = msds_byparticle.mean(axis=XX)
            


            

    




        


    def _conclude(self):

