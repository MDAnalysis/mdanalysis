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


    def __init__(self, u, selection, t0, delta_t, msd_type='xyz', position_treatment='atom', mass_weighted=False)
                 kwargs=None, **basekwargs):
        self.u = u
        super(MeanSquaredDisplacement, self).__init__(self.u.trajectory, **basekwargs)

        self.selection = selection
        self.t0 = t0
        self.delta_t = delta_t
        self.msd_type = msd_type
        self.position_treatment = position_treatment
        self.mass_weighted = mass_weighted

        self.dim_fac = 0
        self.atoms = None
        self.timeseries = None

        
    def _prepare(self):
        parse_msd_type()
        check_masses()
        select_reference_positions()
        
    def check_masses(self):
        self.atoms = u.select_atoms(self.selection)
        if (mass_weighted or position_treatment == 'com'):
           
            #TODO check that all the atoms have masses
        else:
            pass
        
    def select_reference_positions(self):
        if self.position_treatment == 'atom'
            self.timeseries = atoms #TODO work out timeseries
        elif self.position_treatment == 'com':
            
            self.timeseries = calculate_com(atoms) #TODO work out timeseries
        else:
            raise ValueError('invalid position_treatment specified')

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
        

    
    def _run():

    
    def one_time_point(self,):
        for i in self.timeseries:
            displ_sq_sum = 0.0
            for j in self.n_atoms:
                displ = np.abs(r[self._dim]-r0[self._dim])
                displ_sq = displ*displ
                displ_sq_sum += displ_sq
            1.0/self.n_atoms



        


    def _conclude(self):

