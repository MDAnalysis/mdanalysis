# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
from __future__ import division, absolute_import, print_function

import pytest
from six.moves import range

import MDAnalysis as mda
from MDAnalysis.analysis.msd import MeanSquaredDisplacement as MSD

from numpy.testing import (assert_array_less,
                           assert_almost_equal, assert_equal)
import numpy as np
from scipy import fft,ifft

from MDAnalysisTests.datafiles import PSF, DCD, DCD

SELECTION = 'backbone and name CA and resid 1-10'

@pytest.fixture(scope='module')
def u():
    return mda.Universe(PSF, DCD)

@
def test_class_construct(u):
    m = MSD(u, SELECTION, msd_type='xyz', position_treatment='atom', mass_weighted=False)
    m.run()

def unit_step_traj(u):
    nstep = 1000
    x = np.arange(nstep)
    traj = np.hstack([x,y,z])
    m = MSD(u, SELECTION, msd_type='xyz', position_treatment='atom', mass_weighted=False)

