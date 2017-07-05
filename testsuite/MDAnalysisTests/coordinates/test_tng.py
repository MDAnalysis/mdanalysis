# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import division, absolute_import
from six.moves import zip

import MDAnalysis as mda
import numpy as np
import os
import warnings

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_array_almost_equal, dec,
                           assert_almost_equal, )


from MDAnalysisTests.datafiles import (GRO, TNG, PDB)

from MDAnalysisTests.datafiles import (COORDINATES_TOPOLOGY,
                                       COORDINATES_TNG)
from MDAnalysisTests.coordinates.base import (MultiframeReaderTest, BaseReference,
                                              # BaseWriterTest,
                                              assert_timestep_almost_equal)
from MDAnalysisTests import tempdir

# from MDAnalysis.coordinates.TNG import TNGReader

# I want to catch all warnings in the tests. If this is not set at the start it
# could cause test that check for warnings to fail.
warnings.simplefilter('always')




class TRRReference(BaseReference):
    def __init__(self):
        super(TRRReference, self).__init__()
        self.trajectory = COORDINATES_TNG
        self.topology = COORDINATES_TOPOLOGY
        self.changing_dimensions = True
        self.reader = mda.coordinates.TRR.TRRReader
        self.writer = mda.coordinates.TRR.TRRWriter
        self.ext = 'trr'
        self.prec = 3
        self.first_frame.velocities = self.first_frame.positions / 10
        self.first_frame.forces = self.first_frame.positions / 100

        self.second_frame.velocities = self.second_frame.positions / 10
        self.second_frame.forces = self.second_frame.positions / 100

        self.last_frame.velocities = self.last_frame.positions / 10
        self.last_frame.forces = self.last_frame.positions / 100

        self.jump_to_frame.velocities = self.jump_to_frame.positions / 10
        self.jump_to_frame.forces = self.jump_to_frame.positions / 100

    def iter_ts(self, i):
        ts = self.first_frame.copy()
        ts.positions = 2**i * self.first_frame.positions
        ts.velocities = ts.positions / 10
        ts.forces = ts.positions / 100
        ts.time = i
        ts.frame = i
        return ts


class TestTRRReader_2(MultiframeReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = TRRReference()
        super(TestTRRReader_2, self).__init__(reference)
