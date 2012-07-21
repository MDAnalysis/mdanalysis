# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

import MDAnalysis
from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small, PDB_closed, GRO, PDB, XTC, TRR
import MDAnalysis.core.AtomGroup
from MDAnalysis.core.AtomGroup import Atom, AtomGroup
from MDAnalysis import NoDataError

import numpy
from numpy.testing import *
from numpy import array, float32, rad2deg
from nose.plugins.attrib import attr

import os
import tempfile

try:
    from numpy.testing import assert_
except ImportError:
    # missing in numpy 1.2 but needed here:
    # copied code from numpy.testing 1.5
    def assert_(val, msg='') :
        """
        Assert that works in release mode.

        The Python built-in ``assert`` does not work when executing code in
        optimized mode (the ``-O`` flag) - no byte-code is generated for it.

        For documentation on usage, refer to the Python documentation.

        """
        if not val :
            raise AssertionError(msg)

class TestTuple1(TestCase):
    """Tests of Tuple(filename,format)."""
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(GRO, [(PDB,'pdb'),(XTC,'xtc'),(TRR,'trr')])
        self.numframes = self.universe.trajectory.numframes

    def tearDown(self):
        del self.universe
        del self.numframes

    def test_frame(self):
        assert_equal(self.numframes, 21)

class TestTuple2(TestCase):
    """Tests of Tuple(filename,format)."""
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, [(PDB_small,'pdb'),DCD])
        self.numframes = self.universe.trajectory.numframes

    def tearDown(self):
        del self.universe
        del self.numframes

    def test_frame(self):
        assert_equal(self.numframes, 99)

class TestChain(TestCase):
    """Tests of ChainReader with format overwritten."""
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, [PDB_small, PDB_closed],format='pdb')
        self.numframes = self.universe.trajectory.numframes

    def tearDown(self):
        del self.universe
        del self.numframes

    def test_frame(self):
        assert_equal(self.numframes, 2)

