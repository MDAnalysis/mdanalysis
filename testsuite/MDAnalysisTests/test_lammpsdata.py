# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from numpy.testing import (
    assert_,
    assert_equal,
)

import MDAnalysis
from MDAnalysis.tests.datafiles import (
    LAMMPSdata,
)

from MDAnalysisTests.test_topology import _TestTopology



class TestLammpsData(_TestTopology):
    """Tests the reading of lammps .data topology files.

    The reading of coords and velocities is done separately in test_coordinates
    """
    topology = LAMMPSdata
    parser = MDAnalysis.topology.LAMMPSParser.DATAParser
    ref_n_atoms = 18360
    ref_numresidues = 24

    def test_charge(self):
        # No charges were supplied, should default to 0.0
        assert_equal(self.universe.atoms[0].charge, 0.0)

    def test_resid(self):
        assert_equal(len(self.universe.residues[0]), 765)

    # Testing _psf prevent building TGs
    # test length and random item from within
    def test_bonds(self):
        assert_equal(len(self.universe._topology['bonds']), 18336)
        assert_((5684, 5685) in self.universe._topology['bonds'])

    def test_angles(self):
        assert_equal(len(self.universe._topology['angles']), 29904)
        assert_((7575, 7578, 7579) in self.universe._topology['angles'])

    def test_dihedrals(self):
        assert_equal(len(self.universe._topology['dihedrals']), 5712)
        assert_((3210, 3212, 3215, 3218) in self.universe._topology['dihedrals'])

    def test_masses(self):
        assert_equal(self.universe.atoms[0].mass, 0.012)


