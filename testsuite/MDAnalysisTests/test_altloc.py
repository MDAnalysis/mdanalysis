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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from MDAnalysis import Universe
import tempfile
import os
from numpy.testing import *
from MDAnalysisTests.datafiles import PDB_full


class TestAltloc(TestCase):
    def setUp(self):
        self.filename = PDB_full
        fd, self.outfile = tempfile.mkstemp(suffix=".pdb")  # output is always same as input
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

    def test_atomgroups(self):
        u = Universe(self.filename)
        segidB0 = len(u.select_atoms("segid B and (not altloc B)"))
        segidB1 = len(u.select_atoms("segid B and (not altloc A)"))
        assert_equal(segidB0, segidB1)
        altlocB0 = len(u.select_atoms("segid B and (altloc A)"))
        altlocB1 = len(u.select_atoms("segid B and (altloc B)"))
        assert_equal(altlocB0, altlocB1)
        sum = len(u.select_atoms("segid B"))
        assert_equal(sum, segidB0 + altlocB0)

    def test_bonds(self):
        u = Universe(self.filename, guess_bonds=True)
        # need to force topology to load before querying individual atom bonds
        u.build_topology()
        bonds0 = u.select_atoms("segid B and (altloc A)")[0].bonds
        bonds1 = u.select_atoms("segid B and (altloc B)")[0].bonds
        assert_equal(len(bonds0), len(bonds1))

    def test_write_read(self):
        u = Universe(self.filename)
        u.select_atoms("all").write(self.outfile)
        u2 = Universe(self.outfile)
        assert_equal(len(u.atoms), len(u2.atoms))
