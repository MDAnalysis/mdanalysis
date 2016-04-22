# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import MDAnalysis
from MDAnalysis.tests.datafiles import PDBQT_input, PDBQT_querypdb
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch

from numpy.testing import assert_equal, TestCase

import os
from MDAnalysisTests import tempdir


class TestPDBQT(TestCase):
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PDBQT_input)  # PDBQT
        self.query_universe = MDAnalysis.Universe(PDBQT_querypdb)  # PDB file
        self.tempdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tempdir.name, 'test.pdbqt')

    def tearDown(self):
        del self.universe
        del self.query_universe
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

    def test_segid(self):
        sel = self.universe.select_atoms('segid A')
        assert_equal(sel.n_atoms, 909, "failed to select segment A")
        sel = self.universe.select_atoms('segid B')
        assert_equal(sel.n_atoms, 896, "failed to select segment B")

    def test_protein(self):
        sel = self.universe.select_atoms('protein')
        assert_equal(sel.n_atoms, 1805, "failed to select protein")
        assert_equal(sel._atoms, self.universe.atoms._atoms,
                     "selected protein is not the same as auto-generated protein segment A+B")

    def test_backbone(self):
        sel = self.universe.select_atoms('backbone')
        assert_equal(sel.n_atoms, 796)

    def test_neighborhood(self):
        '''test KDTree-based distance search around query atoms

        Creates a KDTree of the protein and uses the coordinates of
        the atoms in the query pdb to create a list of protein
        residues within 4.0A of the query atoms.
        '''
        protein = self.universe.select_atoms("protein")
        ns_protein = AtomNeighborSearch(protein)
        query_atoms = self.query_universe.atoms
        residue_neighbors = ns_protein.search(query_atoms, 4.0)
        assert_equal(len(residue_neighbors), 80)

    def test_writer(self):
        self.universe.A.write(self.outfile)
        uA = MDAnalysis.Universe(self.outfile)
        assert_equal(uA.atoms._atoms, self.universe.A.atoms._atoms,
                     "Writing and reading of chain A does not recover same atoms")
        del uA

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1, "wrong number of frames in pdb")

    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0, "wrong time of the frame")

    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0,
                     "wrong frame number (0-based, should be 0 for single frame readers)")
