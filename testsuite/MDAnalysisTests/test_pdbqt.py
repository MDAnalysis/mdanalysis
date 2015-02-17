# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
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

import MDAnalysis
from MDAnalysis.tests.datafiles import PDBQT_input, PDBQT_querypdb
import MDAnalysis.KDTree.NeighborSearch as kdNS

from numpy.testing import *
from nose.plugins.attrib import attr

import os
import tempfile


class TestPDBQT(TestCase):
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PDBQT_input)  # PDBQT
        self.query_universe = MDAnalysis.Universe(PDBQT_querypdb)  # PDB file
        fd, self.outfile = tempfile.mkstemp(suffix='.pdbqt')
        os.close(fd)

    def tearDown(self):
        del self.universe
        del self.query_universe
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

    def test_segid(self):
        sel = self.universe.selectAtoms('segid A')
        assert_equal(sel.numberOfAtoms(), 909, "failed to select segment A")
        sel = self.universe.selectAtoms('segid B')
        assert_equal(sel.numberOfAtoms(), 896, "failed to select segment B")

    def test_protein(self):
        sel = self.universe.selectAtoms('protein')
        assert_equal(sel.numberOfAtoms(), 1805, "failed to select protein")
        assert_equal(sel._atoms, self.universe.atoms._atoms,
                     "selected protein is not the same as auto-generated protein segment A+B")

    def test_backbone(self):
        sel = self.universe.selectAtoms('backbone')
        assert_equal(sel.numberOfAtoms(), 796)

    def test_neighborhood(self):
        '''test KDTree-based distance search around query atoms

        Creates a KDTree of the protein and uses the coordinates of
        the atoms in the query pdb to create a list of protein
        residues within 4.0A of the query atoms.
        '''
        protein = self.universe.selectAtoms("protein")
        ns_protein = kdNS.AtomNeighborSearch(protein)
        query_atoms = self.query_universe.atoms
        residue_neighbors = ns_protein.search_list(query_atoms, 4.0)
        assert_equal(len(residue_neighbors), 80)

    def test_writer(self):
        self.universe.A.write(self.outfile)
        uA = MDAnalysis.Universe(self.outfile)
        assert_equal(uA.atoms._atoms, self.universe.A.atoms._atoms,
                     "Writing and reading of chain A does not recover same atoms")
        del uA
