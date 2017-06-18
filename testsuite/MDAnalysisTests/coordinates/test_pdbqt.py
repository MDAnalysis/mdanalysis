# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

import MDAnalysis as mda
from MDAnalysisTests.datafiles import PDBQT_input, PDBQT_querypdb
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch

from numpy.testing import (
    assert_,
    assert_equal,
    assert_array_equal,
    assert_warns,
)

import os
from MDAnalysisTests import tempdir, make_Universe


class TestPDBQT(object):
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = mda.Universe(PDBQT_input)  # PDBQT

    def tearDown(self):
        del self.universe

    def test_segid(self):
        sel = self.universe.select_atoms('segid A')
        assert_equal(sel.n_atoms, 909, "failed to select segment A")
        sel = self.universe.select_atoms('segid B')
        assert_equal(sel.n_atoms, 896, "failed to select segment B")

    def test_protein(self):
        sel = self.universe.select_atoms('protein')
        assert_equal(sel.n_atoms, 1805, "failed to select protein")
        assert_array_equal(sel.atoms.ix, self.universe.atoms.ix,
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
        query_universe = mda.Universe(PDBQT_querypdb)  # PDB file

        protein = self.universe.select_atoms("protein")
        ns_protein = AtomNeighborSearch(protein)
        query_atoms = query_universe.atoms
        residue_neighbors = ns_protein.search(query_atoms, 4.0)
        assert_equal(len(residue_neighbors), 80)

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1, "wrong number of frames in pdb")

    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0, "wrong time of the frame")

    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0,
                     "wrong frame number (0-based, should be 0 for single frame readers)")


class TestPDBQTWriter(object):
    def setUp(self):
        self.reqd_attributes = ['names', 'types', 'resids', 'resnames', 'radii', 'charges']
        self.tmpdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tmpdir.name, 'out.pdbqt')

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.tmpdir
        del self.outfile
        del self.reqd_attributes

    def test_roundtrip_writing_coords(self):
        u = mda.Universe(PDBQT_input)
        u.atoms.write(self.outfile)
        u2 = mda.Universe(self.outfile)

        assert_array_equal(u2.atoms.positions, u.atoms.positions,
                           "Round trip does not preserve coordinates")

    def test_roundtrip_formatting(self):
        # Compare formatting of first line
        u = mda.Universe(PDBQT_input)
        u.atoms.write(self.outfile)

        with open(PDBQT_input, 'r') as inf:
            l_ref = inf.readline().strip()
        with open(self.outfile, 'r') as inf:
            inf.readline()  # header
            inf.readline()  # cryst
            l_new = inf.readline().strip()
        assert_(l_ref == l_new)

    @staticmethod
    def assert_writing_warns(u, outfile):
        # write the test universe, and check warning is raised
        assert_warns(UserWarning, u.atoms.write, outfile)

    def test_write_no_charges(self):
        attrs = self.reqd_attributes
        attrs.remove('charges')
        u = make_Universe(attrs, trajectory=True)

        self.assert_writing_warns(u, self.outfile)

        u2 = mda.Universe(self.outfile)

        assert_(all(u2.atoms.charges == 0.0))

    def test_write_no_chainids_with_segids(self):
        attrs = self.reqd_attributes
        attrs.append('segids')
        u = make_Universe(attrs, trajectory=True)

        u.atoms.write(self.outfile)
        u2 = mda.Universe(self.outfile)

        # Should have used last letter of segid as chainid
        assert_(all(u2.atoms[:25].segids == 'A'))
        assert_(all(u2.atoms[25:50].segids == 'B'))

    def test_get_writer(self):
        u = mda.Universe(PDBQT_input)
        w = u.trajectory.Writer(self.outfile)

        assert_(isinstance(w, mda.coordinates.PDBQT.PDBQTWriter))
        w.close()
