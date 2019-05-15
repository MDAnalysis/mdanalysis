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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

import pytest

import MDAnalysis as mda
import numpy as np
from unittest import TestCase

from MDAnalysis.coordinates.chemfiles import ChemfilesReader, ChemfilesWriter

from MDAnalysisTests import datafiles, tempdir
from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest, BaseWriterTest, BaseReference
)


class ChemfilesXYZReference(BaseReference):
    def __init__(self):
        super(ChemfilesXYZReference, self).__init__()
        self.trajectory = datafiles.COORDINATES_XYZ
        self.topology = datafiles.COORDINATES_XYZ
        self.reader = ChemfilesReader
        self.writer = ChemfilesWriter
        self.ext = 'xyz'
        self.volume = 0
        self.dimensions = np.zeros(6)
        self.dimensions[3:] = 90.0


class TestChemfilesReader(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return ChemfilesXYZReference()


class TestChemfilesWriter(BaseWriterTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return ChemfilesXYZReference()

    # Disable 'test_no_container' as it try to open a file for writing without
    # extension.
    def test_no_container(self, ref):
        pass


class TestChemfiles(TestCase):
    def test_read_chemfiles_format(self):
        u = mda.Universe(
            datafiles.LAMMPSdata,
            format="chemfiles",
            topology_format="data",
            chemfiles_format="LAMMPS Data"
        )

        for ts in u.trajectory:
            self.assertEqual(ts.n_atoms, 18364)

    def test_changing_system_size(self):
        outfile = tempdir.TempDir().name + "chemfiles-changing-size.xyz"
        with open(outfile, "w") as fd:
            fd.write(VARYING_XYZ)

        u = mda.Universe(outfile, format="chemfiles", topology_format="XYZ")

        with self.assertRaises(IOError):
            u.trajectory._read_next_timestep()

    def test_wrong_open_mode(self):
        with self.assertRaises(IOError):
            _ = ChemfilesWriter("", mode="r")

    def check_topology(self, reference, file):
        u = mda.Universe(reference)
        atoms = set([
            (atom.name, atom.type, atom.record_type)
            for atom in u.atoms
        ])
        bonds = set([
            (bond.atoms[0].ix, bond.atoms[1].ix)
            for bond in u.bonds
        ])

        check = mda.Universe(file)
        np.testing.assert_equal(
            u.trajectory.ts.positions,
            check.trajectory.ts.positions
        )

        for atom in check.atoms:
            self.assertIn((atom.name, atom.type, atom.record_type), atoms)

        for bond in check.bonds:
            self.assertIn((bond.atoms[0].ix, bond.atoms[1].ix), bonds)

    def test_write_topology(self):
        u = mda.Universe(datafiles.CONECT)
        outfile = tempdir.TempDir().name + "chemfiles-write-topology.pdb"
        with ChemfilesWriter(outfile) as writer:
            writer.write(u)
        self.check_topology(datafiles.CONECT, outfile)

        # Manually setting the topology when creating the ChemfilesWriter
        # (1) from an object
        with ChemfilesWriter(outfile, topology=u) as writer:
            writer.write_next_timestep(u.trajectory.ts)
        self.check_topology(datafiles.CONECT, outfile)

        # (2) from a file
        with ChemfilesWriter(outfile, topology=datafiles.CONECT) as writer:
            writer.write_next_timestep(u.trajectory.ts)
        # FIXME: this does not work, since chemfiles also insert the bonds
        # which are implicit in PDB format (bewteen standard residues), while
        # MDAnalysis only read the explicit CONNECT records.

        # self.check_topology(datafiles.CONECT, outfile)


VARYING_XYZ = """2

A 0 0 0
A 0 0 0
4

A 0 0 0
A 0 0 0
A 0 0 0
A 0 0 0
"""
