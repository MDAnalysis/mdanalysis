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

import numpy as np
import pytest

import MDAnalysis as mda
from MDAnalysis.coordinates.chemfiles import (
    ChemfilesReader, ChemfilesWriter, check_chemfiles_version,
)

from MDAnalysisTests import datafiles
from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest, BaseWriterTest, BaseReference
)
from MDAnalysisTests.coordinates.test_xyz import XYZReference


# skip entire test module if no appropriate chemfiles
chemfiles = pytest.importorskip('chemfiles')
@pytest.mark.skipif(not check_chemfiles_version(),
                    reason="Wrong version of chemfiles")
class TestChemFileXYZ(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture
    def ref():
        base = XYZReference()
        base.writer = ChemfilesWriter
        base.dimensions = np.array([0, 0, 0, 90, 90, 90], dtype=np.float32)

        return base

    @pytest.fixture
    def reader(self, ref):
        reader = ChemfilesReader(ref.trajectory)
        reader.add_auxiliary('lowf', ref.aux_lowf, dt=ref.aux_lowf_dt, initial_time=0, time_selector=None)
        reader.add_auxiliary('highf', ref.aux_highf, dt=ref.aux_highf_dt, initial_time=0, time_selector=None)
        return reader



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


@pytest.mark.skipif(not check_chemfiles_version(),
                    reason="Wrong version of chemfiles")
class TestChemfilesReader(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return ChemfilesXYZReference()


@pytest.mark.skipif(not check_chemfiles_version(),
                    reason="Wrong version of chemfiles")
class TestChemfilesWriter(BaseWriterTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return ChemfilesXYZReference()

    # Disable 'test_no_container' as it try to open a file for writing without
    # extension.
    def test_no_container(self, ref):
        pass

    def test_no_extension_raises(self, ref):
        with pytest.raises(chemfiles.ChemfilesError):
            ref.writer('foo')


@pytest.mark.skipif(not check_chemfiles_version(),
                    reason="Wrong version of chemfiles")
class TestChemfiles(object):
    def test_read_chemfiles_format(self):
        u = mda.Universe(
            datafiles.LAMMPSdata,
            format="chemfiles",
            topology_format="data",
            chemfiles_format="LAMMPS Data"
        )

        for ts in u.trajectory:
            assert ts.n_atoms == 18364

    def test_changing_system_size(self, tmpdir):
        outfile = "chemfiles-changing-size.xyz"
        with tmpdir.as_cwd():
            with open(outfile, "w") as fd:
                fd.write(VARYING_XYZ)

            u = mda.Universe(outfile, format="chemfiles", topology_format="XYZ")

            with pytest.raises(IOError):
                u.trajectory._read_next_timestep()

    def test_wrong_open_mode(self):
        with pytest.raises(IOError):
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
            assert (atom.name, atom.type, atom.record_type) in atoms

        for bond in check.bonds:
            assert (bond.atoms[0].ix, bond.atoms[1].ix) in bonds

    def test_write_topology(self, tmpdir):
        u = mda.Universe(datafiles.CONECT)
        outfile = "chemfiles-write-topology.pdb"
        with tmpdir.as_cwd():
            with ChemfilesWriter(outfile) as writer:
                writer.write(u)
            self.check_topology(datafiles.CONECT, outfile)

            # Manually setting the topology when creating the ChemfilesWriter
            # (1) from an object
            with ChemfilesWriter(outfile, topology=u) as writer:
                writer.write(u)
            self.check_topology(datafiles.CONECT, outfile)

            # (2) from a file
            with ChemfilesWriter(outfile, topology=datafiles.CONECT) as writer:
                writer.write(u)
            # FIXME: this does not work, since chemfiles also insert the bonds
            # which are implicit in PDB format (between standard residues), while
            # MDAnalysis only read the explicit CONNECT records.

            # self.check_topology(datafiles.CONECT, outfile)

    def test_write_velocities(self, tmpdir):
        u = mda.Universe.empty(4, trajectory=True)
        u.add_TopologyAttr('type', values=['H', 'H', 'H', 'H'])

        ts = u.trajectory.ts
        ts.dimensions = [20, 30, 41, 90, 90, 90]
        ts.positions = [
            [1, 1, 1],
            [2, 2, 2],
            [3, 3, 3],
            [4, 4, 4],
        ]
        ts.velocities = [
            [10, 10, 10],
            [20, 20, 20],
            [30, 30, 30],
            [40, 40, 40],
        ]

        outfile = "chemfiles-write-velocities.lmp"
        with tmpdir.as_cwd():
            with ChemfilesWriter(outfile, topology=u, chemfiles_format='LAMMPS Data') as writer:
                writer.write(u)

            with open(outfile) as file:
                content = file.read()
                assert content == EXPECTED_LAMMPS_DATA


VARYING_XYZ = """2

A 0 0 0
A 0 0 0
4

A 0 0 0
A 0 0 0
A 0 0 0
A 0 0 0
"""


EXPECTED_LAMMPS_DATA = """LAMMPS data file -- atom_style full -- generated by chemfiles
4 atoms
0 bonds
0 angles
0 dihedrals
0 impropers
1 atom types
0 bond types
0 angle types
0 dihedral types
0 improper types
0 20.0 xlo xhi
0 30.0 ylo yhi
0 41.0 zlo zhi

# Pair Coeffs
# 1 H

Masses

1 0.0 # H

Atoms # full

1 1 1 0.0 1.0 1.0 1.0 # H
2 2 1 0.0 2.0 2.0 2.0 # H
3 3 1 0.0 3.0 3.0 3.0 # H
4 4 1 0.0 4.0 4.0 4.0 # H

Velocities

1 10.0 10.0 10.0
2 20.0 20.0 20.0
3 30.0 30.0 30.0
4 40.0 40.0 40.0
"""
