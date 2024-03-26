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
import numpy as np
from numpy.testing import assert_allclose
import scipy
import pytest
from packaging.version import Version

import MDAnalysis as mda
from MDAnalysis.coordinates.chemfiles import ChemfilesReader, ChemfilesWriter
from MDAnalysis.coordinates.chemfiles import check_chemfiles_version

from MDAnalysisTests import datafiles
from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest,
    BaseWriterTest,
    BaseReference,
)
from MDAnalysisTests.coordinates.test_xyz import XYZReference


# skip entire test module if no appropriate chemfiles
chemfiles = pytest.importorskip("chemfiles")


@pytest.mark.parametrize("version", ["0.9.3", "0.11.0", "1.1.0"])
def test_version_check(version, monkeypatch):
    monkeypatch.setattr("chemfiles.__version__", version)
    assert not check_chemfiles_version()

    with pytest.raises(RuntimeError, match="Please install Chemfiles > 0.10"):
        ChemfilesReader("")

    with pytest.raises(RuntimeError, match="Please install Chemfiles > 0.10"):
        ChemfilesWriter("")


@pytest.mark.skipif(not check_chemfiles_version(), reason="Wrong version of chemfiles")
class TestChemfileXYZ(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture
    def ref():
        base = XYZReference()
        base.writer = ChemfilesWriter
        base.dimensions = None

        return base

    @pytest.fixture
    def reader(self, ref):
        reader = ChemfilesReader(ref.trajectory)
        reader.add_auxiliary(
            ref.aux_lowf,
            "lowf",
            dt=ref.aux_lowf_dt,
            initial_time=0,
            time_selector=None,
        )
        reader.add_auxiliary(
            ref.aux_highf,
            "highf",
            dt=ref.aux_highf_dt,
            initial_time=0,
            time_selector=None,
        )
        return reader


class ChemfilesXYZReference(BaseReference):
    def __init__(self):
        super(ChemfilesXYZReference, self).__init__()
        self.trajectory = datafiles.COORDINATES_XYZ
        self.topology = datafiles.COORDINATES_XYZ
        self.reader = ChemfilesReader
        self.writer = ChemfilesWriter
        self.ext = "xyz"
        self.volume = 0
        self.dimensions = None


@pytest.mark.skipif(not check_chemfiles_version(), reason="Wrong version of chemfiles")
class TestChemfilesReader(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return ChemfilesXYZReference()

    def test_copy(self, ref):
        # Issue #3664 - test not done in test_copying due to dependencies
        original = ChemfilesReader(
            ref.trajectory,
            convert_units=False,
            dt=2,
            time_offset=10,
            foo="bar",
        )
        copy = original.copy()

        assert original.format not in ("MEMORY", "CHAIN")
        assert original.convert_units is False
        assert copy.convert_units is False
        assert original._ts_kwargs["time_offset"] == 10
        assert copy._ts_kwargs["time_offset"] == 10
        assert original._ts_kwargs["dt"] == 2
        assert copy._ts_kwargs["dt"] == 2

        assert original.ts.data["time_offset"] == 10
        assert copy.ts.data["time_offset"] == 10

        assert original.ts.data["dt"] == 2
        assert copy.ts.data["dt"] == 2

        assert copy._kwargs["foo"] == "bar"

        # check coordinates
        assert original.ts.frame == copy.ts.frame
        assert_allclose(original.ts.positions, copy.ts.positions)

        original.next()
        copy.next()

        assert original.ts.frame == copy.ts.frame
        assert_allclose(original.ts.positions, copy.ts.positions)


@pytest.mark.skipif(not check_chemfiles_version(), reason="Wrong version of chemfiles")
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
            ref.writer("foo")


@pytest.mark.skipif(not check_chemfiles_version(), reason="Wrong version of chemfiles")
class TestChemfiles(object):
    def test_read_chemfiles_format(self):
        u = mda.Universe(
            datafiles.LAMMPSdata,
            format="chemfiles",
            topology_format="data",
            chemfiles_format="LAMMPS Data",
        )

        for ts in u.trajectory:
            assert ts.n_atoms == 18364

        # check that chemfiles format is passed to copies Issue #3664
        new_reader = u.trajectory.copy()
        assert new_reader._format == "LAMMPS Data"
        assert new_reader._kwargs["chemfiles_format"] == "LAMMPS Data"

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
        atoms = set([(atom.name, atom.type, atom.record_type) for atom in u.atoms])
        bonds = set([(bond.atoms[0].ix, bond.atoms[1].ix) for bond in u.bonds])

        check = mda.Universe(file)
        np.testing.assert_equal(
            u.trajectory.ts.positions,
            check.trajectory.ts.positions,
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

    def test_write_atom_group(self, tmpdir):
        u = mda.Universe(datafiles.CONECT)
        group = u.select_atoms("resname ARG")
        with tmpdir.as_cwd():
            outfile = "chemfiles-write-atom-group.pdb"
            with ChemfilesWriter(outfile) as writer:
                writer.write(group)

            check = mda.Universe(outfile)
            assert check.trajectory.ts.n_atoms == group.n_atoms

    def test_write_velocities(self, tmpdir):
        u = mda.Universe.empty(4, trajectory=True)

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

        outfile = "chemfiles-write-velocities.nc"
        with tmpdir.as_cwd():
            with ChemfilesWriter(outfile, topology=u) as writer:
                writer.write(u)

            with scipy.io.netcdf_file(outfile) as file:
                assert np.all(file.variables["coordinates"][0] == ts.positions)
                assert np.all(file.variables["velocities"][0] == ts.velocities)


VARYING_XYZ = """2

A 0 0 0
A 0 0 0
4

A 0 0 0
A 0 0 0
A 0 0 0
A 0 0 0
"""
