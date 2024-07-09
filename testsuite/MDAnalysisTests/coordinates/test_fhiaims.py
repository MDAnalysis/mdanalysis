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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import pytest
from io import StringIO

import MDAnalysis as mda
import numpy as np
from MDAnalysisTests.coordinates.base import BaseWriterTest
from MDAnalysisTests.datafiles import FHIAIMS
from numpy.testing import (assert_equal,
                           assert_almost_equal)


@pytest.fixture(scope='module')
def universe():
    return mda.Universe(FHIAIMS, FHIAIMS)


@pytest.fixture(scope='module')
def universe_from_one_file():
    return mda.Universe(FHIAIMS)


@pytest.fixture(scope='module')
def ref():
    return RefFHIAIMS()


@pytest.fixture(scope='class')
def good_input_natural_units():
    buffer = 'atom 0.1 0.2 0.3 H\natom 0.2 0.4 0.6 H\natom 0.3 0.6 0.9 H'
    return StringIO(buffer)


@pytest.fixture(scope='class')
def good_input_with_velocity():
    buffer = 'atom 0.1 0.1 0.1 H\nvelocity 0.1 0.1 0.1'
    return StringIO(buffer)


class RefFHIAIMS(object):
    from MDAnalysis.coordinates.FHIAIMS import (FHIAIMSReader, FHIAIMSWriter)

    filename, trajectory, topology = [FHIAIMS] * 3
    reader, writer = FHIAIMSReader, FHIAIMSWriter
    pos_atom1 = np.asarray(
        [6.861735,   2.103823,  37.753513], dtype=np.float32)
    dimensions = np.asarray(
        [18.6, 18.6,  55.8, 90., 90., 90.], dtype=np.float32)
    n_atoms = 6
    n_frames = 1
    time = 0.0
    ext = '.in'
    prec = 6
    container_format = True
    changing_dimensions = False


class TestFHIAIMSReader(object):

    prec = 6

    @pytest.fixture(scope='class')
    def bad_input_missing_periodicity(self):
        buffer = 'lattice_vector 1.0 0.0 0.0\nlattice_vector 0.0 1.0 0.0\natom 0.1 0.1 0.1 H'
        return StringIO(buffer)

    @pytest.fixture(scope='class')
    def bad_input_relative_positions_no_box(self):
        buffer = 'atom_frac 0.1 0.1 0.1 H'
        return StringIO(buffer)

    @pytest.fixture(scope='class')
    def bad_input_missing_velocity(self):
        buffer = 'atom 0.1 0.1 0.1 H\natom 0.2 0.2 0.2 H\nvelocity 0.1 0.1 0.1'
        return StringIO(buffer)

    @pytest.fixture(scope='class')
    def bad_input_velocity_wrong_position(self):
        buffer = 'atom 0.1 0.1 0.1 H\natom 0.2 0.2 0.2 H\nvelocity 0.1 0.1 0.1\nvelocity 0.1 0.1 0.1'
        return StringIO(buffer)

    @pytest.fixture(scope='class')
    def bad_input_wrong_input_line(self):
        buffer = 'garbage'
        return StringIO(buffer)

    @pytest.fixture(scope='class')
    def good_input_mixed_units(self):
        buffer = 'lattice_vector 1.0 0.0 0.0\nlattice_vector 0.0 2.0 0.0\nlattice_vector 0.0 0.0 3.0\natom 0.1 0.2 0.3 H\natom_frac 0.2 0.2 0.2 H\natom_frac 0.3 0.3 0.3 H'
        return StringIO(buffer)

    def test_single_file(self, universe, universe_from_one_file):
        assert_almost_equal(universe.atoms.positions, universe_from_one_file.atoms.positions,
                            self.prec, "FHIAIMSReader failed to load universe from single file")

    def test_uses_FHIAIMSReader(self, universe):
        from MDAnalysis.coordinates.FHIAIMS import FHIAIMSReader

        assert isinstance(universe.trajectory,
                          FHIAIMSReader), "failed to choose FHIAIMSReader"

    def test_dimensions(self, ref, universe):
        assert_almost_equal(
            ref.dimensions, universe.dimensions,
            self.prec, "FHIAIMSReader failed to get unitcell dimensions")

    def test_n_atoms(self, ref, universe):
        assert_equal(
            ref.n_atoms, universe.trajectory.n_atoms,
            "FHIAIMSReader failed to get the right number of atoms")

    def test_fhiaims_positions(self, ref, universe):
        # first particle
        assert_almost_equal(ref.pos_atom1,
                            universe.atoms.positions[0],
                            self.prec,
                            "FHIAIMSReader failed to read coordinates properly")

    def test_n_frames(self, ref, universe):
        assert_equal(ref.n_frames, universe.trajectory.n_frames,
                     "wrong number of frames")

    def test_time(self, ref, universe):
        assert_equal(ref.time, universe.trajectory.time,
                     "wrong time of the frame")

    def test_bad_input_missing_periodicity(self, bad_input_missing_periodicity):
        with pytest.raises(ValueError, match="Found partial periodicity"):
            u = mda.Universe(bad_input_missing_periodicity, format="FHIAIMS")

    def test_bad_input_relative_positions_no_box(self, bad_input_relative_positions_no_box):
        with pytest.raises(ValueError, match="Found relative coordinates in FHI-AIMS file without lattice info"):
            u = mda.Universe(
                bad_input_relative_positions_no_box, format="FHIAIMS")

    def test_bad_input_missing_velocity(self, bad_input_missing_velocity):
        with pytest.raises(ValueError, match="Found incorrect number of velocity tags"):
            u = mda.Universe(bad_input_missing_velocity, format="FHIAIMS")

    def test_bad_input_velocity_wrong_position(self, bad_input_velocity_wrong_position):
        with pytest.raises(ValueError, match="Non-conforming line .velocity must follow"):
            u = mda.Universe(
                bad_input_velocity_wrong_position, format="FHIAIMS")

    def test_bad_input_wrong_input_line(self, bad_input_wrong_input_line):
        with pytest.raises(ValueError, match="Non-conforming line"):
            u = mda.Universe(bad_input_wrong_input_line, format="FHIAIMS")

    def test_good_input_with_velocity(self, good_input_with_velocity):
        u = mda.Universe(good_input_with_velocity, format="FHIAIMS")
        assert_almost_equal(u.atoms.velocities[0],
                            np.asarray([0.1, 0.1, 0.1]),
                            self.prec,
                            "FHIAIMSReader failed to read velocities properly")

    def test_mixed_units(self, good_input_natural_units, good_input_mixed_units):
        u_natural = mda.Universe(good_input_natural_units, format="FHIAIMS")
        u_mixed = mda.Universe(good_input_mixed_units, format="FHIAIMS")
        print(u_natural.atoms.positions)
        print(u_mixed.atoms.positions)
        assert_almost_equal(u_natural.atoms.positions,
                            u_mixed.atoms.positions,
                            self.prec,
                            "FHIAIMSReader failed to read positions in lattice units properly")


class TestFHIAIMSWriter(BaseWriterTest):
    prec = 6
    ext = ".in"

    @pytest.fixture
    def outfile(self, tmpdir):
        return str(tmpdir.mkdir("FHIAIMSWriter").join('primitive-fhiaims-writer' + self.ext))

    def test_writer(self, universe, outfile):
        """Test writing from a single frame FHIAIMS file to a FHIAIMS file.
        """
        universe.atoms.write(outfile)
        u = mda.Universe(FHIAIMS, outfile)
        assert_almost_equal(u.atoms.positions,
                            universe.atoms.positions, self.prec,
                            err_msg="Writing FHIAIMS file with FHIAIMSWriter "
                                    "does not reproduce original coordinates")

    def test_writer_with_velocity(self, good_input_with_velocity, outfile):
        """Test writing from a single frame FHIAIMS file to a FHIAIMS file.
        """
        universe_in = mda.Universe(good_input_with_velocity, format="FHIAIMS")
        universe_in.atoms.write(outfile)
        u = mda.Universe(outfile)
        assert_almost_equal(u.atoms.velocities,
                            universe_in.atoms.velocities, self.prec,
                            err_msg="Writing FHIAIMS file with FHIAIMSWriter "
                                    "does not reproduce original velocities")

    def test_writer_no_atom_names(self, u_no_names, outfile):
        u_no_names.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.array(['X'] * u_no_names.atoms.n_atoms)
        assert_equal(u.atoms.names, expected)

    def test_writer_with_n_atoms_none(self, good_input_natural_units, outfile):
        u = mda.Universe(good_input_natural_units, format="FHIAIMS")
        with mda.Writer(outfile, natoms=None) as w:
            w.write(u.atoms)
            with open(outfile, 'r') as fhiaimsfile:
                line = fhiaimsfile.readline().strip()
                assert line.startswith(
                    'atom'), "Line written incorrectly with FHIAIMSWriter"
                assert line.endswith(
                    'H'), "Line written incorrectly with FHIAIMSWriter"
                line = np.asarray(line.split()[1:-1], dtype=np.float32)
                assert_almost_equal(line, [0.1, 0.2, 0.3], self.prec,
                                    err_msg="Writing FHIAIMS file with FHIAIMSWriter "
                                    "does not reproduce original positions")
