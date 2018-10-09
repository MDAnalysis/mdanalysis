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
from six import StringIO
from six.moves import zip
import os
from unittest import TestCase

import MDAnalysis as mda
import numpy as np
from MDAnalysisTests import tempdir, make_Universe
from MDAnalysisTests.coordinates.base import _SingleFrameReader
from MDAnalysisTests.coordinates.reference import (RefAdKSmall, Ref4e43,
                                                   RefAdK)
from MDAnalysisTests.datafiles import (PDB, PDB_small, PDB_multiframe,
                                       XPDB_small, PSF, DCD, CONECT, CRD,
                                       INC_PDB, PDB_xlserial, ALIGN, ENT,
                                       PDB_cm, PDB_cm_gz, PDB_cm_bz2,
                                       PDB_mc, PDB_mc_gz, PDB_mc_bz2)
from numpy.testing import (assert_equal,
                           assert_array_almost_equal,
                           assert_almost_equal)


class TestPDBReader(_SingleFrameReader):
    __test__ = True
    def setUp(self):
        # can lead to race conditions when testing in parallel
        self.universe = mda.Universe(RefAdKSmall.filename)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3

    def test_uses_PDBReader(self):
        from MDAnalysis.coordinates.PDB import PDBReader

        assert isinstance(self.universe.trajectory, PDBReader), "failed to choose PDBReader"

    def test_dimensions(self):
        assert_almost_equal(
            self.universe.trajectory.ts.dimensions, RefAdKSmall.ref_unitcell,
            self.prec,
            "PDBReader failed to get unitcell dimensions from CRYST1")

    def test_ENT(self):
        from MDAnalysis.coordinates.PDB import PDBReader
        self.universe = mda.Universe(ENT)
        assert isinstance(self.universe.trajectory, PDBReader), "failed to choose PDBReader"


class _PDBMetadata(TestCase, Ref4e43):
    __test__ = False

    def setUp(self):
        self.universe = mda.Universe(self.filename)

    def tearDown(self):
        del self.universe

    def test_HEADER(self):
        assert_equal(self.universe.trajectory.header,
                     self.header,
                     err_msg="HEADER record not correctly parsed")

    def test_TITLE(self):
        try:
            title = self.universe.trajectory.title
        except AttributeError:
            raise AssertionError("Reader does not have a 'title' attribute.")
        assert_equal(len(title),
                     len(self.title),
                     err_msg="TITLE does not contain same number of lines")
        for lineno, (parsed, reference) in enumerate(zip(title, self.title),
                                                     start=1):
            assert_equal(parsed,
                         reference,
                         err_msg="TITLE line {0} do not match".format(lineno))

    def test_COMPND(self):
        try:
            compound = self.universe.trajectory.compound
        except AttributeError:
            raise AssertionError(
                "Reader does not have a 'compound' attribute.")
        assert_equal(len(compound),
                     len(self.compnd),
                     err_msg="COMPND does not contain same number of lines")
        for lineno, (parsed, reference) in enumerate(zip(compound,
                                                         self.compnd),
                                                     start=1):
            assert_equal(parsed,
                         reference,
                         err_msg="COMPND line {0} do not match".format(lineno))

    def test_REMARK(self):
        try:
            remarks = self.universe.trajectory.remarks
        except AttributeError:
            raise AssertionError("Reader does not have a 'remarks' attribute.")
        assert_equal(len(remarks),
                     self.num_remarks,
                     err_msg="REMARK does not contain same number of lines")
        # only look at the first 5 entries
        for lineno, (parsed, reference) in enumerate(
                zip(remarks[:self.nmax_remarks],
                    self.remarks[:self.nmax_remarks]),
                start=1):
            assert_equal(parsed,
                         reference,
                         err_msg="REMARK line {0} do not match".format(lineno))


class TestExtendedPDBReader(_SingleFrameReader):
    __test__ = True
    def setUp(self):
        self.universe = mda.Universe(PDB_small,
                                     topology_format="XPDB",
                                     format="XPDB")
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3

    def test_long_resSeq(self):
        # it checks that it can read a 5-digit resid
        self.universe = mda.Universe(XPDB_small, topology_format="XPDB")
        u = self.universe.select_atoms(
            'resid 1 or resid 10 or resid 100 or resid 1000 or resid 10000')
        assert_equal(u[4].resid, 10000, "can't read a five digit resid")


class TestPDBWriter(object):
    # 3 decimals in PDB spec
    # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
    prec = 3
    ext = ".pdb"

    @pytest.fixture
    def universe(self):
        return mda.Universe(PSF, PDB_small)

    @pytest.fixture
    def universe2(self):
        return mda.Universe(PSF, DCD)

    @pytest.fixture
    def outfile(self, tmpdir):
        return str(tmpdir.mkdir("PDBWriter").join('primitive-pdb-writer' + self.ext))

    @pytest.fixture
    def u_no_resnames(self):
        return make_Universe(['names', 'resids'], trajectory=True)

    @pytest.fixture
    def u_no_resids(self):
        return make_Universe(['names', 'resnames'], trajectory=True)

    @pytest.fixture
    def u_no_names(self):
        return make_Universe(['resids', 'resnames'], trajectory=True)


    def test_writer(self, universe, outfile):
        "Test writing from a single frame PDB file to a PDB file." ""
        universe.atoms.write(outfile)
        u = mda.Universe(PSF, outfile)
        assert_almost_equal(u.atoms.positions,
                            universe.atoms.positions, self.prec,
                            err_msg="Writing PDB file with PDBWriter "
                                    "does not reproduce original coordinates")

    def test_writer_no_resnames(self, u_no_resnames, outfile):
        u_no_resnames.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.array(['UNK'] * u_no_resnames.atoms.n_atoms)
        assert_equal(u.atoms.resnames, expected)

    def test_writer_no_resids(self, u_no_resids, outfile):
        u_no_resids.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.ones((25,))
        assert_equal(u.residues.resids, expected)

    def test_writer_no_atom_names(self, u_no_names, outfile):
        u_no_names.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.array(['X'] * u_no_names.atoms.n_atoms)
        assert_equal(u.atoms.names, expected)

    def test_writer_no_altlocs(self, u_no_names, outfile):
        u_no_names.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.array([''] * u_no_names.atoms.n_atoms)
        assert_equal(u.atoms.altLocs, expected)

    def test_writer_no_icodes(self, u_no_names, outfile):
        u_no_names.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.array([''] * u_no_names.atoms.n_atoms)
        assert_equal(u.atoms.icodes, expected)

    def test_writer_no_segids(self, u_no_names, outfile):
        u_no_names.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.array(['SYSTEM'] * u_no_names.atoms.n_atoms)
        assert_equal([atom.segid for atom in u.atoms], expected)

    def test_writer_no_occupancies(self, u_no_names, outfile):
        u_no_names.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.ones(u_no_names.atoms.n_atoms)
        assert_equal(u.atoms.occupancies, expected)

    def test_writer_no_tempfactors(self, u_no_names, outfile):
        u_no_names.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.zeros(u_no_names.atoms.n_atoms)
        assert_equal(u.atoms.tempfactors, expected)

    def test_write_single_frame_Writer(self, universe2, outfile):
        """Test writing a single frame from a DCD trajectory to a PDB using
        MDAnalysis.Writer (Issue 105)"""
        u = universe2
        u.trajectory[50]
        with  mda.Writer(outfile) as W:
            W.write(u.select_atoms('all'))
        u2 = mda.Universe(outfile)
        assert_equal(u2.trajectory.n_frames,
                     1,
                     err_msg="The number of frames should be 1.")

    def test_write_single_frame_AtomGroup(self, universe2, outfile):
        """Test writing a single frame from a DCD trajectory to a PDB using
        AtomGroup.write() (Issue 105)"""
        u = universe2
        u.trajectory[50]
        u.atoms.write(outfile)
        u2 = mda.Universe(PSF, outfile)
        assert_equal(u2.trajectory.n_frames,
                     1,
                     err_msg="Output PDB should only contain a single frame")
        assert_almost_equal(u2.atoms.positions, u.atoms.positions,
                            self.prec, err_msg="Written coordinates do not "
                                               "agree with original coordinates from frame %d" %
                                               u.trajectory.frame)

    def test_check_coordinate_limits_min(self, universe, outfile):
        """Test that illegal PDB coordinates (x <= -999.9995 A) are caught
        with ValueError (Issue 57)"""
        # modify coordinates (universe needs to be a per-function fixture)
        u = universe
        u.atoms[2000].position = [0, -999.9995, 22.8]
        with pytest.raises(ValueError):
            u.atoms.write(outfile)

    def test_check_coordinate_limits_max(self, universe, outfile):
        """Test that illegal PDB coordinates (x > 9999.9995 A) are caught
        with ValueError (Issue 57)"""
        # modify coordinates (universe needs to be a per-function fixture)
        u = universe
        # OB: 9999.99951 is not caught by '<=' ?!?
        u.atoms[1000].position = [90.889, 9999.9996, 12.2]
        with pytest.raises(ValueError):
            u.atoms.write(outfile)

    def test_check_HEADER_TITLE_multiframe(self, universe2, outfile):
        """Check whether HEADER and TITLE are written just once in a multi-
        frame PDB file (Issue 741)"""
        u = universe2
        protein = u.select_atoms("protein and name CA")
        with mda.Writer(outfile, multiframe=True) as pdb:
            for ts in u.trajectory[:5]:
                pdb.write(protein)

        with open(outfile) as f:
            got_header = 0
            got_title = 0
            for line in f:
                if line.startswith('HEADER'):
                    got_header += 1
                    assert got_header <= 1, "There should be only one HEADER."
                elif line.startswith('TITLE'):
                    got_title += 1
                    assert got_title <= 1, "There should be only one TITLE."

    @pytest.mark.parametrize("startframe,maxframes",
                             [(0, 12), (9997, 12)])
    def test_check_MODEL_multiframe(self, universe2, outfile, startframe, maxframes):
        """Check whether MODEL number is in the right column (Issue #1950)"""
        u = universe2
        protein = u.select_atoms("protein and name CA")
        with mda.Writer(outfile, multiframe=True, start=startframe) as pdb:
            for ts in u.trajectory[:maxframes]:
                pdb.write(protein)

        def get_MODEL_lines(filename):
            with open(filename) as pdb:
                for line in pdb:
                    if line.startswith("MODEL"):
                        yield line

        MODEL_lines = list(get_MODEL_lines(outfile))

        assert len(MODEL_lines) == maxframes
        for model, line in enumerate(MODEL_lines, start=startframe+1):
            # test that only the right-most 4 digits are stored (rest must be space)
            # line[10:14] == '9999' or '   1'

            # test appearance with white space
            assert line[5:14] == "{0:>9d}".format(int(str(model)[-4:]))

            # test number (only last 4 digits)
            assert int(line[10:14]) == model % 10000


@pytest.mark.xfail(os.name == 'nt',
                   strict=True,
                   reason="PDB multiframe reading not yet supported "
                          "on Windows.")
class TestMultiPDBReader(TestCase):
    def setUp(self):
        self.multiverse = mda.Universe(PDB_multiframe,
                                       guess_bonds=True)
        self.conect = mda.Universe(CONECT, guess_bonds=True)

    def tearDown(self):
        del self.multiverse
        del self.conect

    def test_n_frames(self):
        assert_equal(self.multiverse.trajectory.n_frames, 24,
                     "Wrong number of frames read from PDB muliple model file")

    def test_n_atoms_frame(self):
        u = self.multiverse
        desired = 392
        for frame in u.trajectory:
            assert_equal(len(u.atoms), desired, err_msg="The number of atoms "
                                                        "in the Universe (%d) does not" " match the number "
                                                        "of atoms in the test case (%d) at frame %d" % (
                                                            len(u.atoms), desired, u.trajectory.frame))

    def test_rewind(self):
        u = self.multiverse
        u.trajectory[11]
        assert_equal(u.trajectory.ts.frame, 11,
                     "Failed to forward to 11th frame (frame index 11)")
        u.trajectory.rewind()
        assert_equal(u.trajectory.ts.frame, 0,
                     "Failed to rewind to 0th frame (frame index 0)")

    def test_iteration(self):
        u = self.multiverse
        frames = []
        for frame in u.trajectory:
            pass
        # should rewind after previous test
        # problem was: the iterator is NoneType and next() cannot be called
        for ts in u.trajectory:
            frames.append(ts)
        assert_equal(
            len(frames), u.trajectory.n_frames,
            "iterated number of frames %d is not the expected number %d; "
            "trajectory iterator fails to rewind" %
            (len(frames), u.trajectory.n_frames))

    def test_slice_iteration(self):
        u = self.multiverse
        frames = []
        for ts in u.trajectory[4:-2:4]:
            frames.append(ts.frame)
        assert_equal(np.array(frames),
                     np.arange(u.trajectory.n_frames)[4:-2:4],
                     err_msg="slicing did not produce the expected frames")

    def test_conect_bonds_conect(self):
        conect = self.conect
        assert_equal(len(conect.atoms), 1890)
        assert_equal(len(conect.bonds), 1922)

        with tempdir.in_tempdir():
            try:
                outfile = 'test-pdb-hbonds.pdb'
                self.conect.atoms.write(outfile, bonds="conect")
                u1 = mda.Universe(outfile, guess_bonds=True)
            finally:
                os.unlink(outfile)
            assert_equal(len(u1.atoms), 1890)
            assert_equal(len(u1.bonds), 1922)


    def test_numconnections(self):
        u = self.multiverse

        # the bond list is sorted - so swaps in input pdb sequence should not
        # be a problem
        desired = [[48, 365],
                   [99, 166],
                   [166, 99],
                   [249, 387],
                   [313, 331],
                   [331, 313, 332, 340],
                   [332, 331, 333, 338, 341],
                   [333, 332, 334, 342, 343],
                   [334, 333, 335, 344, 345],
                   [335, 334, 336, 337],
                   [336, 335],
                   [337, 335, 346, 347, 348], [338, 332, 339, 349],
                   [339, 338],
                   [340, 331],
                   [341, 332],
                   [342, 333],
                   [343, 333],
                   [344, 334],
                   [345, 334],
                   [346, 337],
                   [347, 337],
                   [348, 337],
                   [349, 338],
                   [365, 48],
                   [387, 249]]

        def helper(atoms, bonds):
            """
            Convert a bunch of atoms and bonds into a list of CONECT records
            """
            con = {}

            for bond in bonds:
                a1, a2 = bond[0].index, bond[1].index
                if a1 not in con:
                    con[a1] = []
                if a2 not in con:
                    con[a2] = []
                con[a2].append(a1)
                con[a1].append(a2)

            atoms = sorted([a.index for a in atoms])

            conect = [([a, ] + sorted(con[a])) for a in atoms if a in con]
            conect = [[a + 1 for a in c] for c in conect]

            return conect

        conect = helper(self.multiverse.atoms, [b for b in u.bonds
                                                if not b.is_guessed])
        assert_equal(conect, desired, err_msg="The bond list does not match "
                                              "the test reference; len(actual) is %d, len(desired) "
                                              "is %d" % (len(u._topology.bonds.values), len(desired)))

def test_conect_bonds_all(tmpdir):
    conect = mda.Universe(CONECT, guess_bonds=True)

    assert_equal(len(conect.atoms), 1890)
    assert_equal(len(conect.bonds), 1922)

    outfile = os.path.join(str(tmpdir), 'pdb-connect-bonds.pdb')
    conect.atoms.write(outfile, bonds="all")
    u2 = mda.Universe(outfile, guess_bonds=True)

    assert_equal(len(u2.atoms), 1890)
    assert_equal(len([b for b in u2.bonds if not b.is_guessed]), 1922)

    # assert_equal(len([b for b in conect.bonds if not b.is_guessed]), 1922)

@pytest.mark.xfail(os.name == 'nt',
                   strict=True,
                   reason="PDB multiframe reading not yet supported "
                          "on Windows.")
class TestMultiPDBWriter(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small)
        self.multiverse = mda.Universe(PDB_multiframe)
        self.universe2 = mda.Universe(PSF, DCD)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3
        ext = ".pdb"
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/multiwriter-test-1' + ext
        self.outfile2 = self.tmpdir.name + '/multiwriter-test-2' + ext

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        try:
            os.unlink(self.outfile2)
        except OSError:
            pass
        del self.universe, self.multiverse, self.universe2
        del self.tmpdir

    def test_write_atomselection(self):
        """Test if multiframe writer can write selected frames for an
        atomselection."""
        u = self.multiverse
        group = u.select_atoms('name CA', 'name C')
        desired_group = 56
        desired_frames = 6
        pdb = mda.Writer(self.outfile, multiframe=True, start=12, step=2)
        for ts in u.trajectory[-6:]:
            pdb.write(group)
        pdb.close()
        u2 = mda.Universe(self.outfile)
        assert_equal(len(u2.atoms), desired_group,
                     err_msg="MultiPDBWriter trajectory written for an "
                             "AtomGroup contains %d atoms, it should contain %d" % (
                                 len(u2.atoms), desired_group))

        assert_equal(len(u2.trajectory), desired_frames,
                     err_msg="MultiPDBWriter trajectory written for an "
                             "AtomGroup contains %d frames, it should have %d" % (
                                 len(u.trajectory), desired_frames))

    def test_write_all_timesteps(self):
        """
        Test write_all_timesteps() of the  multiframe writer (selected frames
        for an atomselection)
        """
        u = self.multiverse
        group = u.select_atoms('name CA', 'name C')
        desired_group = 56
        desired_frames = 6

        pdb = mda.Writer(self.outfile, multiframe=True, start=12, step=2)
        pdb.write_all_timesteps(group)
        u2 = mda.Universe(self.outfile)
        assert_equal(len(u2.atoms), desired_group,
                     err_msg="MultiPDBWriter trajectory written for an "
                             "AtomGroup contains %d atoms, it should contain %d" % (
                                 len(u2.atoms), desired_group))

        assert_equal(len(u2.trajectory), desired_frames,
                     err_msg="MultiPDBWriter trajectory written for an "
                             "AtomGroup contains %d frames, it should have %d" % (
                                 len(u.trajectory), desired_frames))

    def test_write_atoms(self):
        u = self.universe2
        W = mda.Writer(self.outfile, multiframe=True)
        # 2 frames expceted
        for ts in u.trajectory[-2:]:
            W.write(u.atoms)
        W.close()
        u0 = mda.Universe(self.outfile)
        assert_equal(u0.trajectory.n_frames,
                     2,
                     err_msg="The number of frames should be 3.")


class TestPDBReaderBig(TestCase, RefAdK):
    def setUp(self):
        self.universe = mda.Universe(PDB)
        self.prec = 6

    def tearDown(self):
        del self.universe


    def test_load_pdb(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from big PDB")
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")


    def test_selection(self):
        na = self.universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size,
                     "Atom selection of last atoms in file")


    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")


    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1,
                     "wrong number of frames")


    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0,
                     "wrong time of the frame")


    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0, "wrong frame number")


    def test_dt(self):
        """testing that accessing universe.trajectory.dt returns the default
        of 1.0 ps"""
        assert_equal(self.universe.trajectory.dt, 1.0)


    def test_coordinates(self):
        A10CA = self.universe.select_atoms('name CA')[10]
        assert_almost_equal(A10CA.position,
                            self.ref_coordinates['A10CA'],
                            self.prec,
                            err_msg="wrong coordinates for A10:CA")


    def test_distances(self):
        NTERM = self.universe.atoms.N[0]
        CTERM = self.universe.atoms.C[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec,
                            err_msg="wrong distance between M1:N and G214:C")


    def test_selection(self):
        na = self.universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size,
                     "Atom selection of last atoms in file")


    def test_unitcell(self):
        assert_array_almost_equal(
            self.universe.coord.dimensions,
            self.ref_unitcell,
            self.prec,
            err_msg="unit cell dimensions (rhombic dodecahedron), issue 60")


    def test_volume(self):
        assert_almost_equal(
            self.universe.coord.volume,
            self.ref_volume,
            0,
            err_msg="wrong volume for unitcell (rhombic dodecahedron)")

    def test_n_residues(self):
        # Should have first 10000 residues, then another 1302
        assert len(self.universe.residues) == 10000 + 1302

    def test_first_residue(self):
        # First residue is a MET, shouldn't be smushed together
        # with a water
        assert len(self.universe.residues[0].atoms) == 19


class TestIncompletePDB(TestCase):
    """Tests for Issue #396

    Reads an incomplete (but still intelligible) PDB file
    """

    def setUp(self):
        self.u = mda.Universe(INC_PDB)

    def tearDown(self):
        del self.u

    def test_natoms(self):
        assert_equal(len(self.u.atoms), 3)

    def test_coords(self):
        assert_array_almost_equal(self.u.atoms.positions,
                                  np.array([[111.2519989, 98.3730011,
                                             98.18699646],
                                            [111.20300293, 101.74199677,
                                             96.43000031], [107.60700226,
                                                            102.96800232,
                                                            96.31600189]],
                                           dtype=np.float32))

    def test_dims(self):
        assert_array_almost_equal(self.u.dimensions,
                                  np.array([216.48899841, 216.48899841,
                                            216.48899841, 90., 90., 90.],
                                           dtype=np.float32))

    def test_names(self):
        assert all(self.u.atoms.names == 'CA')

    def test_residues(self):
        assert_equal(len(self.u.residues), 3)

    def test_resnames(self):
        assert_equal(len(self.u.atoms.resnames), 3)
        assert 'VAL' in self.u.atoms.resnames
        assert 'LYS' in self.u.atoms.resnames
        assert 'PHE' in self.u.atoms.resnames

    def test_reading_trajectory(self):
        for ts in self.u.trajectory:
            pass


class TestPDBXLSerial(TestCase):
    """For Issue #446"""

    def setUp(self):
        self.u = mda.Universe(PDB_xlserial)

    def tearDown(self):
        del self.u

    def test_load(self):
        # Check that universe loads ok, should be 4 atoms
        assert len(self.u.atoms) == 4

    def test_serials(self):
        # These should be none
        assert self.u.atoms[0].id == 99998
        assert self.u.atoms[1].id == 99999
        assert self.u.atoms[2].id == 100000
        assert self.u.atoms[3].id == 100001


# Does not implement Reader.remarks, Reader.header, Reader.title,
# Reader.compounds because the PDB header data in trajectory.metadata are
# already parsed; should perhaps update the PrimitivePDBReader to do the same.


class TestPSF_CRDReader(_SingleFrameReader):
    __test__ = True
    def setUp(self):
        self.universe = mda.Universe(PSF, CRD)
        self.prec = 5  # precision in CRD (at least we are writing %9.5f)


class TestPSF_PDBReader(TestPDBReader):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3

    def test_uses_PDBReader(self):
        from MDAnalysis.coordinates.PDB import PDBReader

        assert isinstance(self.universe.trajectory, PDBReader), "failed to choose PDBReader"


def test_write_occupancies(tmpdir):
    """Tests for Issue #620 Modify occupancies, write out the file and check"""
    u = mda.Universe(PDB_small)
    u.atoms.occupancies = 0.12

    outfile = str(tmpdir.join('occ.pdb'))

    u.atoms.write(outfile)

    u2 = mda.Universe(outfile)

    assert_array_almost_equal(u2.atoms.occupancies, 0.12)


class TestWriterAlignments(object):
    @pytest.fixture(scope='class')
    def writtenstuff(self, tmpdir_factory):
        u = mda.Universe(ALIGN)
        outfile = str(tmpdir_factory.mktemp('pdb').join('nucl.pdb'))
        u.atoms.write(outfile)
        with open(outfile) as fh:
            return fh.readlines()

    def test_atomname_alignment(self, writtenstuff):
        # Our PDBWriter adds some stuff up top, so line 1 happens at [4]
        refs = ("ATOM      1  H5T",
                "ATOM      2  CA ",
                "ATOM      3 CA  ",
                "ATOM      4 H5''",)
        for written, reference in zip(writtenstuff[3:], refs):
            assert_equal(written[:16], reference)

    def test_atomtype_alignment(self, writtenstuff):
        result_line = ("ATOM      1  H5T GUA R   1       7.974   6.430   9.561"
                       "  1.00  0.00      RNAA H\n")
        assert_equal(writtenstuff[3], result_line)


@pytest.mark.parametrize('atom, refname', ((mda.coordinates.PDB.Pair('ASP', 'CA'), ' CA '),  # Regular protein carbon alpha
                                           (mda.coordinates.PDB.Pair('GLU', 'OE1'), ' OE1'),
                                           (mda.coordinates.PDB.Pair('MSE', 'SE'), 'SE  '),  # Selenium like in 4D3L
                                           (mda.coordinates.PDB.Pair('CA', 'CA'), 'CA  '),  # Calcium like in 4D3L
                                           (mda.coordinates.PDB.Pair('HDD', 'FE'), 'FE  '),  # Iron from a heme like in 1GGE
                                           (mda.coordinates.PDB.Pair('PLC', 'P'), ' P  '),  # Lipid phosphorus (1EIN)
))
def test_deduce_PDB_atom_name(atom, refname):
    # The Pair named tuple is used to mock atoms as we only need them to have a
    # ``resname`` and a ``name`` attribute.
    dummy_file = StringIO()
    name = (mda.coordinates.PDB.PDBWriter(dummy_file, n_atoms=1)
            ._deduce_PDB_atom_name(atom.name, atom.resname))
    assert_equal(name, refname)


@pytest.mark.parametrize('pdbfile', [PDB_cm, PDB_cm_bz2, PDB_cm_gz,
                                     PDB_mc, PDB_mc_bz2, PDB_mc_gz])
class TestCrystModelOrder(object):
    """Check offset based reading of pdb files

    Checks
     - len
     - seeking around

    # tests that cryst can precede or follow model header
    # allow frames to follow either of these formats:

    # Case 1 (PDB_mc)
    # MODEL
    # ...
    # ENDMDL
    # CRYST

    # Case 2 (PDB_cm)
    # CRYST
    # MODEL
    # ...
    # ENDMDL
    """
    boxsize = [80, 70, 60]
    position = [10, 20, 30]

    def test_len(self, pdbfile):
        u = mda.Universe(pdbfile)
        assert len(u.trajectory) == 3

    def test_order(self, pdbfile):
        u = mda.Universe(pdbfile)

        for ts, refbox, refpos in zip(
                u.trajectory, self.boxsize, self.position):
            assert_almost_equal(u.dimensions[0], refbox)
            assert_almost_equal(u.atoms[0].position[0], refpos)

    def test_seekaround(self, pdbfile):
        u = mda.Universe(pdbfile)

        for frame in [2, 0, 2, 1]:
            u.trajectory[frame]
            assert_almost_equal(u.dimensions[0], self.boxsize[frame])
            assert_almost_equal(u.atoms[0].position[0], self.position[frame])

    def test_rewind(self, pdbfile):
        u = mda.Universe(pdbfile)

        u.trajectory[2]
        u.trajectory.rewind()
        assert_almost_equal(u.dimensions[0], self.boxsize[0])
        assert_almost_equal(u.atoms[0].position[0], self.position[0])


def test_standalone_pdb():
    # check that PDBReader works without n_atoms kwarg
    r = mda.coordinates.PDB.PDBReader(PDB_cm)

    assert r.n_atoms == 4


def test_write_pdb_zero_atoms():
    # issue 1083
    u = make_Universe(trajectory=True)

    with tempdir.TempDir() as tmpdir:
        outfile = tmpdir + '/out.pdb'

        ag = u.atoms[:0]  # empty ag

        with mda.Writer(outfile, ag.n_atoms) as w:
            with pytest.raises(IndexError):
                w.write(ag)
