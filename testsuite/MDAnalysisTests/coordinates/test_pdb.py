import MDAnalysis as mda
import numpy as np
import os
from six.moves import zip

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_, dec,
                           assert_array_almost_equal,
                           assert_almost_equal, assert_raises)
import tempdir
from unittest import TestCase

from MDAnalysisTests.coordinates.reference import (RefAdKSmall, Ref4e43,
                                                   RefAdK)
from MDAnalysisTests.coordinates.base import _SingleFrameReader
from MDAnalysisTests.datafiles import (PDB, PDB_small, PDB_multiframe,
                                       XPDB_small, PSF, DCD, CONECT, CRD,
                                       INC_PDB, PDB_xlserial,
                                       NUCL)
from MDAnalysisTests.plugins.knownfailure import knownfailure


class TestPDBReader(_SingleFrameReader):
    def setUp(self):
        # use permissive=False instead of changing the global flag as this
        # can lead to race conditions when testing in parallel
        self.universe = mda.Universe(RefAdKSmall.filename, permissive=False)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3

    def test_uses_Biopython(self):
        from MDAnalysis.coordinates.PDB import PDBReader

        assert_(isinstance(self.universe.trajectory, PDBReader),
                "failed to choose Biopython PDBReader")

    @knownfailure("Biopython PDB reader does not parse CRYST1", AssertionError)
    def test_dimensions(self):
        assert_almost_equal(
            self.universe.trajectory.ts.dimensions, RefAdKSmall.ref_unitcell,
            self.prec,
            "Biopython reader failed to get unitcell dimensions from CRYST1")


class _PDBMetadata(TestCase, Ref4e43):
    permissive = True

    def setUp(self):
        self.universe = mda.Universe(self.filename, permissive=self.permissive)

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


class TestPrimitivePDBReader_Metadata(_PDBMetadata):
    permissive = True


class TestPrimitivePDBReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(PDB_small, permissive=True)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3

    def test_missing_natoms(self):
        from MDAnalysis.coordinates.PDB import PrimitivePDBReader

        assert_raises(ValueError, PrimitivePDBReader, 'something.pdb')

    def test_wrong_natoms(self):
        from MDAnalysis.coordinates.PDB import PrimitivePDBReader

        assert_raises(ValueError, PrimitivePDBReader, PDB_small, n_atoms=4000)


class TestExtendedPDBReader(_SingleFrameReader):
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


class TestPSF_PrimitivePDBReader(TestPrimitivePDBReader):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small, permissive=True)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3

    def test_dimensions(self):
        assert_almost_equal(self.universe.trajectory.ts.dimensions,
                            RefAdKSmall.ref_unitcell, self.prec,
                            "Primitive PDB reader failed to get unitcell "
                            "dimensions from CRYST1")


class TestPrimitivePDBWriter(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small, permissive=True)
        self.universe2 = mda.Universe(PSF, DCD, permissive=True)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3
        ext = ".pdb"
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/primitive-pdb-writer' + ext

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe, self.universe2
        del self.tmpdir

    def test_writer(self):
        "Test writing from a single frame PDB file to a PDB file." ""
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(PSF, self.outfile, permissive=True)
        assert_almost_equal(u.atoms.coordinates(),
                            self.universe.atoms.coordinates(), self.prec,
                            err_msg="Writing PDB file with PrimitivePDBWriter "
                            "does not reproduce original coordinates")

    @attr('issue')
    def test_write_single_frame_Writer(self):
        """Test writing a single frame from a DCD trajectory to a PDB using
        MDAnalysis.Writer (Issue 105)"""
        u = self.universe2
        W = mda.Writer(self.outfile)
        u.trajectory[50]
        W.write(u.select_atoms('all'))
        W.close()
        u2 = mda.Universe(self.outfile)
        assert_equal(u2.trajectory.n_frames,
                     1,
                     err_msg="The number of frames should be 1.")

    @attr('issue')
    def test_write_single_frame_AtomGroup(self):
        """Test writing a single frame from a DCD trajectory to a PDB using
        AtomGroup.write() (Issue 105)"""
        u = self.universe2
        u.trajectory[50]
        u.atoms.write(self.outfile)
        u2 = mda.Universe(PSF, self.outfile)
        assert_equal(u2.trajectory.n_frames,
                     1,
                     err_msg="Output PDB should only contain a single frame")
        assert_almost_equal(u2.atoms.coordinates(), u.atoms.coordinates(),
                            self.prec, err_msg="Written coordinates do not "
                            "agree with original coordinates from frame %d" %
                            u.trajectory.frame)

    @attr('issue')
    def test_check_coordinate_limits_min(self):
        """Test that illegal PDB coordinates (x <= -999.9995 A) are caught
        with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(PSF, PDB_small, permissive=True)
        u.atoms[2000].pos[1] = -999.9995
        assert_raises(ValueError, u.atoms.write, self.outfile)
        del u

    @attr('issue')
    def test_check_coordinate_limits_max(self):
        """Test that illegal PDB coordinates (x > 9999.9995 A) are caught
        with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(PSF, PDB_small, permissive=True)
        # OB: 9999.99951 is not caught by '<=' ?!?
        u.atoms[1000].pos[1] = 9999.9996
        assert_raises(ValueError, u.atoms.write, self.outfile)
        del u


class TestMultiPDBReader(TestCase):
    def setUp(self):
        self.multiverse = mda.Universe(PDB_multiframe,
                                       permissive=True,
                                       guess_bonds=True)
        self.multiverse.build_topology()
        self.conect = mda.Universe(CONECT, guess_bonds=True)
        self.conect.build_topology()

    def tearDown(self):
        del self.multiverse
        del self.conect

    @attr('slow')
    def test_n_frames(self):
        assert_equal(self.multiverse.trajectory.n_frames, 24,
                     "Wrong number of frames read from PDB muliple model file")

    @attr('slow')
    def test_n_atoms_frame(self):
        u = self.multiverse
        desired = 392
        for frame in u.trajectory:
            assert_equal(len(u.atoms), desired, err_msg="The number of atoms "
                         "in the Universe (%d) does not" " match the number "
                         "of atoms in the test case (%d) at frame %d" % (
                             len(u.atoms), desired, u.trajectory.frame))

    @attr('slow')
    def test_rewind(self):
        u = self.multiverse
        u.trajectory[11]
        assert_equal(u.trajectory.ts.frame, 11,
                     "Failed to forward to 11th frame (frame index 11)")
        u.trajectory.rewind()
        assert_equal(u.trajectory.ts.frame, 0,
                     "Failed to rewind to 0th frame (frame index 0)")

    @attr('slow')
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

    @attr('slow')
    def test_slice_iteration(self):
        u = self.multiverse
        frames = []
        for ts in u.trajectory[4:-2:4]:
            frames.append(ts.frame)
        assert_equal(np.array(frames),
                     np.arange(u.trajectory.n_frames)[4:-2:4],
                     err_msg="slicing did not produce the expected frames")

    @attr('slow')
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

    @attr('slow')
    def test_conect_bonds_all(self):
        conect = self.conect
        assert_equal(len(conect.atoms), 1890)
        assert_equal(len(conect.bonds), 1922)

        with tempdir.in_tempdir():
            try:
                outfile = 'pdb-connect-bonds.pdb'
                self.conect.atoms.write(outfile, bonds="all")
                u2 = mda.Universe(outfile, guess_bonds=True)
            finally:
                os.unlink(outfile)
            assert_equal(len(u2.atoms), 1890)
            assert_equal(len([b for b in u2.bonds if not b.is_guessed]), 1922)

        #assert_equal(len([b for b in conect.bonds if not b.is_guessed]), 1922)

    @attr('slow')
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
                     "is %d" % (len(u._topology['bonds']), len(desired)))


class TestMultiPDBWriter(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small, permissive=True)
        self.multiverse = mda.Universe(PDB_multiframe, permissive=True)
        self.universe2 = mda.Universe(PSF, DCD, permissive=True)
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

    @attr('slow')
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

    @attr('slow')
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

    @attr('slow')
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

    @dec.slow
    def test_load_pdb(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from big PDB")
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    @dec.slow
    def test_selection(self):
        na = self.universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size,
                     "Atom selection of last atoms in file")

    @dec.slow
    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")

    @dec.slow
    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1,
                     "wrong number of frames")

    @dec.slow
    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0,
                     "wrong time of the frame")

    @dec.slow
    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0, "wrong frame number")

    @dec.slow
    def test_dt(self):
        """testing that accessing universe.trajectory.dt returns the default
        of 1.0 ps"""
        assert_equal(self.universe.trajectory.dt, 1.0)

    @dec.slow
    def test_coordinates(self):
        A10CA = self.universe.SYSTEM.CA[10]
        assert_almost_equal(A10CA.pos,
                            self.ref_coordinates['A10CA'],
                            self.prec,
                            err_msg="wrong coordinates for A10:CA")

    @dec.slow
    def test_distances(self):
        NTERM = self.universe.SYSTEM.N[0]
        CTERM = self.universe.SYSTEM.C[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec,
                            err_msg="wrong distance between M1:N and G214:C")

    @dec.slow
    def test_selection(self):
        na = self.universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size,
                     "Atom selection of last atoms in file")

    @dec.slow
    @attr('issue')
    def test_unitcell(self):
        assert_array_almost_equal(
            self.universe.coord.dimensions,
            self.ref_unitcell,
            self.prec,
            err_msg="unit cell dimensions (rhombic dodecahedron), issue 60")

    @dec.slow
    def test_volume(self):
        assert_almost_equal(
            self.universe.coord.volume,
            self.ref_volume,
            0,
            err_msg="wrong volume for unitcell (rhombic dodecahedron)")

    def test_n_residues(self):
        # Should have first 10000 residues, then another 1302
        assert_(len(self.universe.residues) == 10000 + 1302)

    def test_first_residue(self):
        # First residue is a MET, shouldn't be smushed together
        # with a water
        assert_(len(self.universe.residues[0]) == 19)


class TestIncompletePDB(object):
    """Tests for Issue #396

    Reads an incomplete (but still intelligible) PDB file
    """

    def setUp(self):
        self.u = mda.Universe(INC_PDB)

    def tearDown(self):
        del self.u

    def test_natoms(self):
        assert len(self.u.atoms) == 3

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
        assert len(self.u.residues) == 3

    def test_resnames(self):
        assert len(self.u.atoms.resnames) == 3
        assert 'VAL' in self.u.atoms.resnames
        assert 'LYS' in self.u.atoms.resnames
        assert 'PHE' in self.u.atoms.resnames

    def test_reading_trajectory(self):
        for ts in self.u.trajectory:
            pass

    def test_occupancy(self):
        occupancies = self.u.atoms.occupancies
        assert_array_almost_equal(occupancies, np.ones(len(occupancies)))

    def test_set_occupancy(self):
        for atom in self.u.atoms:
            atom.occupancy = 0
        assert_almost_equal(self.u.atoms.occupancies,
                            np.zeros(self.u.atoms.n_atoms))

    def test_set_occupancies(self):
        self.u.atoms.occupancies = 0.0
        assert_almost_equal(self.u.atoms.occupancies,
                            np.zeros(self.u.atoms.n_atoms))


class TestPDBXLSerial(object):
    """For Issue #446"""

    def setUp(self):
        self.u = mda.Universe(PDB_xlserial)

    def tearDown(self):
        del self.u

    def test_load(self):
        # Check that universe loads ok, should be 4 atoms
        assert_(len(self.u.atoms) == 4)

    def test_serials(self):
        # These should be none
        assert_(self.u.atoms[0].serial == 99998)
        assert_(self.u.atoms[1].serial == 99999)
        assert_(self.u.atoms[2].serial is None)
        assert_(self.u.atoms[3].serial is None)


# Does not implement Reader.remarks, Reader.header, Reader.title,
# Reader.compounds because the PDB header data in trajectory.metadata are
# already parsed; should perhaps update the PrimitivePDBReader to do the same.
# [orbeckst] class TestPDBReader_Metadata(_PDBMetadata): permissive = False


class TestPSF_CRDReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(PSF, CRD)
        self.prec = 5  # precision in CRD (at least we are writing %9.5f)


class TestPSF_PDBReader(TestPDBReader):
    def setUp(self):
        # mda.core.flags['permissive_pdb_reader'] = False
        self.universe = mda.Universe(PSF, PDB_small, permissive=False)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3

    def test_uses_Biopython(self):
        from MDAnalysis.coordinates.PDB import PDBReader

        assert_(isinstance(self.universe.trajectory, PDBReader),
                "failed to choose Biopython PDBReader")


class TestPDBWriterOccupancies(object):
    """Tests for Issue #620"""
    def setUp(self):
        self.tempdir = tempdir.TempDir()
        self.outfile = self.tempdir.name + '/occ.pdb'

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.tempdir
        del self.outfile

    def test_write_occupancies(self):
        """Modify occupancies, write out the file and check"""
        u = mda.Universe(PDB_small)

        u.atoms.occupancies = 0.12

        u.atoms.write(self.outfile)

        u2 = mda.Universe(self.outfile)

        assert_(all(u2.atoms.occupancies == 0.12))

def test_writer_alignments():
    u = mda.Universe(NUCL)

    tmpdir = tempdir.TempDir()
    outfile = tmpdir.name + '/nucl.pdb'

    u.atoms.write(outfile)

    writtenstuff = open(outfile, 'r').readlines()

    # Our PDBWriter adds some stuff up top, so line 1 happens at [4]
    assert_(writtenstuff[4].startswith(
        "ATOM      1  H5T GUA"))
    assert_(writtenstuff[8].startswith(
        "ATOM      5 H5'' GUA"))

    del tmpdir
    del outfile
