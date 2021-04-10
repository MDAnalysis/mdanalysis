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
import os
from io import StringIO

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysisTests import make_Universe
from MDAnalysisTests.coordinates.base import _SingleFrameReader
from MDAnalysisTests.coordinates.reference import (RefAdKSmall,
                                                   RefAdK)
from MDAnalysisTests.datafiles import (PDB, PDB_small, PDB_multiframe,
                                       PDB_full,
                                       XPDB_small, PSF, DCD, CONECT, CRD,
                                       INC_PDB, PDB_xlserial, ALIGN, ENT,
                                       PDB_cm, PDB_cm_gz, PDB_cm_bz2,
                                       PDB_mc, PDB_mc_gz, PDB_mc_bz2,
                                       PDB_CRYOEM_BOX, MMTF_NOCRYST,
                                       PDB_HOLE, mol2_molecule)
from numpy.testing import (assert_equal,
                           assert_array_almost_equal,
                           assert_almost_equal)

IGNORE_NO_INFORMATION_WARNING = 'ignore:Found no information for attr:UserWarning'


@pytest.fixture
def dummy_universe_without_elements():
    n_atoms = 5
    u = make_Universe(size=(n_atoms, 1, 1), trajectory=True)
    u.add_TopologyAttr('resnames', ['RES'])
    u.add_TopologyAttr('names', ['C1', 'O2', 'N3', 'S4', 'NA'])
    u.dimensions = [42, 42, 42, 90, 90, 90]
    return u


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


class TestPDBMetadata(object):
    header = 'HYDROLASE                               11-MAR-12   4E43'
    title = ['HIV PROTEASE (PR) DIMER WITH ACETATE IN EXO SITE AND PEPTIDE '
             'IN ACTIVE', '2 SITE']
    compnd = ['MOL_ID: 1;',
              '2 MOLECULE: PROTEASE;',
              '3 CHAIN: A, B;',
              '4 ENGINEERED: YES;',
              '5 MUTATION: YES;',
              '6 MOL_ID: 2;',
              '7 MOLECULE: RANDOM PEPTIDE;',
              '8 CHAIN: C;',
              '9 ENGINEERED: YES;',
              '10 OTHER_DETAILS: UNKNOWN IMPURITY', ]
    num_remarks = 333
    # only first 5 remarks for comparison
    nmax_remarks = 5
    remarks = [
        '2',
        '2 RESOLUTION.    1.54 ANGSTROMS.',
        '3',
        '3 REFINEMENT.',
        '3   PROGRAM     : REFMAC 5.5.0110',
    ]

    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return mda.Universe(PDB_full)

    def test_HEADER(self, universe):
        assert_equal(universe.trajectory.header,
                     self.header,
                     err_msg="HEADER record not correctly parsed")

    def test_TITLE(self, universe):
        try:
            title = universe.trajectory.title
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

    def test_COMPND(self, universe):
        try:
            compound = universe.trajectory.compound
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

    def test_REMARK(self, universe):
        try:
            remarks = universe.trajectory.remarks
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
    def universe3(self):
        return mda.Universe(PDB)

    @pytest.fixture
    def universe4(self):
        return mda.Universe(PDB_HOLE)

    @pytest.fixture
    def universe5(self):
        return mda.Universe(mol2_molecule)

    @pytest.fixture(params=[
            [PDB_CRYOEM_BOX, np.zeros(6)],
            [MMTF_NOCRYST, None]
        ])
    def universe_and_expected_dims(self, request):
        """
        File with meaningless CRYST1 record and expected dimensions.

        Notes
        -----
        This will need to be made consistent, see Issue #2698
        """
        filein = request.param[0]
        expected_dims = request.param[1]

        return mda.Universe(filein), expected_dims

    @pytest.fixture
    def outfile(self, tmpdir):
        return str(tmpdir.mkdir("PDBWriter").join('primitive-pdb-writer' + self.ext))

    @pytest.fixture
    def u_no_ids(self):
        # The test universe does not have atom ids, but it has everything
        # else the PDB writer expects to avoid issuing warnings.
        universe = make_Universe(
            [
                'names', 'resids', 'resnames', 'altLocs',
                'segids', 'occupancies', 'tempfactors',
            ],
            trajectory=True,
        )
        universe.add_TopologyAttr('icodes', [' '] * len(universe.residues))
        universe.add_TopologyAttr('record_types', ['ATOM'] * len(universe.atoms))
        universe.dimensions = [10, 10, 10, 90, 90, 90]
        return universe

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
        expected = np.array(['X'] * u_no_names.atoms.n_atoms)
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

    def test_write_nodims(self, universe_and_expected_dims, outfile):
        """
        Test :code:`PDBWriter` for universe without cell dimensions.

        Notes
        -----
        Test fix for Issue #2679.
        """

        u, expected_dims = universe_and_expected_dims

        # See Issue #2698
        if expected_dims is None:
            assert u.dimensions is None
        else:
            assert np.allclose(u.dimensions, expected_dims)

        expected_msg = "Unit cell dimensions not found. CRYST1 record set to unitary values."

        with pytest.warns(UserWarning, match=expected_msg):
            u.atoms.write(outfile)

        with pytest.warns(UserWarning, match="Unit cell dimensions will be set to zeros."):
            uout = mda.Universe(outfile)

        assert_almost_equal(
            uout.dimensions, np.zeros(6),
            self.prec,
            err_msg="Problem with default box."
        )

        assert_equal(
            uout.trajectory.n_frames, 1,
            err_msg="Output PDB should only contain a single frame"
        )

        assert_almost_equal(
            u.atoms.positions, uout.atoms.positions,
            self.prec,
            err_msg="Written coordinates do not "
                    "agree with original coordinates from frame %d" %
                    u.trajectory.frame
        )


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

    @pytest.mark.parametrize("bad_chainid",
                             ['@', '', 'AA'])
    def test_chainid_validated(self, universe3, outfile, bad_chainid):
        """
        Check that an atom's chainID is set to 'X' if the chainID
        does not confirm to standards (issue #2224)
        """
        default_id = 'X'
        u = universe3
        u.atoms.chainIDs = bad_chainid
        u.atoms.write(outfile)
        u_pdb = mda.Universe(outfile)
        assert_equal(u_pdb.segments.chainIDs[0][0], default_id)

    def test_stringio_outofrange(self, universe3):
        """
        Check that when StringIO is used, the correct out-of-range error for
        coordinates is raised (instead of failing trying to remove StringIO
        as a file).
        """

        u = universe3

        u.atoms.translate([-9999, -9999, -9999])

        outstring = StringIO()

        errmsg = "PDB files must have coordinate values between"

        with pytest.raises(ValueError, match=errmsg):
            with mda.coordinates.PDB.PDBWriter(outstring) as writer:
                writer.write(u.atoms)

    def test_hetatm_written(self, universe4, tmpdir, outfile):
        """
        Checks that HETATM record types are written.
        """

        u = universe4
        u_hetatms = u.select_atoms("resname ETA and record_type HETATM")
        assert_equal(len(u_hetatms), 8)

        u.atoms.write(outfile)
        written = mda.Universe(outfile)
        written_atoms = written.select_atoms("resname ETA and "
                                             "record_type HETATM")

        assert len(u_hetatms) == len(written_atoms), \
            "mismatched HETATM number"
        assert_almost_equal(u_hetatms.atoms.positions,
                            written_atoms.atoms.positions)

    def test_default_atom_record_type_written(self, universe5, tmpdir,
                                              outfile):
        """
        Checks that ATOM record types are written when there is no
        record_type attribute.
        """

        u = universe5

        expected_msg = ("Found no information for attr: "
                        "'record_types' Using default value of 'ATOM'")

        with pytest.warns(UserWarning, match=expected_msg):
            u.atoms.write(outfile)

        written = mda.Universe(outfile)
        assert len(u.atoms) == len(written.atoms), \
            "mismatched number of atoms"

        atms = written.select_atoms("record_type ATOM")
        assert len(atms.atoms) == len(u.atoms), \
            "mismatched ATOM number"

        hetatms = written.select_atoms("record_type HETATM")
        assert len(hetatms.atoms) == 0, "mismatched HETATM number"

    def test_abnormal_record_type(self, universe5, tmpdir, outfile):
        """
        Checks whether KeyError is raised when record type is
        neither ATOM or HETATM.
        """
        u = universe5
        u.add_TopologyAttr('record_type', ['ABNORM']*len(u.atoms))

        expected_msg = ("Found ABNORM for the record type, but only "
                        "allowed types are ATOM or HETATM")

        with pytest.raises(ValueError, match=expected_msg):
            u.atoms.write(outfile)

    @pytest.mark.filterwarnings(IGNORE_NO_INFORMATION_WARNING)
    def test_no_reindex(self, universe, outfile):
        """
        When setting the `reindex` keyword to False, the atom are
        not reindexed.
        """
        universe.atoms.ids = universe.atoms.ids + 23
        universe.atoms.write(outfile, reindex=False)
        read_universe = mda.Universe(outfile)
        assert np.all(read_universe.atoms.ids == universe.atoms.ids)

    @pytest.mark.filterwarnings(IGNORE_NO_INFORMATION_WARNING)
    def test_no_reindex_bonds(self, universe, outfile):
        """
        When setting the `reindex` keyword to False, the connect
        record match the non-reindexed atoms.
        """
        universe.atoms.ids = universe.atoms.ids + 23
        universe.atoms.write(outfile, reindex=False, bonds='all')
        with open(outfile) as infile:
            for line in infile:
                if line.startswith('CONECT'):
                    assert line.strip() == "CONECT   23   24   25   26   27"
                    break
            else:
                raise AssertError('No CONECT record fond in the output.')

    @pytest.mark.filterwarnings(IGNORE_NO_INFORMATION_WARNING)
    def test_reindex(self, universe, outfile):
        """
        When setting the `reindex` keyword to True, the atom are
        reindexed.
        """
        universe.atoms.ids = universe.atoms.ids + 23
        universe.atoms.write(outfile, reindex=True)
        read_universe = mda.Universe(outfile)
        # AG.ids is 1-based, while AG.indices is 0-based, hence the +1
        assert np.all(read_universe.atoms.ids == universe.atoms.indices + 1)

    def test_no_reindex_missing_ids(self, u_no_ids, outfile):
        """
        When setting `reindex` to False, if there is no AG.ids,
        then an exception is raised.
        """
        # Making sure AG.ids is indeed missing
        assert not hasattr(u_no_ids.atoms, 'ids')
        with pytest.raises(mda.exceptions.NoDataError):
            u_no_ids.atoms.write(outfile, reindex=False)


class TestMultiPDBReader(object):
    @staticmethod
    @pytest.fixture(scope='class')
    def multiverse():
        return mda.Universe(PDB_multiframe, guess_bonds=True)

    @staticmethod
    @pytest.fixture(scope='class')
    def conect():
        return mda.Universe(CONECT, guess_bonds=True)

    def test_n_frames(self, multiverse):
        assert_equal(multiverse.trajectory.n_frames, 24,
                     "Wrong number of frames read from PDB muliple model file")

    def test_n_atoms_frame(self, multiverse):
        u = multiverse
        desired = 392
        for frame in u.trajectory:
            assert_equal(len(u.atoms), desired, err_msg="The number of atoms "
                                                        "in the Universe (%d) does not" " match the number "
                                                        "of atoms in the test case (%d) at frame %d" % (
                                                            len(u.atoms), desired, u.trajectory.frame))

    def test_rewind(self, multiverse):
        u = multiverse
        u.trajectory[11]
        assert_equal(u.trajectory.ts.frame, 11,
                     "Failed to forward to 11th frame (frame index 11)")
        u.trajectory.rewind()
        assert_equal(u.trajectory.ts.frame, 0,
                     "Failed to rewind to 0th frame (frame index 0)")

    def test_iteration(self, multiverse):
        u = multiverse
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

    def test_slice_iteration(self, multiverse):
        u = multiverse
        frames = []
        for ts in u.trajectory[4:-2:4]:
            frames.append(ts.frame)
        assert_equal(np.array(frames),
                     np.arange(u.trajectory.n_frames)[4:-2:4],
                     err_msg="slicing did not produce the expected frames")

    def test_conect_bonds_conect(self, tmpdir, conect):
        assert_equal(len(conect.atoms), 1890)
        assert_equal(len(conect.bonds), 1922)

        outfile = str(tmpdir.join('test-pdb-hbonds.pdb'))
        conect.atoms.write(outfile, bonds="conect")
        u1 = mda.Universe(outfile, guess_bonds=True)

        assert_equal(len(u1.atoms), 1890)
        assert_equal(len(u1.bonds), 1922)

    def test_numconnections(self, multiverse):
        u = multiverse

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

        conect = helper(u.atoms, [b for b in u.bonds if not b.is_guessed])
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


def test_write_bonds_partial(tmpdir):
    u = mda.Universe(CONECT)
    # grab all atoms with bonds
    ag = (u.atoms.bonds.atom1 + u.atoms.bonds.atom2).unique

    outfile = os.path.join(str(tmpdir), 'test.pdb')
    ag.write(outfile)

    u2 = mda.Universe(outfile)

    assert len(u2.atoms.bonds) > 0
    # check bonding is correct in new universe
    for a_ref, atom in zip(ag, u2.atoms):
        assert len(a_ref.bonds) == len(atom.bonds)


class TestMultiPDBWriter(object):
    # 3 decimals in PDB spec
    # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
    prec = 3

    @staticmethod
    @pytest.fixture
    def universe():
        return mda.Universe(PSF, PDB_small)

    @staticmethod
    @pytest.fixture
    def multiverse():
        return mda.Universe(PDB_multiframe)

    @staticmethod
    @pytest.fixture
    def universe2():
        return mda.Universe(PSF, DCD)

    @staticmethod
    @pytest.fixture
    def outfile(tmpdir):
        return os.path.join(str(tmpdir), 'multiwriter-test-1.pdb')

    def test_write_atomselection(self, multiverse, outfile):
        """Test if multiframe writer can write selected frames for an
        atomselection."""
        u = multiverse
        group = u.select_atoms('name CA', 'name C')
        desired_group = 56
        desired_frames = 6
        pdb = mda.Writer(outfile, multiframe=True, start=12, step=2)
        for ts in u.trajectory[-6:]:
            pdb.write(group)
        pdb.close()
        u2 = mda.Universe(outfile)
        assert_equal(len(u2.atoms), desired_group,
                     err_msg="MultiPDBWriter trajectory written for an "
                             "AtomGroup contains %d atoms, it should contain %d" % (
                                 len(u2.atoms), desired_group))

        assert_equal(len(u2.trajectory), desired_frames,
                     err_msg="MultiPDBWriter trajectory written for an "
                             "AtomGroup contains %d frames, it should have %d" % (
                                 len(u.trajectory), desired_frames))

    def test_write_all_timesteps(self, multiverse, outfile):
        """
        Test write_all_timesteps() of the  multiframe writer (selected frames
        for an atomselection)
        """
        u = multiverse
        group = u.select_atoms('name CA', 'name C')
        desired_group = 56
        desired_frames = 6

        with mda.Writer(outfile, multiframe=True, start=12, step=2) as W:
            W.write_all_timesteps(group)
        u2 = mda.Universe(outfile)
        assert_equal(len(u2.atoms), desired_group,
                     err_msg="MultiPDBWriter trajectory written for an "
                             "AtomGroup contains %d atoms, it should contain %d" % (
                                 len(u2.atoms), desired_group))

        assert_equal(len(u2.trajectory), desired_frames,
                     err_msg="MultiPDBWriter trajectory written for an "
                             "AtomGroup contains %d frames, it should have %d" % (
                                 len(u.trajectory), desired_frames))

        with open(outfile, "r") as f:
            lines = f.read()
            assert lines.count("CONECT") == 2  # Expected two CONECT records

    def test_write_loop(self, multiverse, outfile):
        """
        Test write() in a loop with the multiframe writer (selected frames
        for an atomselection)
        """
        u = multiverse
        group = u.select_atoms('name CA', 'name C')
        desired_group = 56
        desired_frames = 6

        with mda.Writer(outfile, multiframe=True) as W:
            for ts in u.trajectory[12::2]:
                W.write(group)

        u2 = mda.Universe(outfile)
        assert_equal(len(u2.atoms), desired_group,
                     err_msg="MultiPDBWriter trajectory written for an "
                             f"AtomGroup contains {len(u2.atoms)} atoms, "
                             f"it should contain {desired_group}")

        assert_equal(len(u2.trajectory), desired_frames,
                     err_msg="MultiPDBWriter trajectory written for an "
                             f"AtomGroup contains {len(u.trajectory)} "
                             f"frames, it should have {desired_frames}")

        with open(outfile, "r") as f:
            lines = f.read()

            # Expected only two CONECT records
            assert lines.count("CONECT") == 2

    def test_write_atoms(self, universe2, outfile):
        u = universe2
        with mda.Writer(outfile, multiframe=True) as W:
            # 2 frames expected
            for ts in u.trajectory[-2:]:
                W.write(u.atoms)

        u0 = mda.Universe(outfile)
        assert_equal(u0.trajectory.n_frames,
                     2,
                     err_msg="The number of frames should be 2.")


class TestPDBReaderBig(RefAdK):
    prec = 6

    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return mda.Universe(PDB)

    def test_load_pdb(self, universe):
        U = universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from big PDB")
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_selection(self, universe):
        na = universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size,
                     "Atom selection of last atoms in file")

    def test_n_atoms(self, universe):
        assert_equal(universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")

    def test_n_frames(self, universe):
        assert_equal(universe.trajectory.n_frames, 1,
                     "wrong number of frames")

    def test_time(self, universe):
        assert_equal(universe.trajectory.time, 0.0,
                     "wrong time of the frame")

    def test_frame(self, universe):
        assert_equal(universe.trajectory.frame, 0, "wrong frame number")

    def test_dt(self, universe):
        """testing that accessing universe.trajectory.dt returns the default
        of 1.0 ps"""
        assert_equal(universe.trajectory.dt, 1.0)

    def test_coordinates(self, universe):
        A10CA = universe.select_atoms('name CA')[10]
        assert_almost_equal(A10CA.position,
                            self.ref_coordinates['A10CA'],
                            self.prec,
                            err_msg="wrong coordinates for A10:CA")

    def test_distances(self, universe):
        NTERM = universe.select_atoms('name N')[0]
        CTERM = universe.select_atoms('name C')[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec,
                            err_msg="wrong distance between M1:N and G214:C")

    def test_selection(self, universe):
        na = universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size,
                     "Atom selection of last atoms in file")

    def test_unitcell(self, universe):
        assert_array_almost_equal(
            universe.dimensions,
            self.ref_unitcell,
            self.prec,
            err_msg="unit cell dimensions (rhombic dodecahedron), issue 60")

    def test_volume(self, universe):
        assert_almost_equal(
            universe.coord.volume,
            self.ref_volume,
            0,
            err_msg="wrong volume for unitcell (rhombic dodecahedron)")

    def test_n_residues(self, universe):
        # Should have first 10000 residues, then another 1302
        assert len(universe.residues) == 10000 + 1302

    def test_first_residue(self, universe):
        # First residue is a MET, shouldn't be smushed together
        # with a water
        assert len(universe.residues[0].atoms) == 19


class TestIncompletePDB(object):
    """Tests for Issue #396

    Reads an incomplete (but still intelligible) PDB file
    """
    @staticmethod
    @pytest.fixture(scope='class')
    def u():
        return mda.Universe(INC_PDB)

    def test_natoms(self, u):
        assert_equal(len(u.atoms), 3)

    def test_coords(self, u):
        assert_array_almost_equal(u.atoms.positions,
                                  np.array([[111.2519989, 98.3730011,
                                             98.18699646],
                                            [111.20300293, 101.74199677,
                                             96.43000031], [107.60700226,
                                                            102.96800232,
                                                            96.31600189]],
                                           dtype=np.float32))

    def test_dims(self, u):
        assert_array_almost_equal(u.dimensions,
                                  np.array([216.48899841, 216.48899841,
                                            216.48899841, 90., 90., 90.],
                                           dtype=np.float32))

    def test_names(self, u):
        assert all(u.atoms.names == 'CA')

    def test_residues(self, u):
        assert_equal(len(u.residues), 3)

    def test_resnames(self, u):
        assert_equal(len(u.atoms.resnames), 3)
        assert 'VAL' in u.atoms.resnames
        assert 'LYS' in u.atoms.resnames
        assert 'PHE' in u.atoms.resnames

    def test_reading_trajectory(self, u):
        counter = 0
        for ts in u.trajectory:
            counter += 1
        assert counter == 2


class TestPDBXLSerial(object):
    """For Issue #446"""
    @staticmethod
    @pytest.fixture(scope='class')
    def u():
        return mda.Universe(PDB_xlserial)

    def test_load(self, u):
        # Check that universe loads ok, should be 4 atoms
        assert len(u.atoms) == 4

    def test_serials(self, u):
        # These should be none
        assert u.atoms[0].id == 99998
        assert u.atoms[1].id == 99999
        assert u.atoms[2].id == 100000
        assert u.atoms[3].id == 100001


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
        # Our PDBWriter adds some stuff up top, so line 1 happens at [9]
        refs = ("ATOM      1  H5T",
                "ATOM      2  CA ",
                "ATOM      3 CA  ",
                "ATOM      4 H5''",)

        for written, reference in zip(writtenstuff[9:], refs):
            assert_equal(written[:16], reference)

    def test_atomtype_alignment(self, writtenstuff):
        result_line = ("ATOM      1  H5T GUA X   1       7.974   6.430   9.561"
                       "  1.00  0.00      RNAA  \n")
        assert_equal(writtenstuff[9], result_line)


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


def test_write_pdb_zero_atoms(tmpdir):
    # issue 1083
    u = make_Universe(trajectory=True)

    with tmpdir.as_cwd():
        outfile = 'out.pdb'

        ag = u.atoms[:0]  # empty ag

        with mda.Writer(outfile, ag.n_atoms) as w:
            with pytest.raises(IndexError):
                w.write(ag)


def test_atom_not_match(tmpdir):
    # issue 1998
    outfile = str(tmpdir.mkdir("PDBReader").join('test_atom_not_match' + ".pdb"))
    u = mda.Universe(PSF, DCD)
    # select two groups of atoms
    protein = u.select_atoms("protein and name CA")
    atoms = u.select_atoms(
        'resid 1 or resid 10 or resid 100 or resid 1000 or resid 10000')
    with mda.Writer(outfile, multiframe=True, n_atoms=10) as pdb:
        # write these two groups of atoms to pdb
        # Then the n_atoms will not match
        pdb.write(protein)
        pdb.write(atoms)
    reader = mda.coordinates.PDB.PDBReader(outfile)
    with pytest.raises(ValueError) as excinfo:
        reader._read_frame(1)
    assert 'Inconsistency in file' in str(excinfo.value)


def test_partially_missing_cryst():
    # issue 2252
    raw = open(INC_PDB, 'r').readlines()
    # mangle the cryst lines so that only box angles are left
    # this mimics '6edu' from PDB
    raw = [line if not line.startswith('CRYST')
           else line[:6] + ' ' * 28 + line[34:]
           for line in raw]

    with pytest.warns(UserWarning):
        u = mda.Universe(StringIO('\n'.join(raw)), format='PDB')

    assert len(u.atoms) == 3
    assert len(u.trajectory) == 2
    assert_array_almost_equal(u.dimensions, 0.0)


@pytest.mark.filterwarnings(IGNORE_NO_INFORMATION_WARNING)
def test_write_no_atoms_elements(dummy_universe_without_elements):
    """
    If no element symbols are provided, the PDB writer guesses.
    """
    destination = StringIO()
    with mda.coordinates.PDB.PDBWriter(destination) as writer:
        writer.write(dummy_universe_without_elements.atoms)
        content = destination.getvalue()
    element_symbols = [
        line[76:78].strip()
        for line in content.splitlines()
        if line[:6] == 'ATOM  '
    ]
    expectation = ['', '', '', '', '']
    assert element_symbols == expectation


@pytest.mark.filterwarnings(IGNORE_NO_INFORMATION_WARNING)
def test_write_atom_elements(dummy_universe_without_elements):
    """
    If element symbols are provided, they are used when writing the file.

    See `Issue 2423 <https://github.com/MDAnalysis/mdanalysis/issues/2423>`_.
    """
    elems = ['S', 'O', '', 'C', 'Na']
    expectation = ['S', 'O', '', 'C', 'NA']
    dummy_universe_with_elements = dummy_universe_without_elements
    dummy_universe_with_elements.add_TopologyAttr('elements', elems)
    destination = StringIO()
    with mda.coordinates.PDB.PDBWriter(destination) as writer:
        writer.write(dummy_universe_without_elements.atoms)
        content = destination.getvalue()
    element_symbols = [
        line[76:78].strip()
        for line in content.splitlines()
        if line[:6] == 'ATOM  '
    ]
    assert element_symbols == expectation


def test_elements_roundtrip(tmpdir):
    """
    Roundtrip test for PDB elements reading/writing.
    """
    u = mda.Universe(CONECT)
    elements = u.atoms.elements

    outfile = os.path.join(str(tmpdir), 'elements.pdb')
    with mda.coordinates.PDB.PDBWriter(outfile) as writer:
        writer.write(u.atoms)

    u_written = mda.Universe(outfile)

    assert_equal(elements, u_written.atoms.elements)


def test_cryst_meaningless_warning():
    # issue 2599
    # FIXME: This message might change with Issue #2698
    with pytest.warns(UserWarning, match="Unit cell dimensions will be set to zeros."):
        mda.Universe(PDB_CRYOEM_BOX)


def test_cryst_meaningless_select():
    # issue 2599
    u = mda.Universe(PDB_CRYOEM_BOX)
    cur_sele = u.select_atoms('around 0.1 (resid 4 and name CA and segid A)')
    assert cur_sele.n_atoms == 0
