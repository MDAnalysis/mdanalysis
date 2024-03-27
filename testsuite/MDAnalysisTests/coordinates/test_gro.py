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
from io import StringIO
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.coordinates.GRO import GROReader, GROWriter
from MDAnalysis.transformations import translate
from MDAnalysisTests import make_Universe
from MDAnalysisTests.coordinates.base import (
    BaseReference, BaseReaderTest, BaseWriterTest,
)
from MDAnalysisTests.coordinates.reference import RefAdK
from MDAnalysisTests.datafiles import (
    COORDINATES_GRO,
    COORDINATES_GRO_INCOMPLETE_VELOCITY,
    COORDINATES_GRO_BZ2,
    GRO,
    GRO_large,
    two_water_gro_multiframe,
    GRO_huge_box,
    PDB_closed,
)
from numpy.testing import (
    assert_almost_equal,
    assert_equal,
)
import pytest


class TestGROReaderOld(RefAdK):
    # lower prec in gro!! (3 decimals nm -> 2 decimals in Angstroem)
    prec = 2

    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(GRO)

    def test_load_gro(self, universe):
        U = universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from small GRO")
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_coordinates(self, universe):
        A10CA = universe.select_atoms('name CA')[10]
        assert_almost_equal(A10CA.position,
                            self.ref_coordinates['A10CA'],
                            self.prec,
                            err_msg="wrong coordinates for A10:CA")

    def test_distances(self, universe):
        # NOTe that the prec is only 1 decimal: subtracting two low precision
        #      coordinates low prec: 9.3455122920041109; high prec (from pdb):
        #      9.3513174
        NTERM = universe.select_atoms('name N')[0]
        CTERM = universe.select_atoms('name C')[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec - 1,
                            err_msg="distance between M1:N and G214:C")

    def test_selection(self, universe):
        na = universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size,
                     "Atom selection of last atoms in file")

    def test_unitcell(self, universe):
        assert_almost_equal(
            universe.trajectory.ts.dimensions,
            self.ref_unitcell,
            self.prec,
            err_msg="unit cell dimensions (rhombic dodecahedron)")


class TestGROReaderNoConversionOld(RefAdK):
    prec = 3

    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(GRO, convert_units=False)

    def test_coordinates(self, universe):
        # note: these are the native coordinates in nm; for the test to succeed
        # we loaded with convert_units=False
        A10CA = universe.select_atoms('name CA')[10]
        # coordinates in nm
        assert_almost_equal(A10CA.position,
                            RefAdK.ref_coordinates['A10CA'] / 10.0,
                            self.prec, err_msg="wrong native coordinates "
                                               "(in nm) for A10:CA")

    def test_distances(self, universe):
        # 3 decimals on nm in gro but we compare to the distance
        # computed from the pdb file, so the effective precision is 2 again.
        # (Otherwise the distance test fails:
        #  Arrays are not almost equal distance between M1:N and G214:C
        #    ACTUAL: 0.93455122920041123
        #    DESIRED: 0.93513173999999988
        NTERM = universe.select_atoms('name N')[0]
        CTERM = universe.select_atoms('name C')[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        # coordinates in nm
        assert_almost_equal(d, RefAdK.ref_distances['endtoend'] / 10.0,
                            self.prec - 1, err_msg="distance between M1:N "
                                                   "and G214:C")

    def test_unitcell(self, universe):
        # lengths in A : convert to nm
        assert_almost_equal(
            universe.trajectory.ts.dimensions[:3],
            self.ref_unitcell[:3] / 10.0,
            self.prec,
            err_msg="unit cell A,B,C (rhombic dodecahedron)")
        # angles should not have changed
        assert_almost_equal(
            universe.trajectory.ts.dimensions[3:],
            self.ref_unitcell[3:],
            self.prec,
            err_msg="unit cell alpha,beta,gamma (rhombic dodecahedron)")

    def test_volume(self, universe):
        # ref lengths in A (which was originally converted from nm)
        assert_almost_equal(
            universe.trajectory.ts.volume,
            self.ref_volume / 1000.,
            3,
            err_msg="wrong volume for unitcell (rhombic dodecahedron)")


class GROReference(BaseReference):
    def __init__(self):
        super(GROReference, self).__init__()
        self.trajectory = COORDINATES_GRO
        self.topology = COORDINATES_GRO
        self.reader = GROReader
        self.writer = GROWriter
        self.ext = 'gro'
        self.n_frames = 1
        self.prec = 4
        self.first_frame.velocities = np.array(
            [[0.0000, 0.100, 0.200],
             [0.300, 0.400, 0.500],
             [0.600, 0.700, 0.800],
             [0.900, 1.000, 1.100],
             [1.200, 1.300, 1.400]],
            dtype=np.float32)
        self.totaltime = 0
        self.container_format = True


class TestGROReader(BaseReaderTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return GROReference()

    def test_time(self, ref, reader):
        u = mda.Universe(ref.topology, ref.trajectory)
        assert_equal(u.trajectory.time, 0.0,
                     "wrong time of the frame")

    def test_full_slice(self, ref, reader):
        u = mda.Universe(ref.topology, ref.trajectory)
        trj_iter = u.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(u.trajectory.n_frames))


class TestGROWriter(BaseWriterTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return GROReference()

    def test_write_velocities(self, ref, tmpdir):
        u = mda.Universe(ref.topology, ref.trajectory)
        with tmpdir.as_cwd():
            outfile = 'write-velocities-test.' + ref.ext
            u.atoms.write(outfile)

            u2 = mda.Universe(outfile)
            assert_almost_equal(u.atoms.velocities,
                                u2.atoms.velocities)

    def test_write_no_resnames(self, u_no_resnames, ref, tmpdir):
        outfile = 'write-no-resnames-test.' + ref.ext
        with tmpdir.as_cwd():
            u_no_resnames.atoms.write(outfile)
            u = mda.Universe(outfile)
            expected = np.array(['UNK'] * u_no_resnames.atoms.n_atoms)
            assert_equal(u.atoms.resnames, expected)

    def test_write_no_resids(self, u_no_resids, ref, tmpdir):
        outfile = 'write-no-resids-test.' + ref.ext
        with tmpdir.as_cwd():
            u_no_resids.atoms.write(outfile)
            u = mda.Universe(outfile)
            expected = np.ones((25,))
            assert_equal(u.residues.resids, expected)

    def test_writer_no_atom_names(self, u_no_names, ref, tmpdir):
        outfile = 'write-no-names-test.' + ref.ext
        with tmpdir.as_cwd():
            u_no_names.atoms.write(outfile)
            u = mda.Universe(outfile)
            expected = np.array(['X'] * u_no_names.atoms.n_atoms)
            assert_equal(u.atoms.names, expected)

    def test_check_coordinate_limits_min(self, ref, tmpdir):
        """Test that illegal GRO coordinates (x <= -999.9995 nm) are caught
        with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(GRO)
        u.atoms[2000].position = [11.589, -999.9995 * 10, 22.2]  # nm -> A
        outfile = 'coordinate-limits-min-test.' + ref.ext
        with tmpdir.as_cwd():
            with pytest.raises(ValueError):
                u.atoms.write(outfile)

    def test_check_coordinate_limits_max(self, ref, tmpdir):
        """Test that illegal GRO coordinates (x > 9999.9995 nm) are caught
        with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(GRO)
        # nm -> A  ; [ob] 9999.9996 not caught
        u.atoms[1000].position = [0, 9999.9999 * 10, 1]
        outfile = 'coordinate-limits-max-test.' + ref.ext
        with tmpdir.as_cwd():
            with pytest.raises(ValueError):
                u.atoms.write(outfile)

    def test_check_coordinate_limits_max_noconversion(self, ref, tmpdir):
        """Test that illegal GRO coordinates (x > 9999.9995 nm) also
        raises exception for convert_units=False"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(GRO, convert_units=False)
        u.atoms[1000].position = [22.2, 9999.9999, 37.89]
        outfile = 'coordinate-limits-max-noconversion-test.' + ref.ext
        with tmpdir.as_cwd():
            with pytest.raises(ValueError):
                u.atoms.write(outfile, convert_units=False)


class GRONoConversionReference(GROReference):
    def __init__(self):
        super(GRONoConversionReference, self).__init__()
        self.first_frame.positions /= 10.0
        self.first_frame.velocities /= 10.0
        self.dimensions[:3] /= 10.0
        self.volume /= 1000


class TestGROReaderNoConversion(BaseReaderTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return GRONoConversionReference()

    @staticmethod
    @pytest.fixture(scope='class')
    def reader(ref):
        reader = ref.reader(ref.trajectory, convert_units=False)
        reader.add_auxiliary('lowf', ref.aux_lowf,
                             dt=ref.aux_lowf_dt,
                             initial_time=0, time_selector=None)
        reader.add_auxiliary('highf', ref.aux_highf,
                             dt=ref.aux_highf_dt,
                             initial_time=0, time_selector=None)
        return reader
    
    @staticmethod
    @pytest.fixture(scope='class')
    def transformed(ref):
        transformed = ref.reader(ref.trajectory, convert_units=False)
        transformed.add_transformations(translate([1,1,1]), translate([0,0,0.33]))
        return transformed


class TestGROWriterNoConversion(BaseWriterTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return GRONoConversionReference()

    @staticmethod
    @pytest.fixture(scope='class')
    def writer(ref):
        return ref.writer(ref.trajectory, convert_units=False)


class GROReaderIncompleteVelocitiesReference(GROReference):
    def __init__(self):
        super(GROReaderIncompleteVelocitiesReference, self).__init__()
        self.trajectory = COORDINATES_GRO_INCOMPLETE_VELOCITY
        self.topology = COORDINATES_GRO_INCOMPLETE_VELOCITY
        self.first_frame.velocities = np.array(
            [[0.0000, 0.100, 0.200],
             [0.000, 0.000, 0.000],
             [0.600, 0.700, 0.800],
             [0.900, 1.000, 1.100],
             [1.200, 1.300, 1.400]],
            dtype=np.float32)


class TestGROReaderIncompleteVelocities(BaseReaderTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return GROReaderIncompleteVelocitiesReference()


class TestGROWriterIncompleteVelocities(BaseWriterTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return GROReaderIncompleteVelocitiesReference()


class GROBZReference(GROReference):
    def __init__(self):
        super(GROBZReference, self).__init__()
        self.trajectory = COORDINATES_GRO_BZ2
        self.topology = COORDINATES_GRO_BZ2


class TestGROBZ2Reader(BaseReaderTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return GROBZReference()


class TestGROBZ2Writer(BaseWriterTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return GROBZReference()


class GROLargeReference(GROReference):
    def __init__(self):
        super(GROLargeReference, self).__init__()
        self.trajectory = GRO_large
        self.topology = GRO_large


class TestGROLargeWriter(BaseWriterTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return GROLargeReference()

    def test_writer_large(self, ref, tmpdir):
        """
        Test that atom numbers are truncated for large
        GRO files (Issue 550).
        """
        outfile = 'outfile1.' + ref.ext
        u = mda.Universe(ref.topology, ref.trajectory)
        with tmpdir.as_cwd():
            u.atoms.write(outfile)

            with open(outfile, 'rt') as mda_output:
                with mda.lib.util.anyopen(ref.topology, 'rt') as expected_output:
                    produced_lines = mda_output.readlines()[1:]
                    expected_lines = expected_output.readlines()[1:]
                    assert_equal(produced_lines,
                                 expected_lines,
                                 err_msg="Writing GRO file with > 100 000 "
                                         "coords does not truncate properly.")

    def test_writer_large_residue_count(self, ref, tmpdir):
        """
        Ensure large residue number truncation for
        GRO files (Issue 886).
        """
        outfile = 'outfile2.' + ref.ext
        u = mda.Universe(ref.topology, ref.trajectory)
        target_resname = u.residues[-1].resname
        resid_value = 9999999
        u.residues[-1].resid = resid_value
        with tmpdir.as_cwd():
            u.atoms.write(outfile)

            with open(outfile, 'rt') as mda_output:
                output_lines = mda_output.readlines()
                produced_resid = output_lines[-2].split(target_resname)[0]
                expected_resid = str(resid_value)[:5]
                assert_equal(produced_resid,
                             expected_resid,
                             err_msg="Writing GRO file with > 99 999 "
                                     "resids does not truncate properly.")


def test_growriter_resid_truncation(tmpdir):
    with tmpdir.as_cwd():
        u = make_Universe(extras=['resids'], trajectory=True)
        u.residues[0].resid = 123456789
        u.atoms.write('out.gro')

        with open('out.gro', 'r') as grofile:
            grofile.readline()
            grofile.readline()
            line = grofile.readline()
        # larger digits should get truncated
        assert line.startswith('56789UNK')

class TestGrowriterReindex(object):
    @pytest.fixture()
    def u(self):
        gro = '''test
1
    2CL      CL20850   0.000   0.000   0.000
7.29748 7.66094 9.82962'''
        u = mda.Universe(StringIO(gro), format='gro')
        u.atoms[0].id = 3
        return u

    def test_growriter_resid_true(self, u, tmpdir):
        with tmpdir.as_cwd():
            u.atoms.write('temp.gro', reindex=True)

            with open('temp.gro', 'r') as grofile:
                grofile.readline()
                grofile.readline()
                line = grofile.readline()
            assert line.startswith('    2CL      CL    1')

    def test_growriter_resid_false(self, u, tmpdir):
        with tmpdir.as_cwd():
            u.atoms.write('temp.gro', reindex=False)
            with open('temp.gro', 'r') as grofile:
                grofile.readline()
                grofile.readline()
                line = grofile.readline()
            assert line.startswith('    2CL      CL    3')

    def test_writer_resid_false(self, u, tmpdir):
        with tmpdir.as_cwd():
            with mda.Writer('temp.gro', reindex=False) as w:
                w.write(u.atoms)
            with open('temp.gro', 'r') as grofile:
                grofile.readline()
                grofile.readline()
                line = grofile.readline()
            assert line.startswith('    2CL      CL    3')

    def test_writer_resid_true(self, u, tmpdir):
        with tmpdir.as_cwd():
            with mda.Writer('temp.gro', reindex=True) as w:
                w.write(u.atoms)
            with open('temp.gro', 'r') as grofile:
                grofile.readline()
                grofile.readline()
                line = grofile.readline()
            assert line.startswith('    2CL      CL    1')


def test_multiframe_gro():
    u = mda.Universe(two_water_gro_multiframe)

    # for now, single frame read
    assert len(u.trajectory) == 1
    assert_equal(u.dimensions, np.array([100, 100, 100, 90, 90, 90], dtype=np.float32))


def test_huge_box_gro():
    u = mda.Universe(GRO_huge_box)

    assert_equal(u.dimensions, np.array([4.e+05, 4.e+05, 4.e+05, 90, 90, 90],
                                        dtype=np.float32))


gro_no_dims = """\
Single Atom no dims GRO
1
    1SOL     OW    1   1.000   1.000   1.000
"""


@pytest.mark.parametrize('dims', [1, 2, 4, 5, 6, 7, 8])
def test_bad_box(dims):
    cell = '   '.join([str(float(i)) for i in range(dims)])
    grofile = gro_no_dims + cell

    errmsg = "GRO unitcell has neither 3 nor 9 entries."
    with pytest.raises(ValueError, match=errmsg):
        u = mda.Universe(StringIO(grofile), format='gro')


def test_gro_empty_box_write_read(tmpdir):
    # Issue #3305 - ensure that read/write deals with None dimensions the same
    u = mda.Universe(PDB_closed)

    # Check lack of dimensions
    assert u.dimensions is None

    with tmpdir.as_cwd():
        wmsg = " setting unit cell to zeroed box"
        with pytest.warns(UserWarning, match=wmsg):
            u.atoms.write('test.gro')

        wmsg = "treating as missing unit cell"
        with pytest.warns(UserWarning, match=wmsg):
            u2 = mda.Universe('test.gro')
        assert u2.dimensions is None

def test_gro_pathlib_singleframereaderbase():
    top = Path(GRO)
    assert isinstance(top, Path)
    u = mda.Universe(top)
    assert u.atoms.n_atoms == 47681

def test_gro_string_singleframereaderbase():
    assert isinstance(GRO, str)
    u = mda.Universe(GRO)
    assert u.atoms.n_atoms == 47681
