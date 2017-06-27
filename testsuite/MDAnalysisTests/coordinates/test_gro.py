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
import numpy as np
from MDAnalysis.coordinates.GRO import GROReader, GROWriter
from MDAnalysisTests import make_Universe, tempdir
from MDAnalysisTests.coordinates.base import (
    BaseReference, BaseReaderTest, BaseWriterTest, BaseTimestepTest,
)
from MDAnalysisTests.coordinates.reference import RefAdK
from MDAnalysisTests.datafiles import (
    COORDINATES_GRO,
    COORDINATES_GRO_INCOMPLETE_VELOCITY,
    COORDINATES_GRO_BZ2,
    GRO,
    GRO_large,
)
from nose.plugins.attrib import attr
from numpy.testing import (
    assert_,
    assert_almost_equal,
    assert_array_almost_equal,
    dec,
    assert_equal,
    assert_raises
)

class TestGROReaderOld(RefAdK):
    def setUp(self):
        self.universe = mda.Universe(GRO)
        self.ts = self.universe.trajectory.ts
        # lower prec in gro!! (3 decimals nm -> 2 decimals in Angstroem)
        self.prec = 2

    def tearDown(self):
        del self.universe
        del self.ts

    def test_flag_convert_lengths(self):
        assert_equal(mda.core.flags['convert_lengths'], True,
                     "MDAnalysis.core.flags['convert_lengths'] should be True "
                     "by default")

    def test_load_gro(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from small GRO")
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_coordinates(self):
        A10CA = self.universe.atoms.CA[10]
        assert_almost_equal(A10CA.position,
                            self.ref_coordinates['A10CA'],
                            self.prec,
                            err_msg="wrong coordinates for A10:CA")

    def test_distances(self):
        # NOTe that the prec is only 1 decimal: subtracting two low precision
        #      coordinates low prec: 9.3455122920041109; high prec (from pdb):
        #      9.3513174
        NTERM = self.universe.atoms.N[0]
        CTERM = self.universe.atoms.C[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec - 1,
                            err_msg="distance between M1:N and G214:C")

    def test_selection(self):
        na = self.universe.select_atoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size,
                     "Atom selection of last atoms in file")

    def test_unitcell(self):
        assert_array_almost_equal(
            self.ts.dimensions,
            self.ref_unitcell,
            self.prec,
            err_msg="unit cell dimensions (rhombic dodecahedron)")


class TestGROReaderNoConversionOld(RefAdK):
    def setUp(self):
        self.universe = mda.Universe(GRO, convert_units=False)
        self.ts = self.universe.trajectory.ts
        self.prec = 3

    def tearDown(self):
        del self.universe
        del self.ts

    def test_coordinates(self):
        # note: these are the native coordinates in nm; for the test to succeed
        # we loaded with convert_units=False
        A10CA = self.universe.atoms.CA[10]
        # coordinates in nm
        assert_almost_equal(A10CA.position, RefAdK.ref_coordinates['A10CA'] / 10.0,
                            self.prec, err_msg="wrong native coordinates "
                                               "(in nm) for A10:CA")

    def test_distances(self):
        # 3 decimals on nm in gro but we compare to the distance
        # computed from the pdb file, so the effective precision is 2 again.
        # (Otherwise the distance test fails:
        #  Arrays are not almost equal distance between M1:N and G214:C
        #    ACTUAL: 0.93455122920041123
        #    DESIRED: 0.93513173999999988
        NTERM = self.universe.atoms.N[0]
        CTERM = self.universe.atoms.C[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        # coordinates in nm
        assert_almost_equal(d, RefAdK.ref_distances['endtoend'] / 10.0,
                            self.prec - 1, err_msg="distance between M1:N "
                                                   "and G214:C")

    def test_unitcell(self):
        # lengths in A : convert to nm
        assert_array_almost_equal(
            self.ts.dimensions[:3],
            self.ref_unitcell[:3] / 10.0,
            self.prec,
            err_msg="unit cell A,B,C (rhombic dodecahedron)")
        # angles should not have changed
        assert_array_almost_equal(
            self.ts.dimensions[3:],
            self.ref_unitcell[3:],
            self.prec,
            err_msg="unit cell alpha,beta,gamma (rhombic dodecahedron)")

    def test_volume(self):
        # ref lengths in A (which was originally converted from nm)
        assert_almost_equal(
            self.ts.volume,
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
    def __init__(self, reference=None):
        if reference is None:
            reference = GROReference()
        super(TestGROReader, self).__init__(reference)

    def test_flag_convert_lengths(self):
        assert_equal(mda.core.flags['convert_lengths'], True,
                     "MDAnalysis.core.flags['convert_lengths'] should be True "
                     "by default")

    def test_time(self):
        u = mda.Universe(self.ref.topology, self.ref.trajectory)
        assert_equal(u.trajectory.time, 0.0,
                     "wrong time of the frame")

    def test_full_slice(self):
        u = mda.Universe(self.ref.topology, self.ref.trajectory)
        trj_iter = u.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(u.trajectory.n_frames))


class TestGROWriter(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROReference()
        super(TestGROWriter, self).__init__(reference)
        self.u_no_resnames = make_Universe(['names', 'resids'],
                                           trajectory=True)
        self.u_no_resids = make_Universe(['names', 'resnames'],
                                         trajectory=True)
        self.u_no_names = make_Universe(['resids', 'resnames'],
                                        trajectory=True)

    def test_write_velocities(self):
        u = mda.Universe(self.ref.topology, self.ref.trajectory)
        outfile = self.tmp_file('write-velocities-test')
        u.atoms.write(outfile)

        u2 = mda.Universe(outfile)
        assert_array_almost_equal(u.atoms.velocities,
                                  u2.atoms.velocities)

    @dec.slow
    def test_write_no_resnames(self):
        outfile = self.tmp_file('write-no-resnames-test')
        self.u_no_resnames.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.array(['UNK'] * self.u_no_resnames.atoms.n_atoms)
        assert_equal(u.atoms.resnames, expected)

    @dec.slow
    def test_write_no_resids(self):
        outfile = self.tmp_file('write-no-resids-test')
        self.u_no_resids.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.ones((25,))
        assert_equal(u.residues.resids, expected)

    @dec.slow
    def test_writer_no_atom_names(self):
        outfile = self.tmp_file('write-no-names-test')
        self.u_no_names.atoms.write(outfile)
        u = mda.Universe(outfile)
        expected = np.array(['X'] * self.u_no_names.atoms.n_atoms)
        assert_equal(u.atoms.names, expected)

    @dec.slow
    @attr('issue')
    def test_check_coordinate_limits_min(self):
        """Test that illegal GRO coordinates (x <= -999.9995 nm) are caught
        with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(GRO)
        u.atoms[2000].position = [11.589, -999.9995 * 10, 22.2]  # nm -> A
        outfile = self.tmp_file('coordinate-limits-min-test')
        assert_raises(ValueError, u.atoms.write, outfile)

    @dec.slow
    @attr('issue')
    def test_check_coordinate_limits_max(self):
        """Test that illegal GRO coordinates (x > 9999.9995 nm) are caught
        with ValueError (Issue 57)"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(GRO)
        # nm -> A  ; [ob] 9999.9996 not caught
        u.atoms[1000].position = [0, 9999.9999 * 10, 1]
        outfile = self.tmp_file('coordinate-limits-max-test')
        assert_raises(ValueError, u.atoms.write, outfile)

    @dec.slow
    def test_check_coordinate_limits_max_noconversion(self):
        """Test that illegal GRO coordinates (x > 9999.9995 nm) also
        raises exception for convert_units=False"""
        # modify coordinates so we need our own copy or we could mess up
        # parallel tests
        u = mda.Universe(GRO, convert_units=False)
        u.atoms[1000].position = [22.2, 9999.9999, 37.89]
        outfile = self.tmp_file('coordinate-limits-max-noconversion-test')
        assert_raises(ValueError, u.atoms.write, outfile, convert_units=False)


class GRONoConversionReference(GROReference):
    def __init__(self):
        super(GRONoConversionReference, self).__init__()
        self.first_frame.positions /= 10.0
        self.first_frame.velocities /= 10.0
        self.dimensions[:3] /= 10.0
        self.volume /= 1000


class TestGROReaderNoConversion(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GRONoConversionReference()
        super(TestGROReaderNoConversion, self).__init__(reference)
        self.reader = self.ref.reader(self.ref.trajectory, convert_units=False)
        self.reader.add_auxiliary('lowf', self.ref.aux_lowf,
                                  dt=self.ref.aux_lowf_dt,
                                  initial_time=0, time_selector=None)
        self.reader.add_auxiliary('highf', self.ref.aux_highf,
                                  dt=self.ref.aux_highf_dt,
                                  initial_time=0, time_selector=None)


class TestGROWriterNoConversion(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GRONoConversionReference()
        super(TestGROWriterNoConversion, self).__init__(reference)
        self.writer = self.ref.writer(self.ref.trajectory, convert_units=False)


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
    def __init__(self, reference=None):
        if reference is None:
            reference = GROReaderIncompleteVelocitiesReference()
        super(TestGROReaderIncompleteVelocities, self).__init__(reference)


class TestGROWriterIncompleteVelocities(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROReaderIncompleteVelocitiesReference()
        super(TestGROWriterIncompleteVelocities, self).__init__(reference)


class GROBZReference(GROReference):
    def __init__(self):
        super(GROBZReference, self).__init__()
        self.trajectory = COORDINATES_GRO_BZ2
        self.topology = COORDINATES_GRO_BZ2


class TestGROBZ2Reader(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROBZReference()
        super(TestGROBZ2Reader, self).__init__(reference)


class TestGROBZ2Writer(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROBZReference()
        super(TestGROBZ2Writer, self).__init__(reference)


class GROLargeReference(GROReference):
    def __init__(self):
        super(GROLargeReference, self).__init__()
        self.trajectory = GRO_large
        self.topology = GRO_large


class TestGROLargeWriter(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROLargeReference()
        super(TestGROLargeWriter, self).__init__(reference)

    @dec.slow
    @attr('issue')
    def test_writer_large(self):
        """
        Test that atom numbers are truncated for large
        GRO files (Issue 550).
        """
        outfile = self.tmp_file('outfile1.gro')
        u = mda.Universe(self.ref.topology, self.ref.trajectory)
        u.atoms.write(outfile)

        with open(outfile, 'rt') as mda_output:
            with mda.lib.util.anyopen(self.ref.topology, 'rt') as expected_output:
                produced_lines = mda_output.readlines()[1:]
                expected_lines = expected_output.readlines()[1:]
                assert_equal(produced_lines,
                             expected_lines,
                             err_msg="Writing GRO file with > 100 000 "
                                     "coords does not truncate properly.")

    @dec.slow
    @attr('issue')
    def test_writer_large_residue_count(self):
        """
        Ensure large residue number truncation for
        GRO files (Issue 886).
        """
        outfile = self.tmp_file('outfile2.gro')
        u = mda.Universe(self.ref.topology, self.ref.trajectory)
        target_resname = u.residues[-1].resname
        resid_value = 9999999
        u.residues[-1].resid = resid_value
        u.atoms.write(outfile)

        with open(outfile, 'rt') as mda_output:
            output_lines = mda_output.readlines()
            produced_resid = output_lines[-2].split(target_resname)[0]
            expected_resid = str(resid_value)[:5]
            assert_equal(produced_resid,
                         expected_resid,
                         err_msg="Writing GRO file with > 99 999 "
                                 "resids does not truncate properly.")

@tempdir.run_in_tempdir()
def test_growriter_resid_truncation():
    u = make_Universe(extras=['resids'], trajectory=True)
    u.residues[0].resid = 123456789
    u.atoms.write('out.gro')

    with open('out.gro', 'r') as grofile:
        grofile.readline()
        grofile.readline()
        line = grofile.readline()
    # larger digits should get truncated
    assert_(line.startswith('56789UNK'))

class TestGROTimestep(BaseTimestepTest):
    Timestep = mda.coordinates.GRO.Timestep
    name = "GRO"
    has_box = True
    set_box = True
    unitcell = np.array([10., 11., 12.,
                         0., 0., 0.,
                         0., 0., 0.])
    uni_args = (GRO,)

    def test_unitcell_set2(self):
        box = np.array([80.017, 80.017, 80.017, 60.00, 60.00, 90.00],
                       dtype=np.float32)

        ref = np.array([80.00515747, 80.00515747, 56.57218552,  # v1x, v2y, v3z
                        0., 0.,  # v1y v1z
                        0., 0.,  # v2x v2y
                        40.00257874, 40.00257874], dtype=np.float32)  # v3x, v3y
        self.ts.dimensions = box
        assert_array_almost_equal(self.ts._unitcell, ref, decimal=2)
