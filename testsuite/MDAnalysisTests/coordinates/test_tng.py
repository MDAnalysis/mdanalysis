# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
from __future__ import division, absolute_import
from six.moves import zip

import MDAnalysis as mda
import numpy as np
import os
import warnings

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_array_almost_equal, dec,
                           assert_almost_equal, )


from MDAnalysisTests.datafiles import (GRO, TNG, PDB)

# from MDAnalysisTests.datafiles import (COORDINATES_XTC, COORDINATES_TOPOLOGY,
#                                        COORDINATES_TRR)
# from MDAnalysisTests.coordinates.base import (MultiframeReaderTest, BaseReference,
#                                               BaseWriterTest,
#                                               assert_timestep_almost_equal)
from MDAnalysisTests import tempdir

# from MDAnalysis.coordinates.TNG import TNGReader

# I want to catch all warnings in the tests. If this is not set at the start it
# could cause test that check for warnings to fail.
warnings.simplefilter('always')


class _GromacsReader(object):
    # This base class assumes same lengths and dt for XTC and TRR test cases!
    filename = None
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.],
                            dtype=np.float32)
    # computed with Gromacs: 362.26999999999998 nm**3 * 1000 A**3/nm**3
    ref_volume = 362270.0

    def setUp(self):
        # loading from GRO is 4x faster than the PDB reader
        self.universe = mda.Universe(GRO, self.filename, convert_units=True)
        self.trajectory = self.universe.trajectory
        self.prec = 3
        self.ts = self.universe.coord
        # dummy output file
        ext = os.path.splitext(self.filename)[1]
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/xdr-reader-test' + ext

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe
        del self.tmpdir

    @dec.slow
    def test_coordinates(self):
        ca_nm = np.array([[6.043369675, 7.385184479, 1.381425762]],
                         dtype=np.float32)
        # coordinates in the base unit (needed for True)
        ca_Angstrom = ca_nm * 10.0
        U = self.universe
        T = U.trajectory
        T.rewind()
        T.next()
        T.next()
        assert_equal(self.ts.frame, 2, "failed to step to frame 3")
        ca = U.select_atoms('name CA and resid 122')
        # low precision match (2 decimals in A, 3 in nm) because the above are
        # the trr coords
        assert_array_almost_equal(ca.positions, ca_Angstrom, 2,
                                  err_msg="coords of Ca of resid 122 do not "
                                  "match for frame 3")

    @dec.slow
    @attr('issue')
    def test_unitcell(self):
        """Test that xtc/trr unitcell is read correctly (Issue 34)"""
        self.universe.trajectory.rewind()
        uc = self.ts.dimensions
        assert_array_almost_equal(
            uc,
            self.ref_unitcell,
            self.prec,
            err_msg="unit cell dimensions (rhombic dodecahedron)")

    @dec.slow
    def test_volume(self):
        # need to reduce precision for test (nm**3 <--> A**3)
        self.universe.trajectory.rewind()
        vol = self.ts.volume
        assert_array_almost_equal(
            vol,
            self.ref_volume,
            0,
            err_msg="unit cell volume (rhombic dodecahedron)")

    @dec.slow
    def test_dt(self):
        assert_almost_equal(self.universe.trajectory.dt,
                            100.0,
                            4,
                            err_msg="wrong timestep dt")

    @dec.slow
    def test_totaltime(self):
        # test_totaltime(): need to reduce precision because dt is only precise
        # to ~4 decimals and accumulating the inaccuracy leads to even lower
        # precision in the totaltime (consequence of fixing Issue 64)
        assert_almost_equal(self.universe.trajectory.totaltime,
                            900.0,
                            3,
                            err_msg="wrong total length of trajectory")

    @dec.slow
    def test_frame(self):
        self.trajectory[4]  # index is 0-based and frames are 0-based
        assert_equal(self.universe.trajectory.frame, 4, "wrong frame number")

    @dec.slow
    def test_time(self):
        self.trajectory[4]
        assert_almost_equal(self.universe.trajectory.time,
                            400.0,
                            3,
                            err_msg="wrong time of frame")


class TestTNGReader(_GromacsReader):
    filename = TNG

    @dec.slow
    def test_velocities(self):
        # frame 0, v in nm/ps
        # from gmxdump -f MDAnalysisTests/data/adk_oplsaa.trr
        #      v[47675]={-7.86469e-01,  1.57479e+00,  2.79722e-01}
        #      v[47676]={ 2.70593e-08,  1.08052e-06,  6.97028e-07}
        v_native = np.array(
            [
                [-7.86469e-01, 1.57479e+00, 2.79722e-01
                 ], [2.70593e-08, 1.08052e-06, 6.97028e-07]
            ],
            dtype=np.float32)

        # velocities in the MDA base unit A/ps (needed for True)
        v_base = v_native * 10.0
        self.universe.trajectory.rewind()
        assert_equal(self.ts.frame, 0, "failed to read frame 1")

        assert_array_almost_equal(
            self.universe.trajectory.ts._velocities[[47675, 47676]], v_base,
            self.prec, err_msg="ts._velocities for indices 47675,47676 do not "
            "match known values")

        assert_array_almost_equal(
            self.universe.atoms.velocities[[47675, 47676]], v_base,
            self.prec, err_msg="velocities for indices 47675,47676 do not "
            "match known values")

        for index, v_known in zip([47675, 47676], v_base):
            assert_array_almost_equal(
                self.universe.atoms[index].velocity,
                v_known,
                self.prec,
                err_msg="atom[{0:d}].velocity does not match known values".format(
                index))


class _XDRNoConversion(object):
    filename = None

    def setUp(self):
        self.universe = mda.Universe(PDB, self.filename, convert_units=False)
        self.ts = self.universe.trajectory.ts

    def tearDown(self):
        del self.universe
        del self.ts

    @dec.slow
    def test_coordinates(self):
        # note: these are the native coordinates in nm
        ca_nm = np.array([[6.043369675, 7.385184479, 1.381425762]],
                         dtype=np.float32)
        U = self.universe
        T = U.trajectory
        T.rewind()
        T.next()
        T.next()
        assert_equal(self.ts.frame, 2, "failed to step to frame 3")
        ca = U.select_atoms('name CA and resid 122')
        # low precision match because we also look at the trr: only 3 decimals
        # in nm in xtc!
        assert_array_almost_equal(ca.positions, ca_nm, 3,
                                  err_msg="native coords of Ca of resid 122 "
                                  "do not match for frame 3 with "
                                  "convert_units=False")


class TestTNGNoConversion(_XDRNoConversion):
    filename = TNG


# class TRRReference(BaseReference):
#     def __init__(self):
#         super(TRRReference, self).__init__()
#         self.trajectory = COORDINATES_TRR
#         self.topology = COORDINATES_TOPOLOGY
#         self.changing_dimensions = True
#         self.reader = mda.coordinates.TRR.TRRReader
#         self.writer = mda.coordinates.TRR.TRRWriter
#         self.ext = 'trr'
#         self.prec = 3
#         self.first_frame.velocities = self.first_frame.positions / 10
#         self.first_frame.forces = self.first_frame.positions / 100

#         self.second_frame.velocities = self.second_frame.positions / 10
#         self.second_frame.forces = self.second_frame.positions / 100

#         self.last_frame.velocities = self.last_frame.positions / 10
#         self.last_frame.forces = self.last_frame.positions / 100

#         self.jump_to_frame.velocities = self.jump_to_frame.positions / 10
#         self.jump_to_frame.forces = self.jump_to_frame.positions / 100

#     def iter_ts(self, i):
#         ts = self.first_frame.copy()
#         ts.positions = 2**i * self.first_frame.positions
#         ts.velocities = ts.positions / 10
#         ts.forces = ts.positions / 100
#         ts.time = i
#         ts.frame = i
#         return ts


# class TestTRRReader_2(MultiframeReaderTest):
#     def __init__(self, reference=None):
#         if reference is None:
#             reference = TRRReference()
#         super(TestTRRReader_2, self).__init__(reference)
