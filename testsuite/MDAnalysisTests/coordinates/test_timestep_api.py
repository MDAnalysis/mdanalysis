# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""Tests for the API definition of Timestep

These come in 2 parts

_TestTimestep tests the functionality of a Timestep in isolation

_TestTimestepInterface tests the Readers are correctly using Timesteps
"""

import numpy as np
from numpy.testing import (TestCase, assert_raises, assert_equal,
                           assert_array_almost_equal, assert_, dec)
from nose.plugins.attrib import attr
from nose.tools import assert_not_equal
from MDAnalysisTests.plugins.knownfailure import knownfailure
from MDAnalysisTests import module_not_found

import MDAnalysis as mda
from MDAnalysis.lib.mdamath import triclinic_vectors
from MDAnalysis import NoDataError
from MDAnalysisTests.datafiles import (
    PSF, DCD, DCD_empty, PDB_small, XPDB_small, PDB_closed, PDB_multiframe,
    PDB, CRD, XTC, TRR, GRO, DMS, CONECT, PDBQT_input,
    XYZ, XYZ_bz2, XYZ_psf, PRM, TRJ, TRJ_bz2, PRMpbc, TRJpbc_bz2, PRMncdf, NCDF,
    PQR,
    PDB_sub_dry, TRR_sub_sol, PDB_sub_sol, TRZ, TRZ_psf, LAMMPSdata,
    LAMMPSdata_mini,
    LAMMPSdata2,LAMMPSdcd2,
    PSF_TRICLINIC, DCD_TRICLINIC, PSF_NAMD_TRICLINIC, DCD_NAMD_TRICLINIC,
    GMS_ASYMOPT, GMS_SYMOPT, GMS_ASYMSURF, XYZ_mini, PFncdf_Top, PFncdf_Trj,
    INPCRD, XYZ_five, mol2_molecules,
    DLP_CONFIG, DLP_HISTORY
    )

from MDAnalysisTests.coordinates.base import BaseTimestepTest


# Can add in custom tests for a given Timestep here!
class TestBaseTimestep(BaseTimestepTest):
    #    def test_nothing(self):
    #        assert_equal(1, 1)
    pass


# using test generator, don't change to TestCase
class TestTimestepEquality(object):
    def _get_pos(self):
        # Get generic reference positions
        return np.arange(30).reshape(10, 3) * 1.234

    def _check_ts(self, a, b, err_msg):
        assert_(a == b, err_msg)
        assert_(b == a, err_msg)

    def test_check_equal(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = self._get_pos()

        ts2 = mda.coordinates.base.Timestep(10)
        ts2.positions = self._get_pos()

        self._check_ts(ts1, ts2, 'Failed on base timestep')

    def test_other_timestep(self):
        # use a subclass to base.Timestep to check it works
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = self._get_pos()

        # can't do TRR Timestep here as it always has vels and forces
        # so isn't actually equal to a position only timestep
        for otherTS in [mda.coordinates.DCD.Timestep,
                        mda.coordinates.TRJ.Timestep,
                        mda.coordinates.DMS.Timestep,
                        mda.coordinates.GRO.Timestep,
                        mda.coordinates.TRZ.Timestep,
                        mda.coordinates.XTC.Timestep]:
            ts2 = otherTS(10)
            ts2.positions = self._get_pos()
            yield self._check_ts, ts1, ts2, "Failed on {0}".format(otherTS)

    def test_trr_timestep(self):
        ts1 = mda.coordinates.base.Timestep(10, velocities=True, forces=True)
        ts1.positions = self._get_pos()
        ts1.velocities = self._get_pos() + 10.0
        ts1.forces = self._get_pos() + 100.0

        ts2 = mda.coordinates.TRR.Timestep(10)
        ts2.positions = self._get_pos()
        ts2.velocities = self._get_pos() + 10.0
        ts2.forces = self._get_pos() + 100.0

        self._check_ts(ts1, ts2, "Failed on TRR Timestep")

    def test_wrong_class(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = self._get_pos()

        b = tuple([0, 1, 2, 3])

        assert_(ts1 != b)
        assert_(b != ts1)

    def test_wrong_frame(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = self._get_pos()

        ts2 = mda.coordinates.base.Timestep(10)
        ts2.positions = self._get_pos()
        ts2.frame = 987

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_wrong_n_atoms(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = self._get_pos()

        ts3 = mda.coordinates.base.Timestep(20)

        assert_(ts1 != ts3)
        assert_(ts3 != ts1)

    def test_wrong_pos(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = self._get_pos()

        ts2 = mda.coordinates.base.Timestep(10)
        ts2.positions = self._get_pos() + 1.0

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_vels(self):
        ts1 = mda.coordinates.base.Timestep(10, velocities=True)
        ts2 = mda.coordinates.base.Timestep(10, velocities=True)

        ts1.velocities = self._get_pos()
        ts2.velocities = self._get_pos()

        assert_(ts1 == ts2)
        assert_(ts2 == ts1)

    def test_check_mismatched_vels(self):
        ts1 = mda.coordinates.base.Timestep(10, velocities=True)
        ts2 = mda.coordinates.base.Timestep(10, velocities=False)

        ts1.velocities = self._get_pos()

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_wrong_vels(self):
        ts1 = mda.coordinates.base.Timestep(10, velocities=True)
        ts2 = mda.coordinates.base.Timestep(10, velocities=True)

        ts1.velocities = self._get_pos()
        ts2.velocities = self._get_pos() + 1.0

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_forces(self):
        ts1 = mda.coordinates.base.Timestep(10, forces=True)
        ts2 = mda.coordinates.base.Timestep(10, forces=True)

        ts1.forces = self._get_pos()
        ts2.forces = self._get_pos()

        assert_(ts1 == ts2)
        assert_(ts2 == ts1)

    def test_check_mismatched_forces(self):
        ts1 = mda.coordinates.base.Timestep(10, forces=True)
        ts2 = mda.coordinates.base.Timestep(10, forces=False)

        ts1.forces = self._get_pos()

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_wrong_forces(self):
        ts1 = mda.coordinates.base.Timestep(10, forces=True)
        ts2 = mda.coordinates.base.Timestep(10, forces=True)

        ts1.forces = self._get_pos()
        ts2.forces = self._get_pos() + 1.0

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)


# TODO: Merge this into generic Reader tests
class BaseTimestepInterfaceTest(object):
    """Test the Timesteps created by Readers

    This checks that Readers are creating correct Timestep objects,
    ensuring a consistent interface when using Timestep objects.

    Failures here are the Reader's fault

    See Issue #250 for discussion
    """
    def tearDown(self):
        del self.ts
        del self.u

    def test_frame(self):
        assert_equal(self.ts.frame, 0)

    def test_dt(self):
        assert_equal(self.u.trajectory.dt, self.ts.dt)


class TestCRD(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(XYZ_five, INPCRD)
        self.ts = u.trajectory.ts


class TestDCD(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(PSF, DCD)
        self.ts = u.trajectory.ts


class TestDLPConfig(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(DLP_CONFIG, format='CONFIG')
        self.ts = u.trajectory.ts


class TestDLPHistory(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(DLP_HISTORY, format='HISTORY')
        self.ts = u.trajectory.ts


class TestDMS(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(DMS)
        self.ts = u.trajectory.ts


class TestGMS(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(GMS_ASYMOPT, GMS_ASYMOPT,
                                  format='GMS', topology_format='GMS')
        self.ts = u.trajectory.ts


class TestGRO(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(GRO)
        self.ts = u.trajectory.ts


class TestINPCRD(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(XYZ_five, INPCRD)
        self.ts = u.trajectory.ts


class TestLAMMPS(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(LAMMPSdata)
        self.ts = u.trajectory.ts


class TestLAMMPSDCD(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(LAMMPSdata2, LAMMPSdcd2,
                                  format='LAMMPS', topology_format='DATA',
                                  timeunit='fs')
        self.ts = u.trajectory.ts


class TestMOL2(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(mol2_molecules)
        self.ts = u.trajectory.ts


class TestPDB(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(PDB_small)
        self.ts = u.trajectory.ts


class TestPDBQT(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(PDBQT_input)
        self.ts = u.trajectory.ts


class TestPQR(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(PQR)
        self.ts = u.trajectory.ts


class TestTRJ(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(PRM, TRJ)
        self.ts = u.trajectory.ts


class TestNCDF(BaseTimestepInterfaceTest):
    @dec.skipif(module_not_found("netCDF4"),
                "Test skipped because netCDF is not available.")
    def setUp(self):
        u = self.u = mda.Universe(PRMncdf, NCDF)
        self.ts = u.trajectory.ts


class TestTRR(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(GRO, TRR)
        self.ts = u.trajectory.ts


class TestTRZ(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(TRZ_psf, TRZ)
        self.ts = u.trajectory.ts


class TestXTC(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(GRO, XTC)
        self.ts = u.trajectory.ts


class TestXYZ(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(XYZ_mini)
        self.ts = u.trajectory.ts


# TODO: Merge these tests into BaseTimestepTest
#       and move to per-format style
class TestTimestepCopy(object):
    """Test copy and copy_slice methods on Timestep"""
    formats = [
        ('DCD', (PSF, DCD)),
        ('DMS', (DMS,)),
        ('GRO', (GRO,)),
        ('PDB', (PDB_small,)),
        ('TRJ', (PRM, TRJ)),
        ('TRR', (GRO, TRR)),
        ('TRZ', (TRZ_psf, TRZ)),
        ('XTC', (PDB, XTC)),
    ]

    def _check_copy(self, name, ref_ts):
        """Check basic copy"""
        ts2 = ref_ts.copy()

        err_msg = ("Timestep copy failed for format {form}"
                   " on attribute {att}")

        # eq method checks:
        # - frame
        # - n_atoms
        # - positions, vels and forces
        assert_(ref_ts == ts2)

        assert_array_almost_equal(ref_ts.dimensions, ts2.dimensions,
                                  decimal=4)

        # Check things not covered by eq
        for d in ref_ts.data:
            assert_(d in ts2.data)
            if isinstance(ref_ts.data[d], np.ndarray):
                assert_array_almost_equal(
                    ref_ts.data[d], ts2.data[d])
            else:
                assert_(ref_ts.data[d] == ts2.data[d])

    def _check_independent(self, name, ts):
        """Check that copies made are independent"""
        ts2 = ts.copy()

        if ts.has_positions:
            self._check_array(ts.positions, ts2.positions)
        if ts.has_velocities:
            self._check_array(ts.velocities, ts2.velocities)
        if ts.has_forces:
            self._check_array(ts.forces, ts2.forces)
        self._check_array(ts.dimensions, ts2.dimensions)

    def _check_array(self, arr1, arr2):
        """Check modifying one array doesn't change other"""
        ref = arr1.copy()
        arr2 += 1.0
        assert_array_almost_equal(ref, arr1)

    def _check_copy_slice_indices(self, name, ts):
        sl = slice(0, len(ts), 3)
        ts2 = ts.copy_slice(sl)
        self._check_slice(ts, ts2, sl)

    def _check_copy_slice_slice(self, name, ts):
        sl = [0, 1, 3, 5, 6, 7]
        ts2 = ts.copy_slice(sl)
        self._check_slice(ts, ts2, sl)

    def _check_slice(self, ts1, ts2, sl):
        if ts1.has_positions:
            assert_array_almost_equal(ts1.positions[sl], ts2.positions)
        if ts1.has_velocities:
            assert_array_almost_equal(ts1.velocities[sl], ts2.velocities)
        if ts1.has_forces:
            assert_array_almost_equal(ts1.forces[sl], ts2.forces)

    def test_copy(self):
        for fname, args in self.formats:
            u = mda.Universe(*args)
            ts = u.trajectory.ts

            yield self._check_copy, fname, ts
            yield self._check_independent, fname, ts
            yield self._check_copy_slice_indices, fname, ts
            yield self._check_copy_slice_slice, fname, ts
