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

import itertools
import numpy as np
from numpy.testing import (TestCase, assert_raises, assert_equal, assert_allclose,
                           assert_array_almost_equal, assert_)
from nose.plugins.attrib import attr
from nose.tools import assert_not_equal
from MDAnalysisTests.plugins.knownfailure import knownfailure

import MDAnalysis as mda
from MDAnalysis.lib.mdamath import triclinic_vectors
from MDAnalysis import NoDataError
from MDAnalysisTests.datafiles import (
    PSF, DCD, DCD_empty, PDB_small, XPDB_small, PDB_closed, PDB_multiframe,
    PDB, CRD, XTC, TRR, GRO, DMS, CONECT, PDBQT_input,
    XYZ, XYZ_bz2, XYZ_psf, PRM, TRJ, TRJ_bz2, PRMpbc, TRJpbc_bz2, PRMncdf, NCDF, PQR,
    PDB_sub_dry, TRR_sub_sol, PDB_sub_sol, TRZ, TRZ_psf, LAMMPSdata, LAMMPSdata_mini,
    LAMMPSdata2,LAMMPSdcd2,
    PSF_TRICLINIC, DCD_TRICLINIC, PSF_NAMD_TRICLINIC, DCD_NAMD_TRICLINIC,
    GMS_ASYMOPT, GMS_SYMOPT, GMS_ASYMSURF, XYZ_mini, PFncdf_Top, PFncdf_Trj,
    INPCRD, XYZ_five, mol2_molecules,
    DLP_CONFIG, DLP_HISTORY
    )


# Subclass this and change values where necessary for each format's Timestep.
class _BaseTimestep(object):
    Timestep = mda.coordinates.base.Timestep  # define the class made in test
    name = "base"  # for error messages only
    size = 10  # size of arrays, 10 is enough to allow slicing etc
    refpos = np.arange(size * 3, dtype=np.float32).reshape(size, 3)  # each coord is unique
    has_box = True
    set_box = True  # whether you can set dimensions info.
    # If you can set box, what the underlying unitcell should be if dimensions are:
    newbox = np.array([10., 11., 12., 90., 90., 90.])
    unitcell = np.array([10., 11., 12., 90., 90., 90.])
    ref_volume = 1320.  # what the volume is after setting newbox


class _TRZTimestep(_BaseTimestep):
    Timestep = mda.coordinates.TRZ.Timestep
    name = "TRZ"
    has_box = True
    set_box = True
    unitcell = np.array([10., 0., 0.,
                         0., 11., 0.,
                         0., 0., 12.])


class _DCDTimestep(_BaseTimestep):
    Timestep = mda.coordinates.DCD.Timestep
    name = "DCD"
    has_box = True
    set_box = True
    unitcell = np.array([10., 90., 11., 90., 90., 12.])
    # Could test setting TS box order and checking that this works


class _DMSTimestep(_BaseTimestep):
    Timestep = mda.coordinates.DMS.Timestep
    name = "DMS"
    has_box = True
    unitcell = {'x':np.array([10., 0, 0]),
                'y':np.array([0, 11., 0]),
                'z':np.array([0, 0, 12.])}


class _GROTimestep(_BaseTimestep):
    Timestep = mda.coordinates.GRO.Timestep
    name = "GRO"
    has_box = True
    set_box = True
    unitcell = np.array([10., 11., 12.,
                         0., 0., 0.,
                         0., 0., 0.])


class _TRJTimestep(_BaseTimestep):
    Timestep = mda.coordinates.TRJ.Timestep
    name = "TRJ"
    has_box = True
    set_box = True
    unitcell = np.array([10., 11., 12., 90., 90., 90.])


class _XTCTimestep(_BaseTimestep):
    Timestep = mda.coordinates.xdrfile.core.Timestep
    name = "XTC"
    has_box = True
    set_box = True
    unitcell = np.array([[10., 0., 0.],
                         [0., 11., 0.],
                         [0., 0., 12.]])


class _TRRTimestep(_XTCTimestep):
    Timestep = mda.coordinates.TRR.Timestep
    name = "TRR"


class _DLPolyTimestep(_XTCTimestep):
    Timestep = mda.coordinates.DLPoly.Timestep
    name = "DLPoly"


class _TestTimestep(TestCase):
    """Test all the base functionality of a Timestep

    All Timesteps must pass these tests!

    These test the Timestep independent of the Reader which it
    comes into contact with.  Failures here are the Timesteps fault.
    """

    def setUp(self):
        self.ts = self.Timestep(self.size)
        self.ts.frame += 1
        self.ts.positions = self.refpos

    def tearDown(self):
        del self.ts

    def test_getitem(self):
        assert_equal(self.ts[1], self.refpos[1])

    def test_getitem_neg(self):
        assert_equal(self.ts[-1], self.refpos[-1])

    def test_getitem_neg_IE(self):
        assert_raises(IndexError, self.ts.__getitem__, -(self.size + 1))

    def test_getitem_pos_IE(self):
        assert_raises(IndexError, self.ts.__getitem__, (self.size + 1))

    def test_getitem_slice(self):
        assert_equal(len(self.ts[:2]), len(self.refpos[:2]))
        assert_allclose(self.ts[:2], self.refpos[:2])

    def test_getitem_slice2(self):
        assert_equal(len(self.ts[1::2]), len(self.refpos[1::2]))
        assert_allclose(self.ts[1::2], self.refpos[1::2])

    def test_getitem_ndarray(self):
        sel = np.array([0, 1, 4])
        assert_equal(len(self.ts[sel]), len(self.refpos[sel]))
        assert_allclose(self.ts[sel], self.refpos[sel])

    def test_getitem_TE(self):
        assert_raises(TypeError, self.ts.__getitem__, 'string')

    def test_len(self):
        assert_equal(len(self.ts), self.size)

    def test_iter(self):
        for a, b in itertools.izip(self.ts, self.refpos):
            assert_allclose(a, b)
        assert_equal(len(list(self.ts)), self.size)

    def test_repr(self):
        assert_equal(type(repr(self.ts)), str)

    # Test copy done as separate test

    # Dimensions has 2 possible cases
    # Timestep doesn't do dimensions, should raise NotImplementedError for .dimension and .volume
    # Timestep does do them, should return values properly
    def test_dimensions(self):
        if self.has_box:
            assert_allclose(self.ts.dimensions, np.zeros(6, dtype=np.float32))
        else:
            assert_raises(NotImplementedError, getattr, self.ts, "dimensions")

    def test_dimensions_set_box(self):
        if self.set_box:
            self.ts.dimensions = self.newbox
            assert_allclose(self.ts._unitcell, self.unitcell)
            assert_allclose(self.ts.dimensions, self.newbox)
        else:
            pass

    def test_volume(self):
        if self.has_box and self.set_box:
            self.ts.dimensions = self.newbox
            assert_equal(self.ts.volume, self.ref_volume)
        elif self.has_box and not self.set_box:
            pass  # How to test volume of box when I don't set unitcell first?
        else:
            assert_raises(NotImplementedError, getattr, self.ts, "volume")

    def test_triclinic_vectors(self):
        assert_allclose(self.ts.triclinic_dimensions, triclinic_vectors(self.ts.dimensions))

    def test_set_triclinic_vectors(self):
        ref_vec = triclinic_vectors(self.newbox)
        self.ts.triclinic_dimensions = ref_vec
        assert_equal(self.ts.dimensions, self.newbox)
        assert_allclose(self.ts._unitcell, self.unitcell)

    @attr('issue')
    def test_coordinate_getter_shortcuts(self):
        """testing that reading _x, _y, and _z works as expected (Issue 224) (TestTimestep)"""
        assert_allclose(self.ts._x, self.ts._pos[:,0])
        assert_allclose(self.ts._y, self.ts._pos[:,1])
        assert_allclose(self.ts._z, self.ts._pos[:,2])

    def test_coordinate_setter_shortcuts(self):
        # Check that _x _y and _z are read only
        for coordinate in ('_x', '_y', '_z'):
            random_positions = np.arange(self.size).astype(np.float32)
            assert_raises(AttributeError, setattr, self.ts, coordinate, random_positions)

    # n_atoms should be a read only property
    # all Timesteps require this attribute
    def test_n_atoms(self):
        assert_equal(self.ts.n_atoms, self.ts._n_atoms)

    def test_n_atoms_readonly(self):
        assert_raises(AttributeError, self.ts.__setattr__, 'n_atoms', 20)

    def test_n_atoms_presence(self):
        assert_equal(hasattr(self.ts, '_n_atoms'), True)

    def test_unitcell_presence(self):
        assert_equal(hasattr(self.ts, '_unitcell'), True)

    def test_data_presence(self):
        assert_equal(hasattr(self.ts, 'data'), True)
        assert_equal(isinstance(self.ts.data, dict), True)

    def test_allocate_velocities(self):
        assert_equal(self.ts.has_velocities, False)
        assert_raises(NoDataError, getattr, self.ts, 'velocities')

        self.ts.has_velocities = True
        assert_equal(self.ts.has_velocities, True)
        assert_equal(self.ts.velocities.shape, (self.ts.n_atoms, 3))

    def test_allocate_forces(self):
        assert_equal(self.ts.has_forces, False)
        assert_raises(NoDataError, getattr, self.ts, 'forces')

        self.ts.has_forces = True
        assert_equal(self.ts.has_forces, True)
        assert_equal(self.ts.forces.shape, (self.ts.n_atoms, 3))

    def test_velocities_remove(self):
        ts = self.Timestep(10, velocities=True)
        ts.frame += 1
        assert_equal(ts.has_velocities, True)

        ts.has_velocities = False
        assert_equal(ts.has_velocities, False)
        assert_raises(NoDataError, getattr, ts, 'velocities')

    def test_forces_remove(self):
        ts = self.Timestep(10, forces=True)
        ts.frame += 1
        assert_equal(ts.has_forces, True)

        ts.has_forces = False
        assert_equal(ts.has_forces, False)
        assert_raises(NoDataError, getattr, ts, 'forces')

    def test_from_coordinates_1(self):
        ts = self.Timestep.from_coordinates(self.refpos)

        assert_equal(len(ts), self.size)
        assert_array_almost_equal(ts.positions, self.refpos)
        assert_raises(NoDataError, getattr, ts, 'velocities')
        assert_raises(NoDataError, getattr, ts, 'forces')

    def test_from_coordinates_2(self):
        refvel = self.refpos + 10

        ts = self.Timestep.from_coordinates(self.refpos, velocities=refvel)

        assert_equal(len(ts), self.size)
        assert_array_almost_equal(ts.positions, self.refpos)
        assert_array_almost_equal(ts.velocities, refvel)
        assert_raises(NoDataError, getattr, ts, 'forces')

    def test_from_coordinates_3(self):
        reffor = self.refpos + 10

        ts = self.Timestep.from_coordinates(self.refpos, forces=reffor)

        assert_equal(len(ts), self.size)
        assert_array_almost_equal(ts.positions, self.refpos)
        assert_raises(NoDataError, getattr, ts, 'velocities')
        assert_array_almost_equal(ts.forces, reffor)

    def test_from_coordinates_4(self):
        refvel = self.refpos + 20
        reffor = self.refpos + 10

        ts = self.Timestep.from_coordinates(self.refpos, velocities=refvel,
                                            forces=reffor)

        assert_equal(len(ts), self.size)
        assert_array_almost_equal(ts.positions, self.refpos)
        assert_array_almost_equal(ts.velocities, refvel)
        assert_array_almost_equal(ts.forces, reffor)

    # Time related tests
    def test_supply_dt(self):
        # Check that this gets stored in data properly
        ts = self.Timestep(20, dt=0.04)

        assert_equal(ts.data['dt'], 0.04)
        assert_equal(ts.dt, 0.04)

    def test_redefine_dt(self):
        ts = self.Timestep(20, dt=0.04)
        assert_equal(ts.data['dt'], 0.04)
        assert_equal(ts.dt, 0.04)
        ts.dt = refdt = 0.46
        assert_equal(ts.data['dt'], refdt)
        assert_equal(ts.dt, refdt)

    def test_delete_dt(self):
        ts = self.Timestep(20, dt=0.04)
        assert_equal(ts.data['dt'], 0.04)
        assert_equal(ts.dt, 0.04)
        del ts.dt
        assert_equal('dt' in ts.data, False)
        assert_equal(ts.dt, 1.0)  # default value

    def test_supply_time_offset(self):
        ts = self.Timestep(20, time_offset=100.0)

        assert_equal(ts.data['time_offset'], 100.0)

    def test_time(self):
        ts = self.Timestep(20)
        ts.frame = 0
        ts.time = reftime = 1234.0

        assert_equal(ts.time, reftime)

    def test_delete_time(self):
        ts = self.Timestep(20)
        ts.frame = 0
        ts.time = reftime = 1234.0

        assert_equal(ts.time, reftime)
        del ts.time
        assert_equal(ts.time, 0.0)  # default to 1.0 (dt) * 0 (frame)

    def test_time_with_offset(self):
        ref_offset = 2345.0
        ts = self.Timestep(20, time_offset=ref_offset)
        ts.frame = 0
        ts.time = reftime = 1234.0

        assert_equal(ts.time, reftime + ref_offset)

    def test_dt(self):
        ref_dt = 45.0
        ts = self.Timestep(20, dt=ref_dt)

        for i in range(10):
            ts.frame = i
            assert_equal(ts.time, i * ref_dt)

    def test_change_dt(self):
        ref_dt = 45.0
        ts = self.Timestep(20, dt=ref_dt)

        for i in range(10):
            ts.frame = i
            assert_equal(ts.time, i * ref_dt)

        ts.dt = ref_dt = 77.0

        for i in range(10):
            ts.frame = i
            assert_equal(ts.time, i * ref_dt)

    def test_dt_with_offset(self):
        ref_dt = 45.0
        ref_offset = 2345.0
        ts = self.Timestep(20, dt=ref_dt, time_offset=ref_offset)

        for i in range(10):
            ts.frame = i
            assert_equal(ts.time, i * ref_dt + ref_offset)

    def test_time_overrides_dt_with_offset(self):
        ref_dt = 45.0
        ref_offset = 2345.0
        ts = self.Timestep(20, dt=ref_dt, time_offset=ref_offset)

        ts.frame = 0
        ts.time = reftime = 456.7

        assert_equal(ts.time, reftime + ref_offset)


# Can add in custom tests for a given Timestep here!
class TestBaseTimestep(_TestTimestep, _BaseTimestep):
    #    def test_nothing(self):
    #        assert_equal(1, 1)
    pass


class TestTRZTimestep(_TestTimestep, _TRZTimestep):
    pass


class TestDCDTimestep(_TestTimestep, _DCDTimestep):
    def test_ts_order_define(self):
        """Check that users can hack in a custom unitcell order"""
        old = self.Timestep._ts_order
        self.ts._ts_order = [0, 2, 5, 1, 3, 4]
        self.ts.dimensions = np.array([10, 11, 12, 80, 85, 90])
        assert_allclose(self.ts._unitcell, np.array([10, 80, 11, 85, 90, 12]))
        self.ts._ts_order = old
        self.ts.dimensions = np.zeros(6)


class TestDMSTimestep(_TestTimestep, _DMSTimestep):
    def test_dimensions_set_box(self):
        self.ts.dimensions = self.newbox
        assert_equal(self.ts.dimensions, self.newbox)
        assert_equal(self.ts._unitcell, self.unitcell)

    def test_set_triclinic_vectors(self):
        ref_vec = triclinic_vectors(self.newbox)
        self.ts.triclinic_dimensions = ref_vec
        assert_equal(self.ts.dimensions, self.newbox)
        assert_equal(self.ts._unitcell, self.unitcell)


class TestGROTimestep(_TestTimestep, _GROTimestep):
    def test_unitcell_set2(self):
        old = self.ts._unitcell
        box = np.array([80.017, 80.017, 80.017, 60.00, 60.00, 90.00], dtype=np.float32)

        ref = np.array([80.00515747, 80.00515747, 56.57218552,  # v1x, v2y, v3z
                        0., 0.,  # v1y v1z
                        0., 0.,  # v2x v2y
                        40.00257874, 40.00257874],dtype=np.float32)  # v3x, v3y
        self.ts.dimensions = box
        assert_array_almost_equal(self.ts._unitcell, ref, decimal=2)

        self.ts._unitcell = old


class TestTRJTimestep(_TestTimestep, _TRJTimestep):
    pass


class TestXTCTimestep(_TestTimestep, _XTCTimestep):
    pass


class TestTRRTimestep(_TestTimestep, _TRRTimestep):
    def test_velocities_remove(self):
        # This test is different because TRR requires that the
        # has flags get updated every frame
        ts = self.Timestep(10, velocities=True)
        ts.frame += 1
        ts.has_velocities = True  # This line is extra for TRR
        assert_equal(ts.has_velocities, True)

        ts.has_velocities = False
        assert_equal(ts.has_velocities, False)
        assert_raises(NoDataError, getattr, ts, 'velocities')

    def test_forces_remove(self):
        # This test is different because TRR requires that the
        # has flags get updated every frame
        ts = self.Timestep(10, forces=True)
        ts.frame += 1
        ts.has_forces = True  # This line is extra for TRR
        assert_equal(ts.has_forces, True)

        ts.has_forces = False
        assert_equal(ts.has_forces, False)
        assert_raises(NoDataError, getattr, ts, 'forces')

    def test_positions_expiry(self):
        assert_equal(self.ts.has_positions, True)
        assert_array_almost_equal(self.ts.positions, self.refpos)
        self.ts.frame += 1
        assert_equal(self.ts.has_positions, False)
        assert_raises(NoDataError, getattr, self.ts, 'positions')

    def test_velocities_expiry(self):
        self.ts.velocities = self.refpos + 10
        assert_equal(self.ts.has_velocities, True)
        assert_array_almost_equal(self.ts.velocities, self.refpos + 10)
        self.ts.frame += 1
        assert_equal(self.ts.has_velocities, False)
        assert_raises(NoDataError, getattr, self.ts, 'velocities')

    def test_forces_expiry(self):
        self.ts.forces = self.refpos + 100
        assert_equal(self.ts.has_forces, True)
        assert_array_almost_equal(self.ts.forces, self.refpos + 100)
        self.ts.frame += 1
        assert_equal(self.ts.has_forces, False)
        assert_raises(NoDataError, getattr, self.ts, 'forces')

    def test_positions_renewal(self):
        assert_equal(self.ts.has_positions, True)
        self.ts.frame += 1
        assert_equal(self.ts.has_positions, False)
        assert_raises(NoDataError, getattr, self.ts, 'positions')

        self.ts.has_positions = True
        assert_equal(self.ts.has_positions, True)
        assert_equal(self.ts.positions, self.refpos)

    def test_velocities_renewal(self):
        self.ts.velocities = self.refpos + 10
        assert_equal(self.ts.has_velocities, True)
        assert_array_almost_equal(self.ts.velocities, self.refpos + 10)
        self.ts.frame += 1
        assert_equal(self.ts.has_velocities, False)
        assert_raises(NoDataError, getattr, self.ts, 'velocities')

        self.ts.velocities = self.refpos + 11
        assert_equal(self.ts.has_velocities, True)
        assert_array_almost_equal(self.ts.velocities, self.refpos + 11)

    def test_forces_renewal(self):
        self.ts.forces = self.refpos + 100
        assert_equal(self.ts.has_forces, True)
        self.ts.frame += 1
        assert_equal(self.ts.has_forces, False)
        assert_raises(NoDataError, getattr, self.ts, 'forces')

        self.ts.forces = self.refpos + 101
        assert_equal(self.ts.has_forces, True)
        assert_array_almost_equal(self.ts.forces, self.refpos + 101)


class TestDLPolyTimestep(_TestTimestep, _DLPolyTimestep):
    pass


class TestTimestepEquality(object):  # using test generator, don't change to TestCase
    def test_check_equal(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = np.arange(30).reshape(10, 3)

        ts2 = mda.coordinates.base.Timestep(10)
        ts2.positions = np.arange(30).reshape(10, 3)

        assert_(ts1 == ts2)
        assert_(ts2 == ts1)

    def _check_ts(self, a, b, err_msg):
        assert_(a == b, err_msg)
        assert_(b == a, err_msg)

    def test_other_timestep(self):
        # use a subclass to base.Timestep to check it works
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = np.arange(30).reshape(10, 3)

        # can't do TRR Timestep here as it always has vels and forces
        # so isn't actually equal to a position only timestep
        for otherTS in [mda.coordinates.DCD.Timestep, mda.coordinates.TRJ.Timestep,
                        mda.coordinates.DMS.Timestep, mda.coordinates.GRO.Timestep,
                        mda.coordinates.TRZ.Timestep, mda.coordinates.XTC.Timestep]:
            ts2 = otherTS(10)
            ts2.positions = np.arange(30).reshape(10, 3)
            yield self._check_ts, ts1, ts2, "Failed on {0}".format(otherTS)

    def test_trr_timestep(self):
        ts1 = mda.coordinates.base.Timestep(10, velocities=True, forces=True)
        ts1.positions = np.arange(30).reshape(10, 3)
        ts1.velocities = np.arange(30).reshape(10, 3) + 10.0
        ts1.forces = np.arange(30).reshape(10, 3) + 100.0

        ts2 = mda.coordinates.TRR.Timestep(10)
        ts2.positions = np.arange(30).reshape(10, 3)
        ts2.velocities = np.arange(30).reshape(10, 3) + 10.0
        ts2.forces = np.arange(30).reshape(10, 3) + 100.0

        self._check_ts(ts1, ts2, "Failed on TRR Timestep")


    def test_wrong_class(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = np.arange(30).reshape(10, 3)

        b = tuple([0, 1, 2, 3])

        assert_(ts1 != b)
        assert_(b != ts1)

    def test_wrong_frame(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = np.arange(30).reshape(10, 3)

        ts2 = mda.coordinates.base.Timestep(10)
        ts2.positions = np.arange(30).reshape(10, 3)
        ts2.frame = 987

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_wrong_n_atoms(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = np.arange(30).reshape(10, 3)

        ts3 = mda.coordinates.base.Timestep(20)

        assert_(ts1 != ts3)
        assert_(ts3 != ts1)

    def test_wrong_pos(self):
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = np.arange(30).reshape(10, 3)

        ts2 = mda.coordinates.base.Timestep(10)
        ts2.positions = np.arange(30).reshape(10, 3) + 1.0

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_vels(self):
        ts1 = mda.coordinates.base.Timestep(10, velocities=True)
        ts2 = mda.coordinates.base.Timestep(10, velocities=True)

        ts1.velocities = np.arange(30).reshape(10, 3)
        ts2.velocities = np.arange(30).reshape(10, 3)

        assert_(ts1 == ts2)
        assert_(ts2 == ts1)

    def test_check_mismatched_vels(self):
        ts1 = mda.coordinates.base.Timestep(10, velocities=True)
        ts2 = mda.coordinates.base.Timestep(10, velocities=False)

        ts1.velocities = np.arange(30).reshape(10, 3)

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_wrong_vels(self):
        ts1 = mda.coordinates.base.Timestep(10, velocities=True)
        ts2 = mda.coordinates.base.Timestep(10, velocities=True)

        ts1.velocities = np.arange(30).reshape(10, 3)
        ts2.velocities = np.arange(30).reshape(10, 3) + 1.0

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_forces(self):
        ts1 = mda.coordinates.base.Timestep(10, forces=True)
        ts2 = mda.coordinates.base.Timestep(10, forces=True)

        ts1.forces = np.arange(30).reshape(10, 3)
        ts2.forces = np.arange(30).reshape(10, 3)

        assert_(ts1 == ts2)
        assert_(ts2 == ts1)

    def test_check_mismatched_forces(self):
        ts1 = mda.coordinates.base.Timestep(10, forces=True)
        ts2 = mda.coordinates.base.Timestep(10, forces=False)

        ts1.forces = np.arange(30).reshape(10, 3)

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_wrong_forces(self):
        ts1 = mda.coordinates.base.Timestep(10, forces=True)
        ts2 = mda.coordinates.base.Timestep(10, forces=True)

        ts1.forces = np.arange(30).reshape(10, 3)
        ts2.forces = np.arange(30).reshape(10, 3) + 1.0

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)


class _TestTimestepInterface(object):
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


class TestCRD(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(XYZ_five, INPCRD)
        self.ts = u.trajectory.ts


class TestDCD(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(PSF, DCD)
        self.ts = u.trajectory.ts


class TestDLPConfig(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(DLP_CONFIG, format='CONFIG')
        self.ts = u.trajectory.ts


class TestDLPHistory(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(DLP_HISTORY, format='HISTORY')
        self.ts = u.trajectory.ts


class TestDMS(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(DMS)
        self.ts = u.trajectory.ts


class TestGMS(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(GMS_ASYMOPT, GMS_ASYMOPT, format='GMS', topology_format='GMS')
        self.ts = u.trajectory.ts


class TestGRO(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(GRO)
        self.ts = u.trajectory.ts


class TestINPCRD(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(XYZ_five, INPCRD)
        self.ts = u.trajectory.ts


class TestLAMMPS(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(LAMMPSdata)
        self.ts = u.trajectory.ts


class TestLAMMPSDCD(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(LAMMPSdata2,LAMMPSdcd2, format='LAMMPS', topology_format='DATA', timeunit='fs')
        self.ts = u.trajectory.ts


class TestMOL2(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(mol2_molecules)
        self.ts = u.trajectory.ts


class TestPDB(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(PDB_small)
        self.ts = u.trajectory.ts


class TestPDBQT(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(PDBQT_input)
        self.ts = u.trajectory.ts


class TestPQR(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(PQR)
        self.ts = u.trajectory.ts


class TestTRJ(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(PRM, TRJ)
        self.ts = u.trajectory.ts


class TestNCDF(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(PRMncdf, NCDF)
        self.ts = u.trajectory.ts


class TestTRR(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(GRO, TRR)
        self.ts = u.trajectory.ts


class TestTRZ(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(TRZ_psf, TRZ)
        self.ts = u.trajectory.ts


class TestXTC(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(GRO, XTC)
        self.ts = u.trajectory.ts


class TestXYZ(_TestTimestepInterface):
    def setUp(self):
        u = self.u = mda.Universe(XYZ_mini)
        self.ts = u.trajectory.ts


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
