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

"""Tests for the API definition of Timestep


_TestTimestepInterface tests the Readers are correctly using Timesteps
"""
import itertools
import numpy as np
from numpy.testing import assert_equal, assert_allclose, assert_array_almost_equal

from MDAnalysis.lib.mdamath import triclinic_vectors
import MDAnalysis as mda
from MDAnalysisTests.datafiles import (PSF, XYZ_five, INPCRD, DCD, DLP_CONFIG,
                                       DLP_HISTORY, DMS, GMS_ASYMOPT, GRO, XTC,
                                       TRR, LAMMPSdata, LAMMPSdata2,
                                       LAMMPSdcd2, mol2_molecules, PDB_small,
                                       PDBQT_input, PQR, PRM, TRJ, PRMncdf,
                                       NCDF, TRZ_psf, TRZ)

from MDAnalysisTests.coordinates.base import assert_timestep_equal, assert_timestep_almost_equal
from MDAnalysis.coordinates.timestep import Timestep
import pytest


class TestTimestep(object):
    """Test all the base functionality of a Timestep

    All Timesteps must pass these tests!

    These test the Timestep independent of the Reader which it
    comes into contact with.  Failures here are the Timesteps fault.
    """
    # define the class made in test
    Timestep = Timestep
    name = "base"  # for error messages only
    size = 10  # size of arrays, 10 is enough to allow slicing etc
    # each coord is unique
    refpos = np.arange(size * 3, dtype=np.float32).reshape(size, 3) * 1.234
    refvel = np.arange(size * 3, dtype=np.float32).reshape(size, 3) * 2.345
    reffor = np.arange(size * 3, dtype=np.float32).reshape(size, 3) * 3.456

    # If you can set box, what the underlying unitcell should be
    # if dimensions are:
    newbox = np.array([10., 11., 12., 90., 90., 90.])
    unitcell = np.array([10., 11., 12., 90., 90., 90.])
    ref_volume = 1320.  # what the volume is after setting newbox
    uni_args = None

    @pytest.fixture()
    def ts(self):
        ts = self.Timestep(self.size)
        ts.frame += 1
        ts.positions = self.refpos
        ts.dimensions = self.newbox

        return ts

    def test_getitem(self, ts):
        assert_equal(ts[1], self.refpos[1])

    def test_getitem_neg(self, ts):
        assert_equal(ts[-1], self.refpos[-1])

    def test_getitem_neg_IE(self, ts):
        with pytest.raises(IndexError):
            ts.__getitem__(-(self.size + 1))


    def test_getitem_pos_IE(self, ts):
        with pytest.raises(IndexError):
            ts.__getitem__((self.size + 1))

    def test_getitem_slice(self, ts):
        assert_equal(len(ts[:2]), len(self.refpos[:2]))
        assert_allclose(ts[:2], self.refpos[:2])

    def test_getitem_slice2(self, ts):
        assert_equal(len(ts[1::2]), len(self.refpos[1::2]))
        assert_allclose(ts[1::2], self.refpos[1::2])

    def test_getitem_ndarray(self, ts):
        sel = np.array([0, 1, 4])
        assert_equal(len(ts[sel]), len(self.refpos[sel]))
        assert_allclose(ts[sel], self.refpos[sel])

    def test_getitem_TE(self, ts):
        with pytest.raises(TypeError):
            ts.__getitem__('string')

    def test_len(self, ts):
        assert_equal(len(ts), self.size)

    def test_iter(self, ts):
        for a, b in zip(ts, self.refpos):
            assert_allclose(a, b)
        assert_equal(len(list(ts)), self.size)

    def test_repr(self, ts):
        assert_equal(type(repr(ts)), str)

    def test_repr_with_box(self, ts):
        assert("with unit cell dimensions" in repr(ts))

    def test_repr_no_box(self, ts):
        ts.dimensions = None
        assert("with unit cell dimensions" not in repr(ts))

    def test_default_dtype_npf32(self, ts):
        assert_equal(ts.dtype, np.float32)


    # Dimensions has 2 possible cases
    # Timestep doesn't do dimensions,
    # returns None for  .dimension and 0 for .volume
    # Timestep does do them, should return values properly
    def test_dimensions(self, ts):
        assert_allclose(ts.dimensions, self.newbox)

    @pytest.mark.parametrize('dtype', (int, np.float32, np.float64))
    def test_dimensions_set_box(self, ts, dtype):
        ts.dimensions = self.newbox.astype(dtype)
        assert ts.dimensions.dtype == np.float32
        assert_allclose(ts.dimensions, self.newbox)

    def test_volume(self, ts):
        ts.dimensions = self.newbox
        assert_equal(ts.volume, self.ref_volume)

    def test_triclinic_vectors(self, ts):
        assert_allclose(ts.triclinic_dimensions,
                        triclinic_vectors(ts.dimensions))

    def test_set_triclinic_vectors(self, ts):
        ref_vec = triclinic_vectors(self.newbox)
        ts.triclinic_dimensions = ref_vec
        assert_equal(ts.dimensions, self.newbox)

    def test_set_dimensions_None(self,ts):
        ts.dimensions = None
        assert(not ts._unitcell.any())

    def test_set_triclinic_dimensions_None(self,ts):
        ts.triclinic_dimensions = None
        assert(not ts._unitcell.any())

    def test_coordinate_getter_shortcuts(self, ts):
        """testing that reading _x, _y, and _z works as expected
        # (Issue 224) (TestTimestep)"""
        assert_equal(ts._x, ts.positions[:, 0])
        assert_equal(ts._y, ts.positions[:, 1])
        assert_equal(ts._z, ts.positions[:, 2])

    def test_coordinate_setter_shortcuts(self, ts):
        # Check that _x _y and _z are read only
        for coordinate in ('_x', '_y', '_z'):
            random_positions = np.arange(self.size).astype(np.float32)
            with pytest.raises(AttributeError):
                setattr(ts, coordinate, random_positions)

    # n_atoms should be a read only property
    # all Timesteps require this attribute
    def test_n_atoms(self, ts):
        assert_equal(ts.n_atoms, self.size)

    def test_n_atoms_readonly(self, ts):
        with pytest.raises(AttributeError):
            ts.__setattr__('n_atoms', 20)

    def test_n_atoms_presence(self, ts):
        assert_equal(hasattr(ts, 'n_atoms'), True)

    def test_unitcell_presence(self, ts):
        assert_equal(hasattr(ts, 'dimensions'), True)

    def test_data_presence(self, ts):
        assert_equal(hasattr(ts, 'data'), True)
        assert_equal(isinstance(ts.data, dict), True)

    def test_allocate_velocities(self, ts):
        assert_equal(ts.has_velocities, False)
        with pytest.raises(mda.NoDataError):
            getattr(ts, 'velocities')

        ts.has_velocities = True
        assert_equal(ts.has_velocities, True)
        assert_equal(ts.velocities.shape, (ts.n_atoms, 3))

    def test_allocate_forces(self, ts):
        assert_equal(ts.has_forces, False)
        with pytest.raises(mda.NoDataError):
            getattr(ts, 'forces')

        ts.has_forces = True
        assert_equal(ts.has_forces, True)
        assert_equal(ts.forces.shape, (ts.n_atoms, 3))

    def test_velocities_remove(self):
        ts = self.Timestep(10, velocities=True)
        ts.frame += 1
        assert_equal(ts.has_velocities, True)

        ts.has_velocities = False
        assert_equal(ts.has_velocities, False)
        with pytest.raises(mda.NoDataError):
            getattr(ts, 'velocities')

    def test_forces_remove(self):
        ts = self.Timestep(10, forces=True)
        ts.frame += 1
        assert_equal(ts.has_forces, True)

        ts.has_forces = False
        assert_equal(ts.has_forces, False)
        with pytest.raises(mda.NoDataError):
            getattr(ts, 'forces')

    def test_check_ts(self):
        with pytest.raises(ValueError):
            self.Timestep.from_coordinates(None, None, None)

    def _from_coords(self, p, v, f):
        posi = self.refpos if p else None
        velo = self.refvel if v else None
        forc = self.reffor if f else None

        ts = self.Timestep.from_coordinates(posi, velo, forc)

        return ts

    @pytest.mark.parametrize('p, v, f', filter(any,
                                               itertools.product([True, False],
                                                                 repeat=3)))
    def test_from_coordinates(self, p, v, f):
        ts = self._from_coords(p, v, f)

        if p:
            assert_array_almost_equal(ts.positions, self.refpos)
        else:
            with pytest.raises(mda.NoDataError):
                getattr(ts, 'positions')
        if v:
            assert_array_almost_equal(ts.velocities, self.refvel)
        else:
            with pytest.raises(mda.NoDataError):
                getattr(ts, 'velocities')
        if f:
            assert_array_almost_equal(ts.forces, self.reffor)
        else:
            with pytest.raises(mda.NoDataError):
                getattr(ts, 'forces')

    def test_from_coordinates_mismatch(self):
        velo = self.refvel[:2]

        with pytest.raises(ValueError):
            self.Timestep.from_coordinates(self.refpos, velo)

    def test_from_coordinates_nodata(self):
        with pytest.raises(ValueError):
            self.Timestep.from_coordinates()

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

    def _check_copy(self, name, ref_ts):
        """Check basic copy"""
        ts2 = ref_ts.copy()

        err_msg = ("Timestep copy failed for format {form}"
                   " on attribute {att}")

        # eq method checks:
        # - frame
        # - n_atoms
        # - positions, vels and forces
        assert ref_ts == ts2

        if not ref_ts.dimensions is None:
            assert_array_almost_equal(ref_ts.dimensions, ts2.dimensions,
                                      decimal=4)
        else:
            assert ref_ts.dimensions == ts2.dimensions

        # Check things not covered by eq
        for d in ref_ts.data:
            assert d in ts2.data
            if isinstance(ref_ts.data[d], np.ndarray):
                assert_array_almost_equal(
                    ref_ts.data[d], ts2.data[d])
            else:
                assert ref_ts.data[d] == ts2.data[d]

    def _check_independent(self, name, ts):
        """Check that copies made are independent"""
        ts2 = ts.copy()

        if ts.has_positions:
            self._check_array(ts.positions, ts2.positions)
        if ts.has_velocities:
            self._check_array(ts.velocities, ts2.velocities)
        if ts.has_forces:
            self._check_array(ts.forces, ts2.forces)
        if not ts.dimensions is None:
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

    def _check_npint_slice(self, name, ts):
        for integers in [np.byte, np.short, np.intc, np.int_, np.longlong,
                         np.intp, np.int8, np.int16, np.int32, np.int64,
                         np.ubyte, np.ushort, np.uintc, np.ulonglong,
                         np.uintp, np.uint8, np.uint16, np.uint32, np.uint64]:
            sl = slice(1, 2, 1)
            ts2 = ts.copy_slice(slice(integers(1), integers(2), integers(1)))
            self._check_slice(ts, ts2, sl)

    def _check_slice(self, ts1, ts2, sl):
        if ts1.has_positions:
            assert_array_almost_equal(ts1.positions[sl], ts2.positions)
        if ts1.has_velocities:
            assert_array_almost_equal(ts1.velocities[sl], ts2.velocities)
        if ts1.has_forces:
            assert_array_almost_equal(ts1.forces[sl], ts2.forces)

    @pytest.mark.parametrize('func', [
        _check_copy,
        _check_independent,
        _check_copy_slice_indices,
        _check_copy_slice_slice,
        _check_npint_slice
    ])
    def test_copy(self, func, ts):
        if self.uni_args is None:
            return
        u = mda.Universe(*self.uni_args)  # pylint: disable=not-an-iterable
        ts = u.trajectory.ts
        func(self, self.name, ts)

    @pytest.fixture(params=filter(any,
                                  itertools.product([True, False], repeat=3)))
    def some_ts(self, request):
        p, v, f = request.param
        return self._from_coords(p, v, f)

    @pytest.mark.parametrize('func', [
        _check_copy,
        _check_independent,
        _check_copy_slice_indices,
        _check_copy_slice_slice,
        _check_npint_slice
    ])
    def test_copy_slice(self, func, some_ts):
        func(self, self.name, some_ts)

    def test_bad_slice(self, some_ts):
        sl = ['this', 'is', 'silly']
        with pytest.raises(TypeError):
            some_ts.copy_slice(sl)

    def test_from_timestep(self, some_ts):
        ts = some_ts
        ts2 = self.Timestep.from_timestep(ts)

        assert_timestep_almost_equal(ts, ts2)

    def _get_pos(self):
        # Get generic reference positions
        return np.arange(30).reshape(10, 3) * 1.234

    @pytest.mark.parametrize('p, v, f', filter(any,
                                               itertools.product([True, False],
                                                                 repeat=3)))
    def test_check_equal(self, p, v, f):
        ts1 = self.Timestep(self.size,
                            positions=p,
                            velocities=v,
                            forces=f)
        ts2 = self.Timestep(self.size,
                            positions=p,
                            velocities=v,
                            forces=f)
        if p:
            ts1.positions = self.refpos.copy()
            ts2.positions = self.refpos.copy()
        if v:
            ts1.velocities = self.refvel.copy()
            ts2.velocities = self.refvel.copy()
        if f:
            ts1.forces = self.reffor.copy()
            ts2.forces = self.reffor.copy()

        assert_timestep_equal(ts1, ts2)

    def test_wrong_class_equality(self):
        ts1 = self.Timestep(self.size)
        ts1.positions = self._get_pos()

        b = tuple([0, 1, 2, 3])

        assert ts1 != b
        assert b != ts1

    def test_wrong_frame_equality(self):
        ts1 = self.Timestep(self.size)
        ts1.positions = self._get_pos()

        ts2 = self.Timestep(self.size)
        ts2.positions = self._get_pos()
        ts2.frame = 987

        assert ts1 != ts2
        assert ts2 != ts1

    def test_wrong_n_atoms_equality(self):
        ts1 = self.Timestep(self.size)
        ts1.positions = self._get_pos()

        ts3 = self.Timestep(self.size * 2)

        assert ts1 != ts3
        assert ts3 != ts1

    def test_wrong_pos_equality(self):
        ts1 = self.Timestep(self.size)
        ts1.positions = self._get_pos()

        ts2 = self.Timestep(self.size)
        ts2.positions = self._get_pos() + 1.0

        assert ts1 != ts2
        assert ts2 != ts1

    def test_no_pos_inequality(self):
        ts1 = self.Timestep(self.size, positions=False)
        ts2 = self.Timestep(self.size)

        assert ts1 != ts2
        assert ts2 != ts1

    def test_check_vels_equality(self):
        ts1 = self.Timestep(self.size, velocities=True)
        ts2 = self.Timestep(self.size, velocities=True)

        ts1.velocities = self._get_pos()
        ts2.velocities = self._get_pos()

        assert ts1 == ts2
        assert ts2 == ts1

    def test_check_mismatched_vels_equality(self):
        ts1 = self.Timestep(self.size, velocities=True)
        ts2 = self.Timestep(self.size, velocities=False)

        ts1.velocities = self._get_pos()

        assert ts1 != ts2
        assert ts2 != ts1

    def test_check_wrong_vels_equality(self):
        ts1 = self.Timestep(self.size, velocities=True)
        ts2 = self.Timestep(self.size, velocities=True)

        ts1.velocities = self._get_pos()
        ts2.velocities = self._get_pos() + 1.0

        assert ts1 != ts2
        assert ts2 != ts1

    def test_check_forces_equality(self):
        ts1 = self.Timestep(self.size, forces=True)
        ts2 = self.Timestep(self.size, forces=True)

        ts1.forces = self._get_pos()
        ts2.forces = self._get_pos()

        assert ts1 == ts2
        assert ts2 == ts1

    def test_check_mismatched_forces_equality(self):
        ts1 = self.Timestep(self.size, forces=True)
        ts2 = self.Timestep(self.size, forces=False)

        ts1.forces = self._get_pos()

        assert ts1 != ts2
        assert ts2 != ts1

    def test_check_wrong_forces_equality(self):
        ts1 = self.Timestep(self.size, forces=True)
        ts2 = self.Timestep(self.size, forces=True)

        ts1.forces = self._get_pos()
        ts2.forces = self._get_pos() + 1.0

        assert ts1 != ts2
        assert ts2 != ts1

    @pytest.mark.parametrize('dim1', [None, [2., 2., 2., 90., 90., 90.]])
    def test_dims_mismatch_inequality(self, dim1):
        ts1 = self.Timestep(self.size)
        ts1.dimensions = dim1
        ts2 = self.Timestep(self.size)
        ts2.dimensions = [1., 1., 1., 90., 90., 90.]

        assert ts1 != ts2
        assert ts2 != ts1


# TODO: Merge this into generic Reader tests
# These tests are all included in BaseReaderTest
# Once Readers use that TestClass, delete this one

class TestBaseTimestepInterface(object):
    """Test the Timesteps created by Readers

    This checks that Readers are creating correct Timestep objects,
    ensuring a consistent interface when using Timestep objects.

    Failures here are the Reader's fault

    See Issue #250 for discussion
    """
    @pytest.fixture(params=(
        (XYZ_five, INPCRD, None, None),
        (PSF, DCD, None, None),
        (DLP_CONFIG, None, 'CONFIG', None),
        (DLP_HISTORY, None, 'HISTORY', None),
        (DMS, None, None, None),
        (GRO, None, None, None),
        (XYZ_five, INPCRD, None, None),
        (LAMMPSdata, None, None, None),
        (mol2_molecules, None, None, None),
        (PDB_small, None, None, None),
        (PQR, None, None, None),
        (PDBQT_input, None, None, None),
        (PRM, TRJ, None, None),
        (GRO, XTC, None, None),
        (TRZ_psf, TRZ, None, None),
        (GRO, TRR, None, None),
        (GMS_ASYMOPT, GMS_ASYMOPT, 'GMS', 'GMS'),
        (LAMMPSdata2, LAMMPSdcd2, 'LAMMPS', 'DATA'),
        (PRMncdf, NCDF, None, None),
    ))
    def universe(self, request):
        topology, trajectory, trajectory_format, topology_format = request.param
        if trajectory_format is not None and topology_format is not None:
            return mda.Universe(topology, trajectory, format=trajectory_format,
                                topology_format=topology_format)

        if trajectory is not None:
            return mda.Universe(topology, trajectory)
        else:
            return mda.Universe(topology, format=trajectory_format)

    def test_frame(self, universe):
        assert_equal(universe.trajectory.ts.frame, 0)

    def test_dt(self, universe):
        assert_equal(universe.trajectory.dt, universe.trajectory.ts.dt)


@pytest.mark.parametrize('uni', [
    [(PSF, DCD), {}],  # base Timestep
    [(DLP_CONFIG,), {'format': 'CONFIG'}],  # DLPoly
    [(DMS,), {}],  # DMS
    [(GRO,), {}],  # GRO
    [(np.zeros((1, 30, 3)),), {}],  # memory
    [(PRM, TRJ), {}],  # TRJ
    [(TRZ_psf, TRZ), {}],  # TRZ
])
def test_atomgroup_dims_access(uni):
    uni_args, uni_kwargs = uni
    # check that AtomGroup.dimensions always returns a copy
    u = mda.Universe(*uni_args, **uni_kwargs, to_guess=())

    ag = u.atoms[:10]

    dims = ag.dimensions
    ts = u.trajectory.ts

    # dimensions from AtomGroup should be equal
    assert_equal(dims, ts.dimensions)
    # but if np array shouldn't be the same object, i.e. it should have been copied
    if dims is not None:
        assert dims is not ts.dimensions


class TimeTimestepNoBox:
    @staticmethod
    @pytest.fixture
    def ts():
        ts = mda.coordinates.base.Timestep(20)

        ts.dimensions = None

        return ts

    def test_Noneness(self, ts):
        assert ts.dimensions is None

    def test_triclinic_Noneness(self, ts):
        assert ts.triclinic_dimensions is None

    def test_set_dimensions(self, ts):
        # None to something
        box = np.array([10, 10, 10, 90, 90, 90])
        ts.dimensions = box

        assert_array_equal(ts.dimensions, box)

    def test_set_triclinic_dimensions(self, ts):
        # None to something via triclinic
        box = np.array([[10, 0, 0], [5, 10, 0], [0, 5, 10]])

        ts.triclinic_dimensions = box

        assert_array_equal(ts.triclinic_dimensions, box)

    def test_set_None(self, ts):
        ts.dimensions = [10, 10, 10, 90, 90, 90]

        assert not ts.dimensions is None

        ts.dimensions = None

        assert ts.dimensions is None
        assert ts.triclinic_dimenions is None

    def test_set_triclinic_None(self, ts):
        ts.dimensions = [10, 10, 10, 90, 90, 90]

        assert not ts.dimensions is None

        ts.triclinic_dimensions = None

        assert ts.dimensions is None
        assert ts.triclinic_dimenions is None

    def test_set_zero_box(self, ts):
        ts.dimensions = np.zeros(6)

        assert ts.dimensions is None
        assert ts.triclinic_dimenions is None

    def test_set_zero_box_triclinic(self, ts):
        ts.triclinic_dimensions = np.zeros((3, 3))

        assert ts.dimensions is None
        assert ts.triclinic_dimenions is None
