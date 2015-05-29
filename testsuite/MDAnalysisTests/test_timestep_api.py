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

"""Tests for the API definition of Timestep, followed by checks for each subclass

"""

import itertools
import numpy as np
from numpy.testing import (TestCase, assert_raises, assert_equal, assert_allclose,
                           assert_array_almost_equal)
from nose.plugins.attrib import attr

from . import knownfailure

import MDAnalysis as mda

from .datafiles import PSF, DCD, DCD_empty, PDB_small, XPDB_small, PDB_closed, PDB_multiframe, \
    PDB, CRD, XTC, TRR, GRO, DMS, CONECT, \
    XYZ, XYZ_bz2, XYZ_psf, PRM, TRJ, TRJ_bz2, PRMpbc, TRJpbc_bz2, PRMncdf, NCDF, PQR, \
    PDB_sub_dry, TRR_sub_sol, PDB_sub_sol, TRZ, TRZ_psf, LAMMPSdata, LAMMPSdata_mini, \
    PSF_TRICLINIC, DCD_TRICLINIC, PSF_NAMD_TRICLINIC, DCD_NAMD_TRICLINIC, \
    GMS_ASYMOPT, GMS_SYMOPT, GMS_ASYMSURF, XYZ_mini, PFncdf_Top, PFncdf_Trj, \
    INPCRD, XYZ_five


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


class _xdrCoreTimestep(_BaseTimestep):
    Timestep = mda.coordinates.xdrfile.core.Timestep
    name = "xdrCore"
    has_box = True
    set_box = True
    unitcell = np.array([[10., 0., 0.],
                         [0., 11., 0.],
                         [0., 0., 12.]])


class _TRRTimestep(_xdrCoreTimestep):
    Timestep = mda.coordinates.TRR.Timestep
    name = "TRR"


class _XTCTimestep(_xdrCoreTimestep):
    Timestep = mda.coordinates.XTC.Timestep
    name = "XTC"


class _TestTimestep(TestCase):
    """Test all the base functionality of a Timestep

    All Timesteps must pass these tests!
    """

    def setUp(self):
        self.ts = self.Timestep(self.size)
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

    def test_copy(self):
        ts2 = self.ts.copy()
        # no __eq__ method
        assert_allclose(ts2._pos, self.ts._pos)
        assert_equal(self.ts is ts2, False)

    def _test_TS_slice(self, ref_TS, TS2, sel):
        per_atom = ['_x', '_y', '_z', '_pos', '_velocities', '_forces']
        ignore = ['numatoms', '_numatoms']

        for att in ref_TS.__dict__:
            try:
                if att in per_atom:
                    assert_equal(ref_TS.__dict__[att][sel], TS2.__dict__[att],
                                 err_msg="Timestep slice failed for format: '%s' on attribute: '%s'"
                                         % (self.name, att))
                elif not att in ignore:
                    assert_equal(ref_TS.__dict__[att], TS2.__dict__[att],
                                 err_msg="Timestep slice failed for format: '%s' on attribute: '%s'"
                                         % (self.name, att))
            except KeyError:
                self.fail("Timestep copy failed for format: '%s' on attribute: '%s'"
                          % (self.name, att))

    def test_copy_slice(self):
        sel = slice(0, self.size, 2)
        ts2 = self.ts.copy_slice(sel)

        self._test_TS_slice(self.ts, ts2, sel)

    def test_copy_slice2(self):
        sel = [0, 1, 3]
        ts2 = self.ts.copy_slice(sel)

        self._test_TS_slice(self.ts, ts2, sel)

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

    @attr('issue')
    def test_coordinate_getter_shortcuts(self):
        """testing that reading _x, _y, and _z works as expected (Issue 224) (TestTimestep)"""
        assert_allclose(self.ts._x, self.ts._pos[:,0])
        assert_allclose(self.ts._y, self.ts._pos[:,1])
        assert_allclose(self.ts._z, self.ts._pos[:,2])

    @attr('issue')
    @knownfailure("Timesteps without setter decorators break the _x, _y, _z views onto"
                  " _pos (Issue 213)", AssertionError, mightpass=True)
    def test_coordinate_setter_shortcuts(self):
        # Check that writing to _x, _y, and _z works as expected (Issues 224 and 213)#
        # For TRRs setting _x, _y, and _z works because the assignment is managed by a decorator.
        # For other formats that don't have such decorators assigning to any of these silently breaks
        #  their being a view of the _pos array. Hopefuly all gets clean after Issue 213 is addressed
        #  and this test can be moved to the general Timestep test class.
        for idx, coordinate in enumerate(('_x', '_y', '_z')):
            random_positions = np.random.random(self.size).astype(np.float32)
            setattr(self.ts, coordinate, random_positions)
            assert_allclose(getattr(self.ts, coordinate), random_positions)
            assert_allclose(self.ts._pos[:,idx], random_positions)

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


class TestxdrCoreTimestep(_TestTimestep, _xdrCoreTimestep):
    pass


class TestXTCTimestep(_TestTimestep, _XTCTimestep):
    pass


class TestTRRTimestep(_TestTimestep, _TRRTimestep):
    pass

class TestTimestep_Copy(TestCase):
    """
    Timestep.copy() method seems to be broken, (Issue 164).  The base.Timestep .copy() method returns a TS of
    class base.Timestep rather than the appropriate subclass.

    This class makes a TS object of the first frame, .copy()'s this as a new object and compares the content
    of the two resulting objects.

    This test class is then subclassed below to try and test all Timestep classes that exist within MDA.
    """

    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.name = 'DCD (base)'

    def tearDown(self):
        del self.universe
        del self.name

    def test_TS_copy(self):
        """
        Checks equality between two Timesteps

        Will check that TS2 has all the same attributes and values for these attributes as ref_TS.
        """
        ref_TS = self.universe.trajectory.ts
        TS2 = ref_TS.copy()

        err_msg = ("Timestep copy failed for format {form}"
                   " on attribute {att}") 
        
        for att in ref_TS.__dict__:
            ref = ref_TS.__dict__[att]

            try:
                if isinstance(ref, np.ndarray):
                    assert_array_almost_equal(ref, TS2.__dict__[att], prec=4, err_msg=err_msg.format(form=self.name, att=att))
                else:
                    assert_equal(ref, TS2.__dict__[att], err_msg=err_msg.format(form=self.name, att=att))
            except KeyError:
                self.fail(err_msg.format(form=self.name, att=att))

    def test_TS_slice(self):
        ref_TS = self.universe.trajectory.ts

        sel = slice(0, 100, 4)
        TS2 = ref_TS.copy_slice(sel)

        self._test_TS_slice(ref_TS, TS2, sel)

    def test_TS_indices(self):
        ref_TS = self.universe.trajectory.ts

        sel = [0, 1, 2, 3, 5, 8, 13, 21, 34, 55]
        TS2 = ref_TS.copy_slice(sel)

        self._test_TS_slice(ref_TS, TS2, sel)

    def _test_TS_slice(self, ref_TS, TS2, sel):
        per_atom = [
            '_x', '_y', '_z', '_pos', '_velocities', '_forces',
            '_tpos', '_tvelocities', '_tforces']
        ignore = ['numatoms', '_numatoms']

        for att in ref_TS.__dict__:
            try:
                if att in per_atom:
                    assert_equal(ref_TS.__dict__[att][sel], TS2.__dict__[att],
                                 err_msg="Timestep slice failed for format: '%s' on attribute: '%s'"
                                         % (self.name, att))
                elif not att in ignore:
                    assert_equal(ref_TS.__dict__[att], TS2.__dict__[att],
                                 err_msg="Timestep slice failed for format: '%s' on attribute: '%s'"
                                         % (self.name, att))
            except KeyError:
                self.fail("Timestep copy failed for format: '%s' on attribute: '%s'"
                          % (self.name, att))


class TestTimestep_Copy_DMS(TestTimestep_Copy):
    def setUp(self):
        self.universe = mda.Universe(DMS)
        self.name = 'DMS'


class TestTimestep_Copy_GRO(TestTimestep_Copy):
    def setUp(self):
        self.universe = mda.Universe(GRO)
        self.name = 'GRO'


class TestTimestep_Copy_PDB(TestTimestep_Copy):
    def setUp(self):
        self.universe = mda.Universe(PDB_small)
        self.name = 'PDB'


class TestTimestep_Copy_TRJ(TestTimestep_Copy):
    def setUp(self):
        self.universe = mda.Universe(PRM, TRJ)
        self.name = 'TRJ'


class TestTimestep_Copy_TRR(TestTimestep_Copy):
    def setUp(self):
        self.universe = mda.Universe(GRO, TRR)
        self.name = 'TRR'


class TestTimestep_Copy_TRZ(TestTimestep_Copy):
    def setUp(self):
        self.universe = mda.Universe(TRZ_psf, TRZ)
        self.name = 'TRZ'


class TestTimestep_Copy_XTC(TestTimestep_Copy):
    def setUp(self):
        self.universe = mda.Universe(PDB, XTC)
        self.name = 'XTC'


