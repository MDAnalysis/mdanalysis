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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import numpy as np
try:
    from numpy import shares_memory
except ImportError:
    shares_memory = False

from numpy.testing import assert_equal, assert_almost_equal
import pytest
import MDAnalysis as mda

from MDAnalysisTests.datafiles import (
    AUX_XVG,
    AUX_EDR,
    CRD,
    DCD,
    DMS,
    DLP_CONFIG,
    DLP_HISTORY,
    FHIAIMS,
    INPCRD,
    GMS_ASYMOPT,
    GRO,
    GSD,
    LAMMPSdata_mini,
    mol2_molecules,
    MMTF,
    NAMDBIN,
    NCDF,
    PDB_small,
    PDBQT_input,
    PQR,
    TRR,
    TRJ,
    TRZ,
    TXYZ,
    XTC,
    XPDB_small,
    XYZ_mini,
)
from MDAnalysis.coordinates.core import get_reader_for
from MDAnalysis.coordinates.GSD import HAS_GSD
from MDAnalysis.auxiliary.EDR import HAS_PYEDR


@pytest.fixture(params=[
    # formatname, filename
    ('CRD', CRD, dict()),
    ('DATA', LAMMPSdata_mini, dict(n_atoms=1)),
    ('DCD', DCD, dict()),
    ('DMS', DMS, dict()),
    ('CONFIG', DLP_CONFIG, dict()),
    ('HISTORY', DLP_HISTORY, dict()),
    ('INPCRD', INPCRD, dict()),
    ('GMS', GMS_ASYMOPT, dict()),
    ('GRO', GRO, dict()),
    ('MMTF', MMTF, dict()),
    ('MOL2', mol2_molecules, dict()),
    ('PDB', PDB_small, dict()),
    ('PQR', PQR, dict()),
    ('PDBQT', PDBQT_input, dict()),
    ('TRR', TRR, dict()),
    ('TRZ', TRZ, dict(n_atoms=8184)),
    ('TRJ', TRJ, dict(n_atoms=252)),
    ('XTC', XTC, dict()),
    ('XPDB', XPDB_small, dict()),
    ('XYZ', XYZ_mini, dict()),
    ('NCDF', NCDF, dict()),
    ('memory', np.arange(60).reshape(2, 10, 3).astype(np.float64), dict()),
    ('CHAIN', [GRO, GRO, GRO], dict()),
    ('FHIAIMS', FHIAIMS, dict()),
    pytest.param(
        ('GSD', GSD, dict()),
        marks=pytest.mark.skipif(not HAS_GSD, reason='gsd not installed')
    ),
    ('NAMDBIN', NAMDBIN, dict()),
    ('TXYZ', TXYZ, dict()),
])
def ref_reader(request):
    fmt_name, filename, extras = request.param

    r = get_reader_for(filename, format=fmt_name)(filename, **extras)
    try:
        yield r
    finally:
        # make sure file handle is closed afterwards
        r.close()


@pytest.fixture(params=[
    # formatname, filename
    ('CRD', CRD, dict()),
    ('DATA', LAMMPSdata_mini, dict(n_atoms=1)),
    ('DCD', DCD, dict()),
    ('DMS', DMS, dict()),
    ('CONFIG', DLP_CONFIG, dict()),
    ('HISTORY', DLP_HISTORY, dict()),
    ('INPCRD', INPCRD, dict()),
    ('GMS', GMS_ASYMOPT, dict()),
    ('GRO', GRO, dict()),
    ('MMTF', MMTF, dict()),
    ('MOL2', mol2_molecules, dict()),
    ('PDB', PDB_small, dict()),
    ('PQR', PQR, dict()),
    ('PDBQT', PDBQT_input, dict()),
    ('TRR', TRR, dict()),
    ('TRZ', TRZ, dict(n_atoms=8184)),
    ('TRJ', TRJ, dict(n_atoms=252)),
    ('XTC', XTC, dict()),
    ('XPDB', XPDB_small, dict()),
    ('XYZ', XYZ_mini, dict()),
    ('NCDF', NCDF, dict(mmap=False)),
    ('memory', np.arange(60).reshape(2, 10, 3).astype(np.float64), dict()),
    ('CHAIN', [GRO, GRO, GRO], dict()),
    ('FHIAIMS', FHIAIMS, dict()),
    pytest.param(
        ('GSD', GSD, dict()),
        marks=pytest.mark.skipif(not HAS_GSD, reason='gsd not installed')
    ),
    ('NAMDBIN', NAMDBIN, dict()),
    ('TXYZ', TXYZ, dict()),
])
def ref_reader_extra_args(request):
    fmt_name, filename, extras = request.param

    r = get_reader_for(filename, format=fmt_name)(
            filename, convert_units=False, dt=2, time_offset=10,
            foo="bar", **extras)
    try:
        yield r
    finally:
        # make sure file handle is closed afterwards
        r.close()


@pytest.fixture()
def original_and_copy(ref_reader):
    new = ref_reader.copy()
    try:
        yield ref_reader, new
    finally:
        new.close()


@pytest.fixture()
def original_and_copy_extra_args(ref_reader_extra_args):
    new = ref_reader_extra_args.copy()
    try:
        yield ref_reader_extra_args, new
    finally:
        new.close()


def test_reader_n_atoms(original_and_copy):
    original, copy = original_and_copy
    assert original.n_atoms == copy.n_atoms


def test_reader_filename(original_and_copy):
    original, copy = original_and_copy
    assert original.filename == copy.filename


def test_reader_copied_extra_attributes(original_and_copy_extra_args):
    # Issue #3664
    original, copy = original_and_copy_extra_args

    # memory reader subclass protoreader directly and
    # therefore don't have convert_units or _ts_kwargs
    if original.__class__.__bases__[0].__name__ != "ProtoReader":
        assert original.format is not 'MEMORY'
        assert original.convert_units is False
        assert copy.convert_units is False
        assert original._ts_kwargs['time_offset'] == 10
        assert copy._ts_kwargs['time_offset'] == 10
        assert original._ts_kwargs['dt'] == 2
        assert copy._ts_kwargs['dt'] == 2

    assert original.ts.data['time_offset'] == 10
    assert copy.ts.data['time_offset'] == 10

    # Issue #3689 XTC and XDR overwrite `dt`
    if original.format not in ('XTC', 'TRR'):
        assert original.ts.data['dt'] == 2
        assert copy.ts.data['dt'] == 2

    assert copy._kwargs['foo'] == 'bar'

    # checking that non-base attributes are also copied (netcdf reader)
    if hasattr(original, "_mmap"):
        assert original._mmap is False
        assert copy._mmap is False


def test_reader_independent_iteration(original_and_copy):
    # check that the two Readers iterate independently
    original, copy = original_and_copy
    if len(original) < 2:
        pytest.skip('Single frame reader')
    # initially at same frame
    assert original.ts.frame == copy.ts.frame

    copy[1]
    assert original.ts.frame == 0
    assert copy.ts.frame == 1


def test_reader_initial_frame_maintained(original_and_copy):
    # check that copy inherits nonzero frame of original
    original, _ = original_and_copy

    if len(original) < 2:
        pytest.skip('Single frame reader')

    # seek
    original[1]

    copy = original.copy()

    assert original.ts.frame == copy.ts.frame
    assert_equal(original.ts.positions, copy.ts.positions)

def test_reader_initial_next(original_and_copy):
    # check that the filehandle (or equivalent) in the copied Reader
    # is identical to the original
    # ie calling next on both should produce identical results
    original, _ = original_and_copy

    if len(original) < 3:
        pytest.skip('Requires 3 frames')

    original[1]

    copy = original.copy()

    assert original.ts.frame == copy.ts.frame
    assert_equal(original.ts.positions, copy.ts.positions)

    original.next()
    copy.next()
    
    assert original.ts.frame == copy.ts.frame
    assert_equal(original.ts.positions, copy.ts.positions)


def test_timestep_copied(ref_reader):
    # modify the positions and dimensions from original
    # then check that the copy gets this (even though not in file)
    ref_reader.ts.positions *= 2
    ref_reader.ts.dimensions = newbox = [4, 5, 6, 70, 80, 90]
    new = ref_reader.copy()

    assert_equal(ref_reader.ts.positions, new.ts.positions)
    assert_almost_equal(new.ts.dimensions, newbox, decimal=4)
    assert ref_reader.ts.positions.dtype == np.float32
    assert new.ts.positions.dtype == np.float32


@pytest.mark.skipif(shares_memory == False,
                    reason='old numpy lacked shares_memory')
def test_positions_share_memory(original_and_copy):
    # check that the memory in Timestep objects is unique
    original, copy = original_and_copy
    assert not np.shares_memory(original.ts.positions, copy.ts.positions)

    original.ts.positions *= 2

    with pytest.raises(AssertionError):
        assert_equal(original.ts.positions, copy.ts.positions)


@pytest.mark.skipif(HAS_PYEDR, reason="pyedr present")
def test_copy_with_auxiliary_no_pyedr():
    # Check that AuxReaders are copied when reader is copied
    u = mda.Universe(XYZ_mini)
    u.trajectory.add_auxiliary(AUX_XVG, "myaux")
    reader = u.trajectory
    copy = reader.copy()
    for auxname in reader._auxs:
        assert reader._auxs[auxname] == copy._auxs[auxname]
        assert reader._auxs[auxname] is not copy._auxs[auxname]


@pytest.mark.skipif(not HAS_PYEDR, reason="pyedr not installed")
def test_copy_with_auxiliary_pyedr():
    # Check that AuxReaders are copied when reader is copied
    u = mda.Universe(XYZ_mini)
    u.trajectory.add_auxiliary(AUX_XVG, "myaux")
    u.trajectory.add_auxiliary(AUX_EDR, {"1": "Bond", "2": "Angle"})

    reader = u.trajectory
    copy = reader.copy()
    for auxname in reader._auxs:
        assert reader._auxs[auxname] == copy._auxs[auxname]
        assert reader._auxs[auxname] is not copy._auxs[auxname]
