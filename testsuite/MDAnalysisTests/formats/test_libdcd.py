# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import pickle

from collections import namedtuple
import os
import sys
import string
import struct
import platform

import hypothesis.strategies as strategies
from hypothesis import example, given
import hypothesis
import numpy as np

from numpy.testing import (assert_array_almost_equal, assert_equal,
                           assert_array_equal, assert_almost_equal)

from MDAnalysis.lib.formats.libdcd import DCDFile, DCD_IS_CHARMM, DCD_HAS_EXTRA_BLOCK

from MDAnalysisTests.datafiles import (
    DCD, DCD_NAMD_TRICLINIC, legacy_DCD_ADK_coords, legacy_DCD_NAMD_coords,
    legacy_DCD_c36_coords, DCD_TRICLINIC)

import pytest


@pytest.mark.parametrize("dcdfile, is_periodic",
                         [(DCD, False), (DCD_NAMD_TRICLINIC, True),
                          (DCD_TRICLINIC, True)])
def test_is_periodic(dcdfile, is_periodic):
    with DCDFile(dcdfile) as f:
        assert f.is_periodic == is_periodic


@pytest.mark.parametrize("dcdfile, natoms", [(DCD, 3341), (DCD_NAMD_TRICLINIC,
                                                           5545),
                                             (DCD_TRICLINIC, 375)])
def test_read_coordsshape(dcdfile, natoms):
    # confirm shape of coordinate data against result from previous
    # MDAnalysis implementation of DCD file handling
    with DCDFile(dcdfile) as dcd:
        dcd_frame = dcd.read()
    xyz = dcd_frame[0]
    assert xyz.shape == (natoms, 3)


@pytest.mark.parametrize(
    "dcdfile, unit_cell",
    [(DCD, [0., 90., 0., 90., 90., 0.]),
     (DCD_NAMD_TRICLINIC, [38.42659378, 0.499563, 38.393102, 0., 0., 44.7598]),
     (DCD_TRICLINIC,
      [30.841836, 14.578635, 31.780088, 9.626323, -2.60815, 32.67009])])
def test_read_unit_cell(dcdfile, unit_cell):
    # confirm unit cell read against result from previous
    # MDAnalysis implementation of DCD file handling
    with DCDFile(dcdfile) as dcd:
        dcd_frame = dcd.read()
    assert_array_almost_equal(dcd_frame.unitcell, unit_cell)


def test_seek_over_max():
    with DCDFile(DCD) as dcd:
        with pytest.raises(EOFError):
            dcd.seek(102)


@pytest.fixture
def dcd():
    with DCDFile(DCD) as dcd:
        yield dcd


def _assert_compare_readers(old_reader, new_reader):
    #  same as next(old_reader)
    frame = old_reader.read()
    #  same as next(new_reader)
    new_frame = new_reader.read()

    assert old_reader.fname == new_reader.fname
    assert old_reader.tell() == new_reader.tell()
    assert_almost_equal(frame.xyz, new_frame.xyz)
    assert_almost_equal(frame.unitcell, new_frame.unitcell)


def test_pickle(dcd):
    mid = len(dcd) // 2
    dcd.seek(mid)
    new_dcd = pickle.loads(pickle.dumps(dcd))
    _assert_compare_readers(dcd, new_dcd)


def test_pickle_last(dcd):
    #  This is the file state when DCDReader is in its last frame.
    #  (Issue #2878)

    dcd.seek(len(dcd) - 1)
    _ = dcd.read()
    new_dcd = pickle.loads(pickle.dumps(dcd))

    assert dcd.fname == new_dcd.fname
    assert dcd.tell() == new_dcd.tell()
    with pytest.raises(StopIteration):
        new_dcd.read()


def test_pickle_closed(dcd):
    dcd.seek(len(dcd) - 1)
    dcd.close()
    new_dcd = pickle.loads(pickle.dumps(dcd))

    assert dcd.fname == new_dcd.fname
    assert dcd.tell() != new_dcd.tell()


def test_pickle_after_read(dcd):
    _ = dcd.read()
    new_dcd = pickle.loads(pickle.dumps(dcd))
    _assert_compare_readers(dcd, new_dcd)


def test_pickle_immediately(dcd):
    new_dcd = pickle.loads(pickle.dumps(dcd))

    assert dcd.fname == new_dcd.fname
    assert dcd.tell() == new_dcd.tell()


@pytest.mark.parametrize("new_frame", (10, 42, 21))
def test_seek_normal(new_frame, dcd):
    # frame seek within range is tested
    dcd.seek(new_frame)
    assert dcd.tell() == new_frame


def test_seek_negative(dcd):
    with pytest.raises(IOError):
        dcd.seek(-78)


def test_iteration(dcd):
    num_iters = 10
    for _ in range(num_iters):
        dcd.__next__()
    assert dcd.tell() == num_iters


def test_open_wrong_mode():
    with pytest.raises(IOError):
        DCDFile(DCD, 'e')


def test_raise_not_existing():
    with pytest.raises(IOError):
        DCDFile('foo')


def test_zero_based_frames_counting(dcd):
    assert dcd.tell() == 0


@pytest.mark.parametrize("dcdfile, natoms", [(DCD, 3341), (DCD_NAMD_TRICLINIC,
                                                           5545),
                                             (DCD_TRICLINIC, 375)])
def test_natoms(dcdfile, natoms):
    with DCDFile(dcdfile) as dcd:
        assert dcd.header['natoms'] == natoms


def test_read_closed(dcd):
    dcd.close()
    with pytest.raises(IOError):
        dcd.read()


@pytest.mark.parametrize("dcdfile, nframes", [(DCD, 98), (DCD_NAMD_TRICLINIC,
                                                          1), (DCD_TRICLINIC,
                                                               10)])
def test_length_traj(dcdfile, nframes):
    with DCDFile(dcdfile) as dcd:
        assert len(dcd) == nframes


def test_read_write_mode_file(tmpdir):
    fname = str(tmpdir.join('foo'))
    with DCDFile(fname, 'w') as f:
        with pytest.raises(IOError):
            f.read()


def test_iterating_twice(dcd):
    with dcd as f:
        for i, _ in enumerate(f):
            assert_equal(i + 1, f.tell())
        # second iteration should work from start again
        for i, _ in enumerate(f):
            assert_equal(i + 1, f.tell())


DCD_HEADER = '''* DIMS ADK SEQUENCE FOR PORE PROGRAM                                            * WRITTEN BY LIZ DENNING (6.2008)                                               *  DATE:     6/ 6/ 8     17:23:56      CREATED BY USER: denniej0                '''
DCD_NAMD_TRICLINIC_HEADER = 'Created by DCD pluginREMARKS Created 06 July, 2014 at 17:29Y5~CORD,'
DCD_TRICLINIC_HEADER = '* CHARMM TRICLINIC BOX TESTING                                                  * (OLIVER BECKSTEIN 2014)                                                       * BASED ON NPTDYN.INP : SCOTT FELLER, NIH, 7/15/95                              * TEST EXTENDED SYSTEM CONSTANT PRESSURE AND TEMPERATURE                        * DYNAMICS WITH WATER BOX.                                                      *  DATE:     7/ 7/14     13:59:46      CREATED BY USER: oliver                  '


@pytest.mark.parametrize("dcdfile, remarks",
                         ((DCD, DCD_HEADER), (DCD_NAMD_TRICLINIC,
                                              DCD_NAMD_TRICLINIC_HEADER),
                          (DCD_TRICLINIC, DCD_TRICLINIC_HEADER)))
def test_header_remarks(dcdfile, remarks):
    # confirm correct header remarks section reading
    with DCDFile(dcdfile) as f:
        assert len(f.header['remarks']) == len(remarks)


@pytest.mark.parametrize("dcdfile, legacy_data, frames",
                         ((DCD, legacy_DCD_ADK_coords, [5, 29]),
                          (DCD_NAMD_TRICLINIC, legacy_DCD_NAMD_coords, [0]),
                          (DCD_TRICLINIC, legacy_DCD_c36_coords, [1, 4])))
def test_read_coord_values(dcdfile, legacy_data, frames):
    # test the actual values of coordinates read in versus
    # stored values read in by the legacy DCD handling framework
    # to reduce repo storage burden, we only compare for a few
    # randomly selected frames
    legacy = np.load(legacy_data)
    with DCDFile(dcdfile) as dcd:
        for index, frame_num in enumerate(frames):
            dcd.seek(frame_num)
            actual_coords = dcd.read()[0]
            desired_coords = legacy[index]
            assert_array_equal(actual_coords, desired_coords)


@pytest.mark.parametrize("dcdfile, legacy_data, frame_idx",
                         ((DCD, legacy_DCD_ADK_coords, [5, 29]),
                          (DCD_NAMD_TRICLINIC, legacy_DCD_NAMD_coords, [0]),
                          (DCD_TRICLINIC, legacy_DCD_c36_coords, [1, 4])))
def test_readframes(dcdfile, legacy_data, frame_idx):
    legacy = np.load(legacy_data)
    with DCDFile(dcdfile) as dcd:
        frames = dcd.readframes()
        xyz = frames.xyz
        assert_equal(len(xyz), len(dcd))
        for index, frame_num in enumerate(frame_idx):
            assert_array_almost_equal(xyz[frame_num], legacy[index])


def test_write_header(tmpdir):
    # test that _write_header() can produce a very crude
    # header for a new / empty file
    testfile = str(tmpdir.join('test.dcd'))
    with DCDFile(testfile, 'w') as dcd:
        dcd.write_header(
            remarks='Crazy!',
            natoms=22,
            istart=12,
            nsavc=10,
            delta=0.02,
            is_periodic=1)

    with DCDFile(testfile) as dcd:
        header = dcd.header
        assert header['remarks'] == 'Crazy!'
        assert header['natoms'] == 22
        assert header['istart'] == 12
        assert header['is_periodic'] == 1
        assert header['nsavc'] == 10
        assert np.allclose(header['delta'], .02)

    # we also check the bytes written directly.
    with open(testfile, 'rb') as fh:
        header_bytes = fh.read()
    # check for magic number
    assert struct.unpack('i', header_bytes[:4])[0] == 84
    # magic number should be written again before remark section
    assert struct.unpack('i', header_bytes[88:92])[0] == 84
    # length of remark section. We hard code this to 244 right now
    assert struct.unpack('i', header_bytes[92:96])[0] == 244
    # say we have 3 block of length 80
    assert struct.unpack('i', header_bytes[96:100])[0] == 3
    # after the remark section the length should be reported again
    assert struct.unpack('i', header_bytes[340:344])[0] == 244
    # this is a magic number as far as I see
    assert struct.unpack('i', header_bytes[344:348])[0] == 4


def test_write_no_header(tmpdir):
    fname = str(tmpdir.join('test.dcd'))
    with DCDFile(fname, 'w') as dcd:
        with pytest.raises(IOError):
            dcd.write(np.ones(3), np.ones(6))


def test_write_header_twice(tmpdir):
    # an IOError should be raised if a duplicate
    # header writing is attempted

    header = {
        "remarks": 'Crazy!',
        "natoms": 22,
        "istart": 12,
        "nsavc": 10,
        "delta": 0.02,
        "is_periodic": 1
    }

    fname = str(tmpdir.join('test.dcd'))
    with DCDFile(fname, 'w') as dcd:
        dcd.write_header(**header)
        with pytest.raises(IOError):
            dcd.write_header(**header)


def test_write_header_wrong_mode(dcd):
    # an exception should be raised on any attempt to use
    # write_header with a DCDFile object in 'r' mode
    with pytest.raises(IOError):
        dcd.write_header(
            remarks='Crazy!',
            natoms=22,
            istart=12,
            nsavc=10,
            delta=0.02,
            is_periodic=1)


def test_write_mode(dcd):
    # ensure that writing of DCD files only occurs with properly
    # opened files
    with pytest.raises(IOError):
        dcd.write(xyz=np.zeros((3, 3)), box=np.zeros(6, dtype=np.float64))


def write_dcd(in_name, out_name, remarks='testing', header=None):
    with DCDFile(in_name) as f_in, DCDFile(out_name, 'w') as f_out:
        if header is None:
            header = f_in.header
        f_out.write_header(**header)
        for frame in f_in:
            f_out.write(xyz=frame.xyz, box=frame.unitcell)


@pytest.mark.xfail((os.name == 'nt'
                    and sys.maxsize <= 2**32) or
                    platform.machine() == 'aarch64',
                   reason="occasional fail on 32-bit windows and ARM")
# occasionally fails due to unreliable test timings
@hypothesis.settings(deadline=None)  # see Issue 3096
@given(remarks=strategies.text(
    alphabet=string.printable, min_size=0,
    max_size=239))  # handle the printable ASCII strings
@example(remarks='')
def test_written_remarks_property(remarks, tmpdir_factory):
    # property based testing for writing of a wide range of string
    # values to REMARKS field
    dcd = DCDFile(DCD)
    dirname = str(id(remarks)) + "_"
    testfile = str(tmpdir_factory.mktemp(dirname).join('test.dcd'))
    header = dcd.header
    header['remarks'] = remarks
    write_dcd(DCD, testfile, header=header)
    expected_remarks = remarks
    with DCDFile(testfile) as f:
        assert f.header['remarks'] == expected_remarks


@pytest.fixture(scope='session')
def written_dcd(tmpdir_factory):
    with DCDFile(DCD) as dcd:
        header = dcd.header
    testfile = tmpdir_factory.mktemp('dcd').join('test.dcd')
    testfile = str(testfile)
    write_dcd(DCD, testfile)
    Result = namedtuple("Result", "testfile, header, orgfile")
    # throw away last char we didn't save due to null termination
    header['remarks'] = header['remarks'][:-1]
    return Result(testfile, header, DCD)


def test_written_header(written_dcd):
    header = written_dcd.header
    with DCDFile(written_dcd.testfile) as dcd:
        dcdheader = dcd.header
        assert dcdheader == header


def test_written_num_frames(written_dcd):
    with DCDFile(written_dcd.testfile) as dcd, DCDFile(
            written_dcd.orgfile) as other:
        assert len(dcd) == len(other)


def test_written_dcd_coordinate_data_shape(written_dcd):
    with DCDFile(written_dcd.testfile) as dcd, DCDFile(
            written_dcd.orgfile) as other:
        for frame, other_frame in zip(dcd, other):
            assert frame.xyz.shape == other_frame.xyz.shape


def test_written_seek(written_dcd):
    # ensure that we can seek properly on written DCD file
    with DCDFile(written_dcd.testfile) as f:
        f.seek(40)
        assert_equal(f.tell(), 40)


def test_written_coord_match(written_dcd):
    with DCDFile(written_dcd.testfile) as test, DCDFile(
            written_dcd.orgfile) as ref:
        for frame, o_frame in zip(test, ref):
            assert_array_almost_equal(frame.xyz, o_frame.xyz)


def test_written_unit_cell(written_dcd):
    with DCDFile(written_dcd.testfile) as test, DCDFile(
            written_dcd.orgfile) as ref:
        for frame, o_frame in zip(test, ref):
            assert_array_almost_equal(frame.unitcell, o_frame.unitcell)


@pytest.mark.parametrize("dtype", (np.int32, np.int64, np.float32, np.float64,
                                   int, float))
def test_write_all_dtypes(tmpdir, dtype):
    fname = str(tmpdir.join('foo.dcd'))
    with DCDFile(fname, 'w') as out:
        natoms = 10
        xyz = np.ones((natoms, 3), dtype=dtype)
        box = np.ones(6, dtype=dtype)
        out.write_header(
            remarks='test',
            natoms=natoms,
            is_periodic=1,
            delta=1,
            nsavc=1,
            istart=1)
        out.write(xyz=xyz, box=box)


@pytest.mark.parametrize("array_like", (np.array, list))
def test_write_array_like(tmpdir, array_like):
    fname = str(tmpdir.join('foo.dcd'))
    with DCDFile(fname, 'w') as out:
        natoms = 10
        xyz = array_like([[1, 1, 1] for i in range(natoms)])
        box = array_like([i for i in range(6)])
        out.write_header(
            remarks='test',
            natoms=natoms,
            is_periodic=1,
            delta=1,
            nsavc=1,
            istart=1)
        out.write(xyz=xyz, box=box)


def test_write_wrong_shape_xyz(tmpdir):
    fname = str(tmpdir.join('foo.dcd'))
    with DCDFile(fname, 'w') as out:
        natoms = 10
        xyz = np.ones((natoms + 1, 3))
        box = np.ones(6)
        out.write_header(
            remarks='test',
            natoms=natoms,
            is_periodic=1,
            delta=1,
            nsavc=1,
            istart=1)
        with pytest.raises(ValueError):
            out.write(xyz=xyz, box=box)


def test_write_wrong_shape_box(tmpdir):
    fname = str(tmpdir.join('foo.dcd'))
    with DCDFile(fname, 'w') as out:
        natoms = 10
        xyz = np.ones((natoms, 3))
        box = np.ones(7)
        out.write_header(
            remarks='test',
            natoms=natoms,
            is_periodic=1,
            delta=1,
            nsavc=1,
            istart=1)
        with pytest.raises(ValueError):
            out.write(xyz=xyz, box=box)


@pytest.mark.parametrize("dcdfile", (DCD, DCD_TRICLINIC, DCD_NAMD_TRICLINIC))
def test_relative_frame_sizes(dcdfile):
    # the first frame of a DCD file should always be >= in size
    # to subsequent frames, as the first frame contains the same
    # atoms + (optional) fixed atoms
    with DCDFile(dcdfile) as dcd:
        first_frame_size = dcd._firstframesize
        general_frame_size = dcd._framesize

    assert first_frame_size >= general_frame_size


@pytest.mark.parametrize("dcdfile", (DCD, DCD_TRICLINIC, DCD_NAMD_TRICLINIC))
def test_file_size_breakdown(dcdfile):
    # the size of a DCD file is equivalent to the sum of the header
    # size, first frame size, and (N - 1 frames) * size per general
    # frame

    expected = os.path.getsize(dcdfile)
    with DCDFile(dcdfile) as dcd:
        actual = dcd._header_size + dcd._firstframesize + (
            (dcd.n_frames - 1) * dcd._framesize)
    assert actual == expected


@pytest.mark.parametrize("dcdfile", (DCD, DCD_TRICLINIC, DCD_NAMD_TRICLINIC))
def test_nframessize_int(dcdfile):
    # require that the (nframessize / framesize) value used by DCDFile
    # is an integer (because nframessize / framesize + 1 = total frames,
    # which must also be an int)
    filesize = os.path.getsize(dcdfile)
    with DCDFile(dcdfile) as dcd:
        nframessize = filesize - dcd._header_size - dcd._firstframesize
        assert float(nframessize) % float(dcd._framesize) == 0


@pytest.mark.parametrize(
    "slice, length",
    [([None, None, None], 98), ([0, None, None], 98), ([None, 98, None], 98),
     ([None, None, 1], 98), ([None, None, -1], 98), ([2, 6, 2], 2),
     ([0, 10, None], 10), ([2, 10, None], 8), ([0, 1, 1], 1), ([1, 1, 1], 0),
     ([1, 2, 1], 1), ([1, 2, 2], 1), ([1, 4, 2], 2), ([1, 4, 4], 1), ([
         0, 5, 5
     ], 1), ([3, 5, 1], 2), ([4, 0, -1], 4), ([5, 0, -2], 3), ([5, 0, -4], 2)])
def test_readframes_slices(slice, length, dcd):
    start, stop, step = slice
    allframes = dcd.readframes().xyz
    frames = dcd.readframes(start=start, stop=stop, step=step)
    xyz = frames.xyz
    assert len(xyz) == length
    assert_array_almost_equal(xyz, allframes[start:stop:step])


@pytest.mark.parametrize("order, shape", (
    ('fac', (98, 3341, 3)),
    ('fca', (98, 3, 3341)),
    ('afc', (3341, 98, 3)),
    ('acf', (3341, 3, 98)),
    ('caf', (3, 3341, 98)),
    ('cfa', (3, 98, 3341)), ))
def test_readframes_order(order, shape, dcd):
    x = dcd.readframes(order=order).xyz
    assert x.shape == shape


@pytest.mark.parametrize("indices", [[1, 2, 3, 4], [5, 10, 15, 19],
                                     [9, 4, 2, 0, 50]])
def test_readframes_atomindices(indices, dcd):
    allframes = dcd.readframes(order='afc').xyz
    frames = dcd.readframes(indices=indices, order='afc')
    xyz = frames.xyz
    assert len(xyz) == len(indices)
    assert_array_almost_equal(xyz, allframes[indices])


def test_write_random_unitcell(tmpdir):
    testname = str(tmpdir.join('test.dcd'))
    rstate = np.random.RandomState(1178083)
    random_unitcells = rstate.uniform(high=80, size=(98, 6)).astype(np.float64)

    with DCDFile(DCD) as f_in, DCDFile(testname, 'w') as f_out:
        header = f_in.header
        header['is_periodic'] = True
        f_out.write_header(**header)
        for index, frame in enumerate(f_in):
            f_out.write(xyz=frame.xyz, box=random_unitcells[index])

    with DCDFile(testname) as test:
        for index, frame in enumerate(test):
            assert_array_almost_equal(frame.unitcell, random_unitcells[index])
