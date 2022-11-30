# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
import pickle

import numpy as np

from numpy.testing import (assert_almost_equal, assert_array_almost_equal,
                           assert_array_equal, assert_equal)

from MDAnalysis.lib.formats.libmdaxdr import TRRFile, XTCFile

from MDAnalysisTests.datafiles import TRR_multi_frame, XTC_multi_frame

import pytest

XTC_OFFSETS = np.array([0, 104, 208, 312, 416, 520, 624, 728, 832, 936])
TRR_OFFSETS = np.array([0, 480, 960, 1440, 1920, 2400, 2880, 3360, 3840, 4320])


@pytest.fixture
def xtc():
    with XTCFile(XTC_multi_frame) as f:
        yield f


@pytest.fixture
def trr():
    with TRRFile(TRR_multi_frame) as f:
        yield f


@pytest.mark.parametrize('fname, xdr', ((XTC_multi_frame, XTCFile),
                                        (TRR_multi_frame, TRRFile)),
                         indirect=True)
class TestCommonAPI(object):
    @staticmethod
    @pytest.fixture
    def fname(request):
        return request.param

    @staticmethod
    @pytest.fixture
    def xdr(request):
        return request.param

    @staticmethod
    @pytest.fixture
    def reader(fname, xdr):
        with xdr(fname) as f:
            yield f

    def test_natoms(self, reader):
        assert reader.n_atoms == 10

    def test_over_seek(self, reader):
        with pytest.raises(EOFError):
            reader.seek(100)

    def test_read_after_iteration(self, reader):
        for frame in reader:
            pass
        with pytest.raises(EOFError):
            reader.read()

    def test_zero_based_frame_tell(self, reader):
        assert reader.tell() == 0

    def test_read_closed(self, reader):
        reader.close()
        with pytest.raises(IOError):
            reader.read()

    def test_iter(self, reader):
        for i, frame in enumerate(reader):
            assert reader.tell() == i + 1

    def test_double_iter(self, reader):
        for i, frame in enumerate(reader):
            assert reader.tell() == i + 1
        for i, frame in enumerate(reader):
            assert reader.tell() == i + 1

    def test_seek(self, reader):
        reader.seek(4)
        assert reader.tell() == 4

    def test_seek_negative(self, reader):
        with pytest.raises(IOError):
            reader.seek(-4)

    def test_len(self, reader):
        assert len(reader) == 10
        # ensure we didn't go to another frame
        assert reader.tell() == 0

    def test_raise_not_existing(self, xdr, fname):
        with pytest.raises(IOError):
            xdr('foo')

    def test_open_wrong_mode(self, xdr, fname):
        with pytest.raises(IOError):
            xdr('foo', 'e')

    def test_read_write_mode_file(self, xdr, tmpdir, fname):
        fname = str(tmpdir.join('foo'))
        with xdr(fname, 'w') as f:
            with pytest.raises(IOError):
                f.read()

    @staticmethod
    def _assert_compare_readers(old_reader, new_reader):
        frame = old_reader.read()
        new_frame = new_reader.read()

        assert old_reader.fname == new_reader.fname
        assert old_reader.tell() == new_reader.tell()

        assert_equal(old_reader.offsets, new_reader.offsets)
        assert_almost_equal(frame.x, new_frame.x)
        assert_almost_equal(frame.box, new_frame.box)
        assert frame.step == new_frame.step
        assert_almost_equal(frame.time, new_frame.time)

    def test_pickle(self, reader):
        mid = len(reader) // 2
        reader.seek(mid)
        new_reader = pickle.loads(pickle.dumps(reader))
        self._assert_compare_readers(reader, new_reader)

    def test_pickle_last_frame(self, reader):
        #  This is the file state when XDRReader is in its last frame.
        #  (Issue #2878)
        reader.seek(len(reader) - 1)
        _ = reader.read()
        new_reader = pickle.loads(pickle.dumps(reader))

        assert reader.fname == new_reader.fname
        assert reader.tell() == new_reader.tell()
        with pytest.raises(StopIteration):
            new_reader.read()

    def test_pickle_closed(self, reader):
        reader.seek(len(reader) - 1)
        reader.close()
        new_reader = pickle.loads(pickle.dumps(reader))

        assert reader.fname == new_reader.fname
        assert reader.tell() != new_reader.tell()

    def test_pickle_immediately(self, reader):
        new_reader = pickle.loads(pickle.dumps(reader))

        assert reader.fname == new_reader.fname
        assert reader.tell() == new_reader.tell()



@pytest.mark.parametrize("xdrfile, fname, offsets",
                         ((XTCFile, XTC_multi_frame, XTC_OFFSETS),
                          (TRRFile, TRR_multi_frame, TRR_OFFSETS)),
                         indirect=True)
class TestOffsets(object):
    @staticmethod
    @pytest.fixture
    def fname(request):
        return request.param

    @staticmethod
    @pytest.fixture
    def xdrfile(request):
        return request.param

    @staticmethod
    @pytest.fixture
    def offsets(request):
        return request.param

    @staticmethod
    @pytest.fixture
    def reader(xdrfile, fname):
        with xdrfile(fname) as f:
            yield f

    def test_offset(self, reader, offsets):
        assert_array_equal(reader.offsets, offsets)

    def test_set_offsets(self, reader, offsets):
        reader.set_offsets(np.arange(len(offsets)))
        assert_array_equal(reader.offsets, np.arange(len(offsets)))

    def test_fail_wrong_offsets(self, reader, offsets):
        reader.set_offsets(np.arange(len(offsets)))
        reader.seek(5)
        with pytest.raises(IOError):
            reader.read()

    def test_seek_overflow(self, reader, offsets):
        with pytest.raises(OverflowError):
            reader._bytes_seek(2**65)

    def test_seek_tell_all(self, reader, offsets):
        for frame, offset in enumerate(offsets):
            reader.seek(frame)
            assert reader._bytes_tell() == offset

    def test_seek_tell_all_bytes(self, reader, offsets):
        for offset in offsets:
            reader._bytes_seek(offset)
            assert reader._bytes_tell() == offset

    def test_seek_tell_largefile(self, reader, offsets):
        # Seeking/telling can be done on offsets larger than the file.
        # Filesize won't change unless a write is done at the offset.
        # We attempt to go just beyond a 4GB limit.
        big_offset = 2**32 + 1000
        # It should be safe to leave the read head so far away
        # because the file is open read-only.
        reader._bytes_seek(big_offset)
        assert reader._bytes_tell() == big_offset


@pytest.mark.parametrize("xdrfile, fname", ((XTCFile, XTC_multi_frame),
                                            (TRRFile, TRR_multi_frame)))
def test_steps(xdrfile, fname):
    with xdrfile(fname) as f:
        for i, frame in enumerate(f):
            assert frame.step == i


@pytest.mark.parametrize("xdrfile, fname", ((XTCFile, XTC_multi_frame),
                                            (TRRFile, TRR_multi_frame)))
def test_time(xdrfile, fname):
    with xdrfile(fname) as f:
        for i, frame in enumerate(f):
            assert frame.time == i * .5


def test_box_xtc(xtc):
    box = np.eye(3) * 20
    for frame in xtc:
        assert_array_almost_equal(frame.box, box, decimal=3)


def test_xyz_xtc(xtc):
    ones = np.ones(30).reshape(10, 3)
    for i, frame in enumerate(xtc):
        assert_array_almost_equal(frame.x, ones * i, decimal=3)


def test_box_trr(trr):
    box = np.eye(3) * 20
    for frame in trr:
        assert_array_almost_equal(frame.box, box)


def test_xyz_trr(trr):
    ones = np.ones(30).reshape(10, 3)
    for i, frame in enumerate(trr):
        assert_array_almost_equal(frame.x, ones * i)


def test_velocities_trr(trr):
    ones = np.ones(30).reshape(10, 3)
    for i, frame in enumerate(trr):
        assert_array_almost_equal(frame.v, ones * i + 10)


def test_forces_trr(trr):
    ones = np.ones(30).reshape(10, 3)
    for i, frame in enumerate(trr):
        assert_array_almost_equal(frame.f, ones * i + 20)


def test_lmbda_trr(trr):
    for i, frame in enumerate(trr):
        assert_almost_equal(frame.lmbda, .01 * i)


@pytest.fixture
def written_xtc(tmpdir, xtc):
    fname = str(tmpdir.join("foo.xtc"))
    with XTCFile(fname, 'w') as f:
        for frame in xtc:
            f.write(*frame)
    with XTCFile(fname) as f:
        yield f


def test_written_step_xtc(written_xtc):
    for i, frame in enumerate(written_xtc):
        assert frame.step == i


def test_written_time_xtc(written_xtc):
    for i, frame in enumerate(written_xtc):
        assert frame.time == i * .5


def test_written_prec_xtc(written_xtc):
    for frame in written_xtc:
        assert frame.prec == 1000.0


def test_written_box_xtc(written_xtc):
    box = np.eye(3) * 20
    for frame in written_xtc:
        assert_array_almost_equal(frame.box, box, decimal=3)


def test_written_xyx_xtc(written_xtc):
    ones = np.ones(30).reshape(10, 3)
    for i, frame in enumerate(written_xtc):
        assert_array_almost_equal(frame.x, ones * i, decimal=3)


@pytest.mark.parametrize(
    'dtype', (np.int32, np.int64, np.float32, np.float64, int, float))
def test_write_xtc_dtype(tmpdir, dtype, xtc):
    fname = str(tmpdir.join("foo.xtc"))
    with XTCFile(fname, 'w') as f:
        for frame in xtc:
            x = frame.x.astype(dtype)
            box = frame.box.astype(dtype)
            f.write(
                xyz=x,
                box=box,
                step=frame.step,
                time=frame.time,
                precision=frame.prec)


@pytest.mark.parametrize('array_like', (np.array, list))
def test_write_xtc_array_like(tmpdir, array_like, xtc):
    fname = str(tmpdir.join("foo.xtc"))
    with XTCFile(fname, 'w') as f:
        for frame in xtc:
            x = array_like(frame.x)
            box = array_like(frame.box)
            f.write(
                xyz=x,
                box=box,
                step=frame.step,
                time=frame.time,
                precision=frame.prec)


def test_write_prec(tmpdir, xtc):
    outname = str(tmpdir.join('foo.xtc'))
    with XTCFile(outname, 'w') as f_out:
        assert f_out.n_atoms == 0
        frame = xtc.read()
        f_out.write(frame.x, frame.box, frame.step, frame.time, 100.0)
    xtc = XTCFile(outname)
    assert len(xtc) == 1
    frame = xtc.read()
    assert frame.prec == 100.0


def test_different_box_xtc(tmpdir, xtc):
    """test if we can write different box-sizes for different frames.
    """
    orig_box = None
    fname = str(tmpdir.join('foo.xtc'))
    with XTCFile(fname, 'w') as f_out:
        assert f_out.n_atoms == 0
        frame = xtc.read()
        f_out.write(frame.x, frame.box, frame.step, frame.time, frame.prec)
        orig_box = frame.box.copy()
        box = frame.box.copy() + 1
        f_out.write(frame.x, box, frame.step, frame.time, frame.prec)

    with XTCFile(fname) as xtc:
        assert len(xtc) == 2
        frame_1 = xtc.read()
        frame_2 = xtc.read()
        assert_array_almost_equal(frame_1.box, orig_box)
        assert_array_almost_equal(frame_1.box + 1, frame_2.box)


def test_write_different_x_xtc(tmpdir, xtc):
    outname = str(tmpdir.join('foo.xtc'))
    with XTCFile(outname, 'w') as f_out:
        assert f_out.n_atoms == 0
        frame = xtc.read()
        f_out.write(frame.x, frame.box, frame.step, frame.time, frame.prec)
        x = np.ones((xtc.n_atoms - 1, 3))
        with pytest.raises(IOError):
            f_out.write(x, frame.box, frame.step, frame.time, frame.prec)


def test_write_different_prec(tmpdir, xtc):
    outname = str(tmpdir.join('foo.xtc'))
    with XTCFile(outname, 'w') as f_out:
        assert f_out.n_atoms == 0
        frame = xtc.read()
        f_out.write(frame.x, frame.box, frame.step, frame.time, frame.prec)
        with pytest.raises(IOError):
            f_out.write(frame.x, frame.box, frame.step, frame.time, 10000.0)


@pytest.fixture
def written_trr(tmpdir, trr):
    fname = str(tmpdir.join("foo.trr"))
    with TRRFile(fname, 'w') as f:
        for frame in trr:
            natoms = frame.x.shape[0]
            f.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)
    with TRRFile(fname) as f:
        yield f


def test_written_step_trr(written_trr):
    for i, frame in enumerate(written_trr):
        assert frame.step == i


def test_written_time_trr(written_trr):
    for i, frame in enumerate(written_trr):
        assert frame.time == i * .5


def test_written_box_trr(written_trr):
    box = np.eye(3) * 20
    for frame in written_trr:
        assert_array_almost_equal(frame.box, box)


def test_written_xyx_trr(written_trr):
    ones = np.ones(30).reshape(10, 3)
    for i, frame in enumerate(written_trr):
        assert_array_almost_equal(frame.x, ones * i)


def test_written_velocities_trr(written_trr):
    ones = np.ones(30).reshape(10, 3)
    for i, frame in enumerate(written_trr):
        assert_array_almost_equal(frame.v, ones * i + 10)


def test_written_forces_trr(written_trr):
    ones = np.ones(30).reshape(10, 3)
    for i, frame in enumerate(written_trr):
        assert_array_almost_equal(frame.f, ones * i + 20)


@pytest.mark.parametrize(
    'dtype', (np.int32, np.int64, np.float32, np.float64, int, float))
def test_write_trr_dtype(tmpdir, dtype, trr):
    fname = str(tmpdir.join("foo.trr"))
    with TRRFile(fname, 'w') as fout:
        for frame in trr:
            natoms = frame.x.shape[0]
            x = frame.x.astype(dtype)
            v = frame.v.astype(dtype)
            f = frame.f.astype(dtype)
            box = frame.box.astype(dtype)
            fout.write(x, v, f, box, frame.step, frame.time, frame.lmbda,
                       natoms)


@pytest.mark.parametrize('array_like', (np.array, list))
def test_write_trr_array_like(tmpdir, array_like, trr):
    fname = str(tmpdir.join("foo.trr"))
    with TRRFile(fname, 'w') as fout:
        for frame in trr:
            natoms = frame.x.shape[0]
            x = array_like(frame.x)
            v = array_like(frame.v)
            f = array_like(frame.f)
            box = array_like(frame.box)
            fout.write(x, v, f, box, frame.step, frame.time, frame.lmbda,
                       natoms)


def test_write_different_box_trr(tmpdir, trr):
    """test if we can write different box-sizes for different frames.
    """
    orig_box = None
    fname = str(tmpdir.join('foo.trr'))
    with TRRFile(fname, 'w') as f_out:
        assert f_out.n_atoms == 0
        frame = trr.read()
        natoms = frame.x.shape[0]
        f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)
        orig_box = frame.box.copy()
        box = frame.box.copy() + 1
        f_out.write(frame.x, frame.v, frame.f, box, frame.step, frame.time,
                    frame.lmbda, natoms)

    with TRRFile(fname) as trr:
        assert len(trr) == 2
        frame_1 = trr.read()
        frame_2 = trr.read()
        assert_array_almost_equal(frame_1.box, orig_box)
        assert_array_almost_equal(frame_1.box + 1, frame_2.box)


@pytest.fixture
def trr_writer(tmpdir, trr):
    outname = str(tmpdir.join('foo.trr'))
    with TRRFile(outname, 'w') as f_out:
        assert f_out.n_atoms == 0
        frame = trr.read()
        natoms = frame.x.shape[0]
        f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)
        yield frame, f_out


def test_write_different_natoms(trr_writer):
    frame, writer = trr_writer
    natoms = frame.x.shape[0]
    with pytest.raises(IOError):
        writer.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                     frame.time, frame.lmbda, natoms - 1)


def test_write_different_x_trr(trr_writer):
    frame, writer = trr_writer
    natoms = frame.x.shape[0]
    x = np.ones((natoms - 1, 3))
    with pytest.raises(IOError):
        writer.write(x, frame.v, frame.f, frame.box, frame.step, frame.time,
                     frame.lmbda, natoms)


def test_write_different_v(trr_writer):
    frame, writer = trr_writer
    natoms = frame.x.shape[0]
    v = np.ones((natoms - 1, 3))
    with pytest.raises(IOError):
        writer.write(frame.x, v, frame.f, frame.box, frame.step, frame.time,
                     frame.lmbda, natoms)


def test_write_different_f(trr_writer):
    frame, writer = trr_writer
    natoms = frame.x.shape[0]
    f = np.ones((natoms - 1, 3))
    with pytest.raises(IOError):
        writer.write(frame.x, frame.v, f, frame.box, frame.step, frame.time,
                     frame.lmbda, natoms)
