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
from __future__ import print_function

from nose.tools import raises
from numpy.testing import assert_equal, assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal

from MDAnalysis.lib.formats.libmdaxdr import TRRFile
from MDAnalysisTests.datafiles import TRR_single_frame
from MDAnalysisTests.datafiles import TRR_multi_frame

from MDAnalysis.lib.formats.libmdaxdr import XTCFile
from MDAnalysisTests.datafiles import XTC_single_frame
from MDAnalysisTests.datafiles import XTC_multi_frame

from tempdir import run_in_tempdir
import numpy as np


class XDRFormatBaseTest(object):

    def test_n_atoms(self):
        f = self.xdrfile(self.single_frame)
        assert_equal(f.n_atoms, 10)

    @raises(IOError)
    def test_raise_not_existing(self):
        self.xdrfile('foo')

    @raises(ValueError)
    def test_open_wrong_mode(self):
        self.xdrfile('foo', 'e')

    @raises(IOError)
    def test_over_seek(self):
        with self.xdrfile(self.multi_frame) as f:
            f.seek(100)

    @raises(RuntimeError)
    @run_in_tempdir()
    def test_read_write_mode_file(self):
        with self.xdrfile('foo', 'w') as f:
            f.read()

    @raises(RuntimeError)
    def test_read_closed(self):
        f = self.xdrfile(self.multi_frame)
        f.close()
        f.read()

    def test_iter(self):
        with self.xdrfile(self.multi_frame) as f:
            for frame in f:
                pass

    def test_tell(self):
        f = self.xdrfile(self.multi_frame)
        assert_equal(f.tell(), 0)
        for i, frame in enumerate(f):
            assert_equal(f.tell(), i + 1)

    def test_seek(self):
        f = self.xdrfile(self.multi_frame)
        f.seek(4)
        assert_equal(f.tell(), 4)

    def test_offset(self):
        f = self.xdrfile(self.multi_frame)
        assert_array_equal(f.offsets, self.offsets)

    @raises(IOError)
    def test_set_offsets(self):
        f = self.xdrfile(self.multi_frame)
        f.set_offsets(np.arange(len(self.offsets)))
        assert_array_equal(f.offsets, np.arange(len(self.offsets)))
        f.seek(6)
        f.read()

    @raises(OverflowError)
    def test_seek_overflow(self):
        f = self.xdrfile(self.multi_frame)
        f._bytes_seek(2**65)

    def test_seek_tell_all(self):
        f = self.xdrfile(self.multi_frame)
        for frame, offset in enumerate(self.offsets):
            f.seek(frame)
            assert_equal(f._bytes_tell(), offset)

    def test_seek_tell_all_bytes(self):
        f = self.xdrfile(self.multi_frame)
        for offset in self.offsets:
            f._bytes_seek(offset)
            assert_equal(f._bytes_tell(), offset)

    def test_seek_tell_largefile(self):
        # Seeking/telling can be done on offsets larger than the file.
        # Filesize won't change unless a write is done at the offset.
        # We attempt to go just beyond a 4GB limit.
        big_offset = 2**32+1000
        f = self.xdrfile(self.multi_frame)
        f._bytes_seek(big_offset)
        assert_equal(f._bytes_tell(), big_offset)
        # It should be safe to leave the read head so far away
        # because the file is open read-only.


class Test_XTCFile(XDRFormatBaseTest):
    xdrfile = XTCFile
    single_frame = XTC_single_frame
    multi_frame = XTC_multi_frame
    offsets = np.array([0, 104, 208, 312, 416, 520, 624, 728, 832, 936])


def test_xtc():
    f = XTCFile(XTC_single_frame)
    assert_equal(f.n_atoms, 10)
    xyz, box, step, time, prec = f.read()
    # xtc only saves with 3 decimal places precision
    assert_array_almost_equal(xyz.flat, np.arange(30), decimal=3)
    assert_array_almost_equal(box, np.eye(3) * 20, decimal=3)
    assert_equal(step, 0)
    assert_equal(time, 10.0)
    assert_equal(prec, 1000.0)


def test_read_multi_frame_xtc():
    with XTCFile(XTC_multi_frame) as f:
        f.seek(5)
        assert_equal(f.tell(), 5)
        assert_equal(len(f), 10)
        assert_equal(f.tell(), 5)
        f.seek(0)
        assert_equal(f.tell(), 0)
        ones = np.ones(30).reshape(10, 3)
        box_compare = np.eye(3) * 20
        for i, frame in enumerate(f):
            xyz, box, step, time, prec = frame
            assert_array_almost_equal(xyz, ones * i, decimal=3)
            assert_array_almost_equal(box, box_compare, decimal=3)
            assert_equal(step, i)
            assert_equal(time, i * 0.5)
            assert_equal(prec, 1000.0)


@run_in_tempdir()
def test_write_xtc():
    with XTCFile(XTC_multi_frame) as f_in, XTCFile('foo.xtc', 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        for frame in f_in:
            f_out.write(*frame)

    with XTCFile('foo.xtc') as f:
        assert_equal(len(f), 10)
        ones = np.ones(30).reshape(10, 3)
        box_compare = np.eye(3) * 20
        for i, frame in enumerate(f):
            xyz, box, step, time, prec = frame
            assert_array_almost_equal(xyz, ones * i, decimal=3)
            assert_array_almost_equal(box, box_compare, decimal=3)
            assert_equal(step, i)
            assert_equal(time, i * 0.5)
            assert_equal(prec, 1000.0)


@run_in_tempdir()
def test_different_box_xtc():
    """test if we can write different box-sizes for different frames.
    """
    fname = 'foo.xtc'
    orig_box = None
    with XTCFile(XTC_multi_frame) as f_in, XTCFile(fname, 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        frame = f_in.read()
        f_out.write(frame.x, frame.box, frame.step,
                    frame.time, frame.prec)
        orig_box = frame.box.copy()
        box = frame.box.copy() + 1
        f_out.write(frame.x, box, frame.step,
                    frame.time, frame.prec)

    with XTCFile(fname) as xtc:
        assert_equal(len(xtc), 2)
        frame_1 = xtc.read()
        frame_2 = xtc.read()

        assert_array_almost_equal(frame_1.box, orig_box)
        assert_array_almost_equal(frame_1.box + 1, frame_2.box)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_x_xtc():
    with XTCFile(XTC_multi_frame) as f_in, XTCFile('foo.xtc', 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        frame = f_in.read()
        f_out.write(frame.x, frame.box, frame.step,
                    frame.time, frame.prec)
        x = np.ones((f_in.n_atoms - 1, 3))
        f_out.write(x, frame.box, frame.step,
                    frame.time, frame.prec)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_prec():
    mf_xtc = XTCFile(XTC_multi_frame)
    with XTCFile('foo.xtc', 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        frame = mf_xtc.read()
        f_out.write(frame.x, frame.box, frame.step,
                    frame.time, frame.prec)
        f_out.write(frame.x, frame.box, frame.step,
                    frame.time, 10000.0)


@run_in_tempdir()
def test_write_prec():
    mf_xtc = XTCFile(XTC_multi_frame)
    with XTCFile('foo.xtc', 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        frame = mf_xtc.read()
        f_out.write(frame.x, frame.box, frame.step,
                    frame.time, 100.0)

    xtc = XTCFile('foo.xtc')
    assert_equal(len(xtc), 1)
    frame = xtc.read()
    assert_equal(frame.prec, 100.0)


class Test_TRRFile(XDRFormatBaseTest):
    xdrfile = TRRFile
    single_frame = TRR_single_frame
    multi_frame = TRR_multi_frame
    offsets = np.array([0, 480, 960, 1440, 1920, 2400, 2880,
                        3360, 3840, 4320])


def test_trr():
    f = TRRFile(TRR_single_frame)
    assert_equal(f.n_atoms, 10)
    frame = f.read()
    # trr only saves with 3 decimal places precision
    assert_array_almost_equal(frame.x.flat, np.arange(30))
    assert_array_almost_equal(frame.v.flat, np.arange(30))
    assert_array_almost_equal(frame.f.flat, np.arange(30))
    assert_array_almost_equal(frame.box, np.eye(3) * 20)
    assert_equal(frame.step, 10)
    assert_equal(frame.time, 5.0)
    assert_equal(frame.lmbda, 0.5)


def test_read_multi_frame_trr():
    with TRRFile(TRR_multi_frame) as f:
        f.seek(5)
        assert_equal(f.tell(), 5)
        assert_equal(len(f), 10)
        assert_equal(f.tell(), 5)
        f.seek(0)
        assert_equal(f.tell(), 0)
        ones = np.ones(30).reshape(10, 3)
        box_compare = np.eye(3) * 20
        for i, frame in enumerate(f):
            assert_array_almost_equal(frame.x, ones * i)
            assert_array_almost_equal(frame.v, ones * i + 10)
            assert_array_almost_equal(frame.f, ones * i + 20)
            assert_array_almost_equal(frame.box, box_compare)
            assert_equal(frame.step, i)
            assert_equal(frame.time, i * 0.5)
            assert_almost_equal(frame.lmbda, .01 * i)


@run_in_tempdir()
def test_write_trr():
    with TRRFile(TRR_multi_frame) as f_in, TRRFile('foo.trr', 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        for frame in f_in:
            natoms = frame.x.shape[0]
            f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                        frame.time, frame.lmbda, natoms)

    with TRRFile('foo.trr') as f:
        assert_equal(len(f), 10)
        ones = np.ones(30).reshape(10, 3)
        box_compare = np.eye(3) * 20
        for i, frame in enumerate(f):
            assert_array_almost_equal(frame.x, ones * i)
            assert_array_almost_equal(frame.v, ones * i + 10)
            assert_array_almost_equal(frame.f, ones * i + 20)
            assert_array_almost_equal(frame.box, box_compare)
            assert_equal(frame.step, i)
            assert_equal(frame.time, i * 0.5)
            assert_almost_equal(frame.lmbda, .01 * i)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_natoms():
    with TRRFile(TRR_multi_frame) as f_in, TRRFile('foo.trr', 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        frame = f_in.read()
        natoms = frame.x.shape[0]
        f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)
        f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms - 1)


@run_in_tempdir()
def test_write_different_box_trr():
    """test if we can write different box-sizes for different frames.
    """
    fname = 'foo.trr'
    orig_box = None
    with TRRFile(TRR_multi_frame) as f_in, TRRFile(fname, 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        frame = f_in.read()
        natoms = frame.x.shape[0]
        f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)
        orig_box = frame.box.copy()
        box = frame.box.copy() + 1
        f_out.write(frame.x, frame.v, frame.f, box, frame.step,
                    frame.time, frame.lmbda, natoms)

    with TRRFile(fname) as trr:
        assert_equal(len(trr), 2)
        frame_1 = trr.read()
        frame_2 = trr.read()

        assert_array_almost_equal(frame_1.box, orig_box)
        assert_array_almost_equal(frame_1.box + 1, frame_2.box)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_x_trr():
    with TRRFile(TRR_multi_frame) as f_in, TRRFile('foo.trr', 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        frame = f_in.read()
        natoms = frame.x.shape[0]
        f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)
        x = np.ones((natoms - 1, 3))
        f_out.write(x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_v():
    with TRRFile(TRR_multi_frame) as f_in, TRRFile('foo.trr', 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        frame = f_in.read()
        natoms = frame.x.shape[0]
        f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)
        v = np.ones((natoms - 1, 3))
        f_out.write(frame.x, v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_f():
    with TRRFile(TRR_multi_frame) as f_in, TRRFile('foo.trr', 'w') as f_out:
        assert_equal(f_out.n_atoms, 0)
        frame = f_in.read()
        natoms = frame.x.shape[0]
        f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)
        f = np.ones((natoms - 1, 3))
        f_out.write(frame.x, frame.v, f, frame.box, frame.step,
                    frame.time, frame.lmbda, natoms)
