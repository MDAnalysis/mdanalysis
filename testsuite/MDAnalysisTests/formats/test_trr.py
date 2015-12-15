from MDAnalysis.lib.formats.xdrlib import TRRFile
from MDAnalysisTests.datafiles import TRR_single_frame
from MDAnalysisTests.datafiles import TRR_multi_frame

from nose.tools import raises
from numpy.testing import assert_array_almost_equal, assert_almost_equal

import numpy as np
from tempdir import run_in_tempdir
from .base import XDRFormatBaseTest


class Test_TRRFile(XDRFormatBaseTest):
    xdrfile = TRRFile
    single_frame = TRR_single_frame
    multi_frame = TRR_multi_frame
    offsets = np.array([0, 480, 960, 1440, 1920, 2400, 2880,
                        3360, 3840, 4320])


def test_trr():
    f = TRRFile(TRR_single_frame)
    assert f.n_atoms == 10
    frame = f.read()
    # trr only saves with 3 decimal places precision
    assert_array_almost_equal(frame.x.flat, np.arange(30))
    assert_array_almost_equal(frame.v.flat, np.arange(30))
    assert_array_almost_equal(frame.f.flat, np.arange(30))
    assert_array_almost_equal(frame.box, np.eye(3) * 20)
    assert frame.step == 10
    assert frame.time == 5.0
    assert frame.lmbda == .5


def test_read_multi_frame_trr():
    with TRRFile(TRR_multi_frame) as f:
        f.seek(5)
        assert f.tell() == 5
        assert len(f) == 10
        assert f.tell() == 5
        f.seek(0)
        assert f.tell() == 0
        ones = np.ones(30).reshape(10, 3)
        box_compare = np.eye(3) * 20
        for i, frame in enumerate(f):
            assert_array_almost_equal(frame.x, ones * i)
            assert_array_almost_equal(frame.v, ones * i + 10)
            assert_array_almost_equal(frame.f, ones * i + 20)
            assert_array_almost_equal(frame.box, box_compare)
            assert frame.step == i
            assert frame.time == i * .5
            assert_almost_equal(frame.lmbda, .01 * i)


@run_in_tempdir()
def test_write_trr():
    with TRRFile(TRR_multi_frame) as f_in:
        with TRRFile('foo.trr', 'w') as f_out:
            assert 0 == f_out.n_atoms
            for frame in f_in:
                natoms = frame.x.shape[0]
                f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                            frame.time, frame.lmbda, natoms)

    with TRRFile('foo.trr') as f:
        assert len(f) == 10
        ones = np.ones(30).reshape(10, 3)
        box_compare = np.eye(3) * 20
        for i, frame in enumerate(f):
            assert_array_almost_equal(frame.x, ones * i)
            assert_array_almost_equal(frame.v, ones * i + 10)
            assert_array_almost_equal(frame.f, ones * i + 20)
            assert_array_almost_equal(frame.box, box_compare)
            assert frame.step == i
            assert frame.time == i * .5
            assert_almost_equal(frame.lmbda, .01 * i)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_natoms():
    with TRRFile(TRR_multi_frame) as f_in:
        with TRRFile('foo.trr', 'w') as f_out:
            assert 0 == f_out.n_atoms
            frame = f_in.read()
            natoms = frame.x.shape[0]
            f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                        frame.time, frame.lmbda, natoms)
            f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                        frame.time, frame.lmbda, natoms - 1)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_box():
    with TRRFile(TRR_multi_frame) as f_in:
        with TRRFile('foo.trr', 'w') as f_out:
            assert 0 == f_out.n_atoms
            frame = f_in.read()
            natoms = frame.x.shape[0]
            f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                        frame.time, frame.lmbda, natoms)
            box = frame.box + 1
            f_out.write(frame.x, frame.v, frame.f, box, frame.step,
                        frame.time, frame.lmbda, natoms)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_x():
    with TRRFile(TRR_multi_frame) as f_in:
        with TRRFile('foo.trr', 'w') as f_out:
            assert 0 == f_out.n_atoms
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
    with TRRFile(TRR_multi_frame) as f_in:
        with TRRFile('foo.trr', 'w') as f_out:
            assert 0 == f_out.n_atoms
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
    with TRRFile(TRR_multi_frame) as f_in:
        with TRRFile('foo.trr', 'w') as f_out:
            assert 0 == f_out.n_atoms
            frame = f_in.read()
            natoms = frame.x.shape[0]
            f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                        frame.time, frame.lmbda, natoms)
            f = np.ones((natoms - 1, 3))
            f_out.write(frame.x, frame.v, f, frame.box, frame.step,
                        frame.time, frame.lmbda, natoms)
