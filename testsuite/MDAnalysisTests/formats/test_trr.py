from MDAnalysis.lib.formats import trr
from MDAnalysisTests.datafiles import TRR_single_frame
from MDAnalysisTests.datafiles import TRR_multi_frame

from nose.tools import raises
from numpy.testing import assert_array_almost_equal, assert_almost_equal

import numpy as np

from tempdir import run_in_tempdir


def test_n_atoms_trr():
    f = trr.TRRFile(TRR_single_frame)
    assert f.n_atoms == 10


def test_trr():
    f = trr.TRRFile(TRR_single_frame)
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


@raises(IOError)
def test_raise_not_existing():
    trr.TRRFile('foo')


@raises(ValueError)
def test_open_wrong_mode():
    trr.TRRFile('foo', 'e')


@raises(RuntimeError)
def test_over_seek_trr():
    with trr.TRRFile(TRR_multi_frame) as f:
        f.seek(100)


@raises(RuntimeError)
@run_in_tempdir()
def test_read_write_mode_file():
    with trr.TRRFile('foo', 'w') as f:
        f.read()


@raises(RuntimeError)
def test_read_closed():
    f = trr.TRRFile(TRR_multi_frame)
    f.close()
    f.read()


def test_read_multi_frame_trr():
    with trr.TRRFile(TRR_multi_frame) as f:
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
    with trr.TRRFile(TRR_multi_frame) as f_in:
        with trr.TRRFile('foo.trr', 'w') as f_out:
            assert 0 == f_out.n_atoms
            for frame in f_in:
                natoms = frame.x.shape[0]
                f_out.write(frame.x, frame.v, frame.f, frame.box, frame.step,
                            frame.time, frame.lmbda, natoms)

    with trr.TRRFile('foo.trr') as f:
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
