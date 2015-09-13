from MDAnalysis.lib.formats import xtc
from MDAnalysisTests.datafiles import XTC_single_frame
from MDAnalysisTests.datafiles import XTC_multi_frame

from nose.tools import raises
from numpy.testing import assert_array_almost_equal

import numpy as np

from tempdir import run_in_tempdir


def test_n_atoms_xtc():
    f = xtc.XTCFile(XTC_single_frame)
    assert f.n_atoms == 10


def test_xtc():
    f = xtc.XTCFile(XTC_single_frame)
    assert f.n_atoms == 10
    xyz, box, step, time, prec = f.read()
    # xtc only saves with 3 decimal places precision
    assert_array_almost_equal(xyz.flat, np.arange(30), decimal=3)
    assert_array_almost_equal(box, np.eye(3) * 20, decimal=3)
    assert step == 0
    assert time == 10.0
    assert prec == 1000.0


@raises(IOError)
def test_raise_not_existing():
    xtc.XTCFile('foo')


@raises(ValueError)
def test_open_wrong_mode():
    xtc.XTCFile('foo', 'e')


@raises(RuntimeError)
def test_over_seek_xtc():
    with xtc.XTCFile(XTC_multi_frame) as f:
        f.seek(100)


@raises(RuntimeError)
@run_in_tempdir()
def test_read_write_mode_file():
    with xtc.XTCFile('foo', 'w') as f:
        f.read()


@raises(RuntimeError)
def test_read_closed():
    f = xtc.XTCFile(XTC_multi_frame)
    f.close()
    f.read()


def test_read_multi_frame_xtc():
    with xtc.XTCFile(XTC_multi_frame) as f:
        f.seek(5)
        assert f.tell() == 5
        assert len(f) == 10
        assert f.tell() == 5
        f.seek(0)
        assert f.tell() == 0
        ones = np.ones(30).reshape(10, 3)
        box_compare = np.eye(3) * 20
        for i, frame in enumerate(f):
            print i
            xyz, box, step, time, prec = frame
            assert_array_almost_equal(xyz, ones * i, decimal=3)
            assert_array_almost_equal(box, box_compare, decimal=3)
            assert step == i
            assert time == i * .5
            assert prec == 1000.0


@run_in_tempdir()
def test_write_xtc():
    with xtc.XTCFile(XTC_multi_frame) as f_in:
        with xtc.XTCFile('foo.xtc', 'w') as f_out:
            assert 0 == f_out.n_atoms
            for frame in f_in:
                f_out.write(*frame)

    with xtc.XTCFile('foo.xtc') as f:
        assert len(f) == 10
        ones = np.ones(30).reshape(10, 3)
        box_compare = np.eye(3) * 20
        for i, frame in enumerate(f):
            print i
            xyz, box, step, time, prec = frame
            assert_array_almost_equal(xyz, ones * i, decimal=3)
            assert_array_almost_equal(box, box_compare, decimal=3)
            assert step == i
            assert time == i * .5
            assert prec == 1000.0
