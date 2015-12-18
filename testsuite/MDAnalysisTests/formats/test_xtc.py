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
from MDAnalysis.lib.formats.xdrlib import XTCFile
from MDAnalysisTests.datafiles import XTC_single_frame
from MDAnalysisTests.datafiles import XTC_multi_frame

from nose.tools import raises
from numpy.testing import assert_array_almost_equal

import numpy as np
from tempdir import run_in_tempdir
from .base import XDRFormatBaseTest


class Test_XTCFile(XDRFormatBaseTest):
    xdrfile = XTCFile
    single_frame = XTC_single_frame
    multi_frame = XTC_multi_frame
    offsets = np.array([0, 104, 208, 312, 416, 520, 624, 728, 832, 936])


def test_xtc():
    f = XTCFile(XTC_single_frame)
    assert f.n_atoms == 10
    xyz, box, step, time, prec = f.read()
    # xtc only saves with 3 decimal places precision
    assert_array_almost_equal(xyz.flat, np.arange(30), decimal=3)
    assert_array_almost_equal(box, np.eye(3) * 20, decimal=3)
    assert step == 0
    assert time == 10.0
    assert prec == 1000.0


def test_read_multi_frame_xtc():
    with XTCFile(XTC_multi_frame) as f:
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
    with XTCFile(XTC_multi_frame) as f_in:
        with XTCFile('foo.xtc', 'w') as f_out:
            assert 0 == f_out.n_atoms
            for frame in f_in:
                f_out.write(*frame)

    with XTCFile('foo.xtc') as f:
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


@raises(ValueError)
@run_in_tempdir()
def test_write_different_box():
    with XTCFile(XTC_multi_frame) as f_in:
        with XTCFile('foo.xtc', 'w') as f_out:
            assert 0 == f_out.n_atoms
            frame = f_in.read()
            f_out.write(frame.x, frame.box, frame.step,
                        frame.time, frame.prec)
            box = frame.box + 1
            f_out.write(frame.x, box, frame.step,
                        frame.time, frame.prec)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_x():
    with XTCFile(XTC_multi_frame) as f_in:
        with XTCFile('foo.xtc', 'w') as f_out:
            assert 0 == f_out.n_atoms
            frame = f_in.read()
            f_out.write(frame.x, frame.box, frame.step,
                        frame.time, frame.prec)
            x = np.ones((f_in.n_atoms - 1, 3))
            f_out.write(x, frame.box, frame.step,
                        frame.time, frame.prec)


@raises(ValueError)
@run_in_tempdir()
def test_write_different_prec():
    with XTCFile(XTC_multi_frame) as f_in:
        with XTCFile('foo.xtc', 'w') as f_out:
            assert 0 == f_out.n_atoms
            frame = f_in.read()
            f_out.write(frame.x, frame.box, frame.step,
                        frame.time, frame.prec)
            f_out.write(frame.x, frame.box, frame.step,
                        frame.time, frame.prec - 900)
