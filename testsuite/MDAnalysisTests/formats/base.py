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
from nose.tools import raises
from numpy.testing import assert_equal, assert_array_equal

from tempdir import run_in_tempdir
import numpy as np


class XDRFormatBaseTest():

    def test_n_atoms(self):
        f = self.xdrfile(self.single_frame)
        assert f.n_atoms == 10

    @raises(IOError)
    def test_raise_not_existing(self):
        self.xdrfile('foo')

    @raises(ValueError)
    def test_open_wrong_mode(self):
        self.xdrfile('foo', 'e')

    @raises(RuntimeError)
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

    @raises(RuntimeError)
    def test_set_offsets(self):
        f = self.xdrfile(self.multi_frame)
        f.set_offsets(np.arange(len(self.offsets)))
        assert_array_equal(f.offsets, np.arange(len(self.offsets)))
        f.seek(6)
        f.read()
