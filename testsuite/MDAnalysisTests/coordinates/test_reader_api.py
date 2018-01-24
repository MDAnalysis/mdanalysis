# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

import numpy as np
from MDAnalysis.coordinates.base import (
    Timestep,
    SingleFrameReaderBase,
    ReaderBase
)
from numpy.testing import assert_equal

import pytest
"""
Isolate the API definitions of Readers independent of implementations
"""


class AmazingMultiFrameReader(ReaderBase):
    format = 'AmazingMulti'

    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.n_frames = 10
        self.n_atoms = 10
        self._auxs = {}
        # ts isn't a real timestep, but just an integer
        # whose value represents the frame number (0 based)
        self.ts = Timestep(self.n_atoms)
        self.ts.frame = -1
        self._read_next_timestep()

    def _read_next_timestep(self):
        self.ts.frame += 1
        if (self.ts.frame + 1) > self.n_frames:
            raise IOError
        else:
            return self.ts

    def _read_frame(self, frame):
        self.ts.frame = frame

        return self.ts

    def _reopen(self):
        self.ts.frame = -1


class AmazingReader(SingleFrameReaderBase):
    format = 'Amazing'

    # have to hack this in to get the base class to "work"
    def _read_first_frame(self):
        self.n_atoms = 10
        self.ts = Timestep(self.n_atoms)
        self.ts.frame = 0


class _TestReader(object):
    __test__ = False
    """Basic API readers"""

    @pytest.fixture()
    def reader(self):
        return self.readerclass('test.txt')

    @pytest.fixture()
    def ts(self, reader):
        return reader.ts

    def test_required_attributes(self, reader):
        """Test that Reader has the required attributes"""
        for attr in ['filename', 'n_atoms', 'n_frames', 'ts',
                     'units', 'format']:
            assert_equal(hasattr(reader, attr), True,
                         "Missing attr: {0}".format(attr))

    def test_iter(self, reader):
        l = [ts for ts in reader]

        assert_equal(len(l), self.n_frames)

    def test_close(self):
        sfr = self.readerclass('text.txt')

        ret = sfr.close()
        # Check that method works?
        assert_equal(ret, None)

    def test_rewind(self, reader):
        ret = reader.rewind()

        assert_equal(ret, None)
        assert_equal(reader.ts.frame, 0)

    def test_context(self):
        with self.readerclass('text.txt') as sfr:
            l = sfr.ts.frame

        assert_equal(l, 0)

    def test_len(self, reader):
        l = len(reader)

        assert_equal(l, self.n_frames)

    def test_raises_StopIteration(self, reader):
        reader[-1]
        with pytest.raises(StopIteration):
            next(reader)


class _Multi(_TestReader):
    n_frames = 10
    n_atoms = 10
    readerclass = AmazingMultiFrameReader
    reference = np.arange(10)


class TestMultiFrameReader(_Multi):

    __test__ = True

    @pytest.mark.parametrize('start, stop, step', [
        (None, None, None),  # blank slice
        (None, 5, None),  # set end point
        (2, None, None),  # set start point
        (2, 5, None),  # start & end
        (None, None, 2),  # set skip
        (None, None, -1),  # backwards skip
        (0, 10, 1),
        (0, 10, 2),
        (None, 20, None),  # end beyond real end
        (None, 20, 2),  # with skip
        (0, 5, 2),
        (5, None, -1),
        (None, 5, -1),
        (100, 10, 1),
        (-10, None, 1),
        (100, None, -1),  # beyond real end
        (100, 5, -20),
        (5, 1, 1),  # Stop less than start
        (1, 5, -1),  # Stop less than start
        (-100, None, None),
        (100, None, None),  # Outside of range of trajectory
        (-2, 10, -2)
    ])
    def test_slice(self, start, stop, step, reader):
        """Compare the slice applied to trajectory, to slice of list"""
        res = [ts.frame for ts in reader[start:stop:step]]
        ref = self.reference[start:stop:step]

        assert_equal(res, ref)

    def test_slice_VE_1(self, reader):
        def sl():
            return list(reader[::0])

        with pytest.raises(ValueError):
            sl()

    def test_slice_TE_1(self, reader):
        def sl():
            return list(reader[1.2:2.5:0.1])

        with pytest.raises(TypeError):
            sl()


    @pytest.mark.parametrize('slice_cls', [list, np.array])
    @pytest.mark.parametrize('sl', [
        [0, 1, 4, 5],
        [5, 1, 6, 2, 7, 3, 8],
        [0, 1, 1, 1, 0, 0, 2, 3, 4],
        [True, False, True, False, True, False, True, False, True, False],
        [True, True, False, False, True, True, False, True, False, True],
        [True, True, True, True, True, True, True, True, True, True],
        [False, False, False, False, False, False, False, False, False, False],
    ])
    def test_getitem(self, slice_cls, sl, reader):
        sl = slice_cls(sl)
        res = [ts.frame for ts in reader[sl]]

        sl = np.asarray(sl)
        ref = self.reference[sl]

        assert_equal(res, ref)

    def test_list_TE(self, reader):
        def sl():
            return list(reader[[0, 'a', 5, 6]])

        with pytest.raises(TypeError):
            sl()


    def test_array_TE(self, reader):
        def sl():
            return list(reader[np.array([1.2, 3.4, 5.6])])

        with pytest.raises(TypeError):
            sl()


class _Single(_TestReader):
    n_frames = 1
    n_atoms = 10
    readerclass = AmazingReader


class TestSingleFrameReader(_Single):
    __test__ = True
    def test_next(self, reader):
        with pytest.raises(StopIteration):
            reader.next()

    # Getitem tests
    # only 0 & -1 should work
    # others should get IndexError
    def _check_get_results(self, l, ts):
        assert_equal(len(l), 1)
        assert_equal(ts in l, True)

    def test_getitem(self, reader, ts):
        fr = [reader[0]]

        self._check_get_results(fr, ts)

    def test_getitem_2(self, reader, ts):
        fr = [reader[-1]]

        self._check_get_results(fr, ts)

    def test_getitem_IE(self, reader):
        with pytest.raises(IndexError):
            reader.__getitem__(1)

    def test_getitem_IE_2(self, reader):
        with pytest.raises(IndexError):
            reader.__getitem__(-2)

    # Slicing should still work!
    def test_slice_1(self, reader, ts):
        l = list(reader[::])
        self._check_get_results(l, ts)

    def test_slice_2(self, reader, ts):
        l = list(reader[::-1])
        self._check_get_results(l, ts)

    def test_reopen(self, reader, ts):
        reader._reopen()
        assert_equal(ts.frame, 0)

    def test_rewind(self, reader, ts):
        reader.rewind()
        assert_equal(ts.frame, 0)

    def test_read_frame(self, reader):
        with pytest.raises(IndexError):
            reader._read_frame(1)
