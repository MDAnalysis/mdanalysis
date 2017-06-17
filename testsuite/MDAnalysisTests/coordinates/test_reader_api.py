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
from numpy.testing import assert_equal, assert_raises

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
    """Basic API readers"""

    def setUp(self):
        self.reader = self.readerclass('test.txt')
        self.ts = self.reader.ts

    def test_required_attributes(self):
        """Test that Reader has the required attributes"""
        for attr in ['filename', 'n_atoms', 'n_frames', 'ts',
                     'units', 'format']:
            assert_equal(hasattr(self.reader, attr), True,
                         "Missing attr: {0}".format(attr))

    def test_iter(self):
        l = [ts for ts in self.reader]

        assert_equal(len(l), self.n_frames)

    def test_close(self):
        sfr = self.readerclass('text.txt')

        ret = sfr.close()
        # Check that method works?
        assert_equal(ret, None)

    def test_rewind(self):
        ret = self.reader.rewind()

        assert_equal(ret, None)
        assert_equal(self.reader.ts.frame, 0)

    def test_context(self):
        with self.readerclass('text.txt') as sfr:
            l = sfr.ts.frame

        assert_equal(l, 0)

    def test_len(self):
        l = len(self.reader)

        assert_equal(l, self.n_frames)

    def test_raises_StopIteration(self):
        self.reader[-1]

        assert_raises(StopIteration, next, self.reader)


class _Multi(_TestReader):
    n_frames = 10
    n_atoms = 10
    readerclass = AmazingMultiFrameReader
    reference = np.arange(10)


class TestMultiFrameReader(_Multi):
    def _check_slice(self, start, stop, step):
        """Compare the slice applied to trajectory, to slice of list"""
        res = [ts.frame for ts in self.reader[start:stop:step]]
        ref = self.reference[start:stop:step]

        assert_equal(res, ref)

    def test_slices(self):
        for start, stop, step in [
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
            (5, 1, 1), # Stop less than start
            (1, 5, -1), # Stop less than start
            (-100, None, None),
            (100, None, None), # Outside of range of trajectory
            (-2, 10, -2)
        ]:
            yield self._check_slice, start, stop, step

    def test_slice_VE_1(self):
        def sl():
            return list(self.reader[::0])

        assert_raises(ValueError, sl)

    def test_slice_TE_1(self):
        def sl():
            return list(self.reader[1.2:2.5:0.1])

        assert_raises(TypeError, sl)

    def _check_getitem(self, sl):
        res = [ts.frame for ts in self.reader[sl]]

        sl = np.asarray(sl)
        ref = self.reference[sl]

        assert_equal(res, ref)

    def test_getitem_list_ints(self):
        for sl in (
                [0, 1, 4, 5],
                np.array([0, 1, 4, 5]),
                [5, 1, 6, 2, 7, 3, 8],
                np.array([5, 1, 6, 2, 7, 3, 8]),
                [0, 1, 1, 1, 0, 0, 2, 3, 4],
                np.array([0, 1, 1, 1, 0, 0, 2, 3, 4]),
        ):
            yield self._check_getitem, sl

    def test_list_TE(self):
        def sl():
            return list(self.reader[[0, 'a', 5, 6]])

        assert_raises(TypeError, sl)

    def test_array_TE(self):
        def sl():
            return list(self.reader[np.array([1.2, 3.4, 5.6])])

        assert_raises(TypeError, sl)

    def test_bool_slice(self):
        t = True
        f = False
        for sl in (
                [t, f, t, f, t, f, t, f, t, f],
                [t, t, f, f, t, t, f, t, f, t],
                [t, t, t, t, t, t, t, t, t, t],
                [f, f, f, f, f, f, f, f, f, f],
        ):
            yield self._check_getitem, sl
            yield self._check_getitem, np.array(sl, dtype=np.bool)


class _Single(_TestReader):
    n_frames = 1
    n_atoms = 10
    readerclass = AmazingReader


class TestSingleFrameReader(_Single):
    def test_next(self):
        assert_raises(StopIteration, self.reader.next)

    # Getitem tests
    # only 0 & -1 should work
    # others should get IndexError
    def _check_get_results(self, l):
        assert_equal(len(l), 1)
        assert_equal(self.ts in l, True)

    def test_getitem(self):
        fr = [self.reader[0]]

        self._check_get_results(fr)

    def test_getitem_2(self):
        fr = [self.reader[-1]]

        self._check_get_results(fr)

    def test_getitem_IE(self):
        assert_raises(IndexError, self.reader.__getitem__, 1)

    def test_getitem_IE_2(self):
        assert_raises(IndexError, self.reader.__getitem__, -2)

    # Slicing should still work!
    def test_slice_1(self):
        l = list(self.reader[::])
        self._check_get_results(l)

    def test_slice_2(self):
        l = list(self.reader[::-1])
        self._check_get_results(l)

    def test_reopen(self):
        self.reader._reopen()
        assert_equal(self.ts.frame, 0)

    def test_rewind(self):
        self.reader.rewind()
        assert_equal(self.ts.frame, 0)

    def test_read_frame(self):
        assert_raises(IndexError, self.reader._read_frame, 1)
