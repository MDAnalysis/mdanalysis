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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import numpy as np
from MDAnalysis.coordinates.base import (
    Timestep,
    SingleFrameReaderBase,
    ReaderBase
)
from numpy.testing import assert_equal, assert_allclose

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
        self._transformations = []
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
        if frame < 0:
            frame = self.n_frames + frame
        if not (0 <= frame < self.n_frames):
            raise IOError
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

    @pytest.mark.parametrize('order', ['turnip', 'abc'])
    def test_timeseries_raises_unknown_order_key(self, reader, order):
        with pytest.raises(ValueError, match="Unrecognized order key"):
            reader.timeseries(order=order)

    @pytest.mark.parametrize('order', ['faac', 'affc', 'afcc', ''])
    def test_timeseries_raises_incorrect_order_key(self, reader, order):
        with pytest.raises(ValueError, match="Repeated or missing keys"):
            reader.timeseries(order=order)

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
        (None, -1, -1),
        (10, 0, -1),
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
        (-2, 10, -2),
        (0, 0, 1),  # empty
        (10, 1, 2), # empty
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

    @pytest.mark.parametrize('sl', [
        [0, 1, 2, 3],  # ordered list of indices without duplicates
        [1, 3, 4, 2, 9],  # disordered list of indices without duplicates
        [0, 1, 1, 2, 2, 2],  # ordered list with duplicates
        [-1, -2, 3, -1, 0],  # disordered list with duplicates
        [True, ] * 10,
        [False, ] * 10,
        [True, False, ] * 5,
        slice(None, None, None),
        slice(0, 10, 1),
        slice(None, None, -1),
        slice(10, 0, -1),
        slice(2, 7, 2),
        slice(7, 2, -2),
        slice(7, 2, 1),  # empty
        slice(0, 0, 1),  # empty
    ])
    def test_getitem_len(self, sl, reader):
        traj_iterable = reader[sl]
        if not isinstance(sl, slice):
            sl = np.array(sl)
        ref = self.reference[sl]
        assert len(traj_iterable) == len(ref)

    @pytest.mark.parametrize('iter_type', (list, np.array))
    def test_getitem_len_empty(self, reader, iter_type):
        # Indexing a numpy array with an empty array tends to break.
        traj_iterable = reader[iter_type([])]
        assert len(traj_iterable) == 0

    # All the sl1 slice must be 5 frames long so that the sl2 can be a mask
    @pytest.mark.parametrize('sl1', [
        [0, 1, 2, 3, 4],
        [1, 1, 1, 1, 1],
        [True, False, ] * 5,
        slice(None, None, 2),
        slice(None, None, -2),
    ])
    @pytest.mark.parametrize('sl2', [
        [0, -1, 2],
        [-1,-1, -1],
        [True, False, True, True, False],
        np.array([True, False, True, True, False]),
        slice(None, None, None),
        slice(None, 3, None),
        slice(4, 0, -1),
        slice(None, None, -1),
        slice(None, None, 2),
    ])
    def test_double_getitem(self, sl1, sl2, reader):
        traj_iterable = reader[sl1][sl2]
        # Old versions of numpy do not behave the same when indexing with a
        # list or with an array.
        if not isinstance(sl1, slice):
            sl1 = np.asarray(sl1)
        if not isinstance(sl2, slice):
            sl2 = np.asarray(sl2)
        print(sl1, sl2, type(sl1), type(sl2))
        ref = self.reference[sl1][sl2]
        res = [ts.frame for ts in traj_iterable]
        assert_equal(res, ref)
        assert len(traj_iterable) == len(ref)

    @pytest.mark.parametrize('sl1', [
        [0, 1, 2, 3, 4],
        [1, 1, 1, 1, 1],
        [True, False, ] * 5,
        slice(None, None, 2),
        slice(None, None, -2),
        slice(None, None, None),
    ])
    @pytest.mark.parametrize('idx2', [0, 2, 4, -1, -2, -4])
    def test_double_getitem_int(self, sl1, idx2, reader):
        ts = reader[sl1][idx2]
        # Old versions of numpy do not behave the same when indexing with a
        # list or with an array.
        if not isinstance(sl1, slice):
            sl1 = np.asarray(sl1)
        ref = self.reference[sl1][idx2]
        assert ts.frame == ref

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

    @pytest.mark.parametrize('sl1', [
        [0, 1, 2, 3, 4],
        [1, 1, 1, 1, 1],
        [True, False, ] * 5,
        slice(None, None, 2),
        slice(None, None, -2),
    ])
    @pytest.mark.parametrize('idx2', [5, -6])
    def test_getitem_IE(self, sl1, idx2, reader):
        partial_reader = reader[sl1]
        with pytest.raises(IndexError):
            partial_reader[idx2]


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

    def test_iter_rewind(self, reader):
        # Issue #3423
        # positions should be zero to start with
        assert_allclose(reader.ts.positions, np.zeros((10, 3)))

        # modify positions in place
        reader.ts.positions = np.ones((10, 3))
        assert_allclose(reader.ts.positions, np.ones((10, 3)))

        # iterate the reader
        for ts in reader:
            assert_allclose(ts.positions, np.zeros((10, 3)))

        assert_allclose(reader.ts.positions, np.zeros((10, 3)))
