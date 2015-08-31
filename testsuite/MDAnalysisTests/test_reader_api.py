# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from MDAnalysis.coordinates.base import Timestep, SingleFrameReader, Reader

from numpy.testing import *

"""
Isolate the API definitions of Readers independent of implementations
"""

class AmazingMultiFrameReader(Reader):
    format = 'AmazingMulti'

    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.n_frames = 10
        self.n_atoms = 10
        # ts isn't a real timestep, but just an integer
        # whose value represents the frame number (0 based)
        self.ts = Timestep(self.n_atoms)
        self.ts.frame = -1
        self._read_next_timestep()

    def __iter__(self):
        self._reopen()
        while True:
            try:
                yield self._read_next_timestep()
            except IOError:
                break

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


class AmazingReader(SingleFrameReader):
    format = 'Amazing'
    # have to hack this in to get the base class to "work"
    def _read_first_frame(self):
        self.n_atoms = 10
        self.ts = Timestep(self.n_atoms)
        self.ts.frame = 0


class _TestReader(TestCase):
    """Basic API for readers"""
    def setUp(self):
        self.reader = self.readerclass('test.txt')
        self.ts = self.reader.ts
    
    def test_required_attributes(self):
        """Test that Reader has the required attributes"""
        for attr in ['filename', 'n_atoms', 'n_frames', 'ts',
                     'units', 'format']:
            assert_equal(hasattr(self.reader, attr), True, "Missing attr: {0}".format(attr))
        
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


class _Multi(object):
    n_frames = 10
    n_atoms = 10
    readerclass = AmazingMultiFrameReader
    reference = [i for i in range(10)]
   

class TestMultiFrameReader(_Multi, _TestReader):
    def _check_slice(self, sl):
        """Compare the slice applied to trajectory, to slice of list"""
        res = [ts.frame for ts in self.reader[sl]]
        ref = self.reference[sl]

        assert_equal(res, ref)
    
    def test_slice_1(self):
        sl = slice(0, 10, 1)
        self._check_slice(sl)

    def test_slice_2(self):
        sl = slice(0, 10, 2)
        self._check_slice(sl)

    def test_slice_3(self):
        """Upper bound below traj length"""
        sl = slice(0, 5, 2)
        self._check_slice(sl)

    def test_slice_4(self):
        """Upper bound above traj length"""
        sl = slice(0, 20, 2)
        self._check_slice(sl)

    def test_slice_5(self):
        """Reverse order"""
        sl = slice(0, 10, -1)
        self._check_slice(sl)

    def test_slice_IE_1(self):
        """Stop less than start"""
        def sl():
            return list(self.reader[5:1:1])
        assert_raises(IndexError, sl)

    def test_slice_IE_2(self):
        """Outside of range of trajectory"""
        def sl():
            return list(self.reader[100:])
        assert_raises(IndexError, sl)

    def test_slice_IE_3(self):
        def sl():
            return list(self.reader[-100:])
        assert_raises(IndexError, sl)

    def test_slice_VE_1(self):
        def sl():
            return list(self.reader[::0])
        assert_raises(ValueError, sl)

    def test_slice_TE_1(self):
        def sl():
            return list(self.reader[1.2:2.5:0.1])
        assert_raises(TypeError, sl)


class _Single(TestCase):
    n_frames = 1
    n_atoms = 10
    readerclass = AmazingReader


class TestSingleFrameReader(_Single, _TestReader):
    def test_next(self):
        assert_raises(IOError, self.reader.next)

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
