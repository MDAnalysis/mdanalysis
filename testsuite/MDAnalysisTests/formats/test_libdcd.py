from __future__ import print_function

from nose.tools import raises
from numpy.testing import assert_equal, assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal

from MDAnalysis.lib.formats.libdcd import DCDFile
from MDAnalysisTests.datafiles import DCD

from unittest import TestCase
from MDAnalysisTests.tempdir import run_in_tempdir
import numpy as np

class DCDReadFrameTest(TestCase):

    def setUp(self):
        self.dcdfile = DCDFile(DCD)

    def tearDown(self):
        del self.dcdfile

    def test_read_coords(self):
        # confirm shape of coordinate data against result from previous
        # MDAnalysis implementation of DCD file handling
        dcd_frame = self.dcdfile.read()
        xyz = dcd_frame[0]
        assert_equal(xyz.shape, (3341, 3))

    def test_read_unit_cell(self):
        # confirm unit cell read against result from previous
        # MDAnalysis implementation of DCD file handling
        dcd_frame = self.dcdfile.read()
        unitcell = dcd_frame[1]
        expected = np.array([  0.,   0.,   0.,  90.,  90.,  90.],
                            dtype=np.float32)
        assert_equal(unitcell, expected)

    def test_seek_over_max(self):
        # should raise IOError if beyond 98th frame
        with self.assertRaises(IOError):
            self.dcdfile.seek(102)

    def test_seek_normal(self):
        # frame seek within range is tested
        new_frame = 91
        self.dcdfile.seek(new_frame)
        assert_equal(self.dcdfile.tell(), new_frame)

    def test_seek_negative(self):
        # frame seek with negative number
        with self.assertRaises(IOError):
            self.dcdfile.seek(-78)

    def test_iteration(self):
        self.dcdfile.__next__()
        self.dcdfile.__next__()
        self.dcdfile.__next__()
        expected_frame = 3
        assert_equal(self.dcdfile.tell(), expected_frame)

    def test_zero_based_frames(self):
        expected_frame = 0
        assert_equal(self.dcdfile.tell(), expected_frame)

    def test_length_traj(self):
        expected = 98
        assert_equal(len(self.dcdfile), expected)

    def test_context_manager(self):
        frame = 22
        with self.dcdfile as f:
            f.seek(frame)
            assert_equal(f.tell(), frame)

    @raises(IOError)
    def test_open_wrong_mode(self):
        DCDFile('foo', 'e')

    @raises(IOError)
    def test_raise_not_existing(self):
        DCDFile('foo')

    def test_n_atoms(self):
        assert_equal(self.dcdfile.n_atoms, 3341)

    @raises(IOError)
    @run_in_tempdir()
    def test_read_write_mode_file(self):
        with DCDFile('foo', 'w') as f:
            f.read()

    @raises(RuntimeError)
    def test_read_closed(self):
        self.dcdfile.close()
        self.dcdfile.read()

    def test_iteration_2(self):
        with self.dcdfile as f:
            for frame in f:
                pass
