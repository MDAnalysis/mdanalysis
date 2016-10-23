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
