import numpy as np
from numpy.testing import (
    assert_,
    assert_array_almost_equal,
)

from MDAnalysisTests.datafiles import MMTF, MMTF_gz

from MDAnalysis.coordinates.MMTF import MMTFReader


class TestMMTFReader(object):
    def setUp(self):
        self.r = MMTFReader(MMTF)

    def tearDown(self):
        del self.r

    def test_read_frame_size(self):
        assert_(self.r.ts.n_atoms == 512)

    def test_read_positions(self):
        assert_array_almost_equal(self.r.ts.positions[0],
                                  np.array([-0.798, 12.632, 23.231]),
                                  decimal=4)
        assert_array_almost_equal(self.r.ts.positions[-1],
                                  np.array([10.677, 15.517, 11.1]),
                                  decimal=4)

    def test_velocities(self):
        assert_(not self.r.ts.has_velocities)

    def test_forces(self):
        assert_(not self.r.ts.has_forces)

    def test_len(self):
        # should be single frame
        assert_(len(self.r) == 1)

class TestMMTFReaderGZ(object):
    def setUp(self):
        self.r = MMTFReader(MMTF_gz)

    def tearDown(self):
        del self.r

    def test_read_frame_size(self):
        assert_(self.r.ts.n_atoms == 1140)

    def test_read_positions(self):
        assert_array_almost_equal(self.r.ts.positions[0],
                                  np.array([38.428, 16.440, 28.841]),
                                  decimal=4)
        assert_array_almost_equal(self.r.ts.positions[-1],
                                  np.array([36.684, 27.024, 20.468]),
                                  decimal=4)

    def test_velocities(self):
        assert_(not self.r.ts.has_velocities)

    def test_forces(self):
        assert_(not self.r.ts.has_forces)

    def test_len(self):
        # should be single frame
        assert_(len(self.r) == 1)
