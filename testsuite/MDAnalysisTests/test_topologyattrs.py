"""Tests for MDAnalysis.core.topologyattrs objects.

"""
from numpy.testing import (
    assert_,
    assert_array_equal,
)
import numpy as np

import MDAnalysis.core.topologyattrs as tpattrs


class TestTopologyAttr(object):
    """Test the common elements to all TopologyAttrs.

    """
    ref = np.array([7, 3, 69, 9993])
    
    def setUp(self):
        self.topattr = tpattrs.TopologyAttr(ref)

    def tearDown(self):
        del self.topattr


