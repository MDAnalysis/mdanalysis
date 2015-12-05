"""Tests for MDAnalysis.core.topologyattrs objects.

"""
from numpy.testing import (
    assert_,
    assert_array_equal,
)
import numpy as np

import MDAnalysis.core.topologyattrs as tpattrs
from MDAnalysis.core.topology import Topology


class TopologyAttrMixin(object):
    """Mixin to test the common elements to all TopologyAttrs.

    """
    # Reference data
    Ridx = np.array([0, 0, 2, 2, 1, 1, 3, 3, 1, 2])
    Sidx = np.array([0, 1, 1, 0])

    def setUp(self):
        self.top = Topology(10, 4, 2, attrs=[self.attrclass(self.values)],
                            atom_resindex=self.Ridx,
                            residue_segindex=self.Sidx)
        
        self.attr = getattr(self.top, self.attrclass.attrname)

    def tearDown(self):
        del self.top

    def test_len(self):
        assert len(self.attr) == len(self.attr.values)


class TestAtomAttr(TopologyAttrMixin):
    """Test atom-level TopologyAttrs.

    """
    values = np.array([7, 3, 69, 9993, 84, 194, 263, 501, 109, 5873])
    attrclass = tpattrs.AtomAttr

    def test_get_atoms(self):
        assert_array_equal(self.attr.get_atoms([2, 1]), np.array([69, 3]))

    def test_set_atoms(self):

        # broadcasting
        self.attr.set_atoms([2, 1], 87)
        assert_array_equal(self.attr.get_atoms([2, 1]), np.array([87, 87]))

        # set with array
        self.attr.set_atoms([3, 7], np.array([23, 504]))
        assert_array_equal(self.attr.get_atoms([3, 7]), np.array([23, 504]))

    def test_get_residues(self):
        """Unless overriden by child class, this should yield values for all
        atoms in residues.

        """
        assert_array_equal(self.attr.get_residues([2, 1]), 
                           np.array([69, 9993, 5873, 84, 194, 109]))

    def test_get_segments(self):
        """Unless overriden by child class, this should yield values for all
        atoms in segments.

        """
        assert_array_equal(self.attr.get_segments([1]), 
                           np.array([84, 194, 109, 69, 9993, 5873]))


class TestResidueAttr(TopologyAttrMixin):
    """Test residue-level TopologyAttrs.

    """
    values = np.array([15.2, 395.6, 0.1, 9.8])
    attrclass = tpattrs.ResidueAttr

    def test_get_atoms(self):
        assert_array_equal(self.attr.get_atoms([7, 3, 9]), np.array([9.8, 0.1, 0.1]))

    def test_get_residues(self):
        assert_array_equal(self.attr.get_residues([1, 2, 1, 3]), 
                           np.array([395.6, 0.1, 395.6, 9.8]))

    def test_set_residues(self):

        # broadcasting
        self.attr.set_residues([3, 0, 2], 74)
        assert_array_equal(self.attr.get_residues([3, 0, 2]), np.array([74, 74, 74]))

        # set with array
        self.attr.set_residues([3, 0, 1], np.array([23, 504, 0.0002]))
        assert_array_equal(self.attr.get_residues([3, 0, 1]), np.array([23, 504, 0.0002]))

    def test_get_segments(self):
        """Unless overriden by child class, this should yield values for all
        atoms in segments.

        """
        assert_array_equal(self.attr.get_segments([0, 1, 1]), 
                           np.array([15.2, 9.8, 395.6, 0.1, 395.6, 0.1]))


class TestSegmentAttr(TopologyAttrMixin):
    """Test segment-level TopologyAttrs.

    """
    values = np.array([-0.19, 500])
    attrclass = tpattrs.SegmentAttr

    def test_get_atoms(self):
        assert_array_equal(self.attr.get_atoms([2, 4, 1]), np.array([500, 500, -0.19]))

    def test_get_residues(self):
        assert_array_equal(self.attr.get_residues([1, 2, 1, 3]), 
                           np.array([500, 500, 500, -0.19]))

    def test_get_segments(self):
        """Unless overriden by child class, this should yield values for all
        atoms in segments.

        """
        assert_array_equal(self.attr.get_segments([1, 0, 0]), 
                           np.array([500, -0.19, -0.19]))

    def test_set_segments(self):

        # broadcasting
        self.attr.set_segments([1, 0], -108)
        assert_array_equal(self.attr.get_segments([0, 1]), np.array([-108, -108]))

        # set with array
        self.attr.set_segments([0, 1], np.array([23, -0.0002]))
        assert_array_equal(self.attr.get_segments([1, 0, 1]), 
                np.array([-0.0002, 23, -0.0002]))
