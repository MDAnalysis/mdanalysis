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
    values = np.array([7, 3, 69, 9993, 84, 194, 263, 501, 109, 5873])

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
