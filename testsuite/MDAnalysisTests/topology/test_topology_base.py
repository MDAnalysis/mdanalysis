import numpy as np
from numpy.testing import (
    assert_array_equal,
)

from MDAnalysis.topology.base import squash_by

class TestSquash(object):
    atom_resids = np.array([2, 2, 1, 1, 5, 5, 4, 4])
    atom_resnames = np.array(['A', 'A', 'B', 'B', 'C', 'C', 'D', 'D'],
                             dtype=object)

    def test_squash(self):
        atom_residx, resids, (resnames,) = squash_by(
            self.atom_resids, self.atom_resnames)

        assert_array_equal(atom_residx, np.array([1, 1, 0, 0, 3, 3, 2, 2]))
        assert_array_equal(resids, np.array([1, 2, 4, 5]))
        assert_array_equal(resnames, np.array(['B', 'A', 'D', 'C']))
