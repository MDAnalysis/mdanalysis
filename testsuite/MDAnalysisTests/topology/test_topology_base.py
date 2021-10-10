import numpy as np
from numpy.testing import assert_equal

from MDAnalysis.topology.base import squash_by, change_squash


class TestSquash(object):
    atom_resids = np.array([2, 2, 1, 1, 5, 5, 4, 4])
    atom_resnames = np.array(['A', 'A', 'B', 'B', 'C', 'C', 'D', 'D'],
                             dtype=object)

    def test_squash(self):
        atom_residx, resids, (resnames, ) = squash_by(self.atom_resids,
                                                      self.atom_resnames)

        assert_equal(atom_residx, np.array([1, 1, 0, 0, 3, 3, 2, 2]))
        assert_equal(resids, np.array([1, 2, 4, 5]))
        assert_equal(resnames, np.array(['B', 'A', 'D', 'C']))


class TestChangeSquash(object):
    def test_resid_squash(self):
        # System of 6 atoms in 3 residues
        # Residues 1 & 2 are Segid A, Residue 3 is Segid B
        # Resid 2 is repeated twice! Should be detected as 2 distinct residues
        resids = np.array([2, 2, 3, 3, 2, 2])
        resnames = np.array(['RsA', 'RsA', 'RsB', 'RsB', 'RsC', 'RsC'])
        segids = np.array(['A', 'A', 'A', 'A', 'B', 'B'])

        residx, (new_resids, new_resnames, new_segids) = change_squash(
            (resids, ), (resids, resnames, segids))

        assert_equal(residx, np.array([0, 0, 1, 1, 2, 2]))
        assert_equal(new_resids, np.array([2, 3, 2]))
        assert_equal(new_resnames, np.array(['RsA', 'RsB', 'RsC']))
        assert_equal(new_segids, np.array(['A', 'A', 'B']))

    def test_segid_squash(self):
        segids = np.array(['A', 'A', 'B'])

        segidx, (new_segids, ) = change_squash((segids, ), (segids, ))

        assert_equal(segidx, np.array([0, 0, 1]))
        assert_equal(new_segids, np.array(['A', 'B']))
