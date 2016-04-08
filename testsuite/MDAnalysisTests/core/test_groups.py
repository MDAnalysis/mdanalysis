import numpy as np
from numpy.testing import (
    assert_,
)

from MDAnalysisTests.core.groupbase import make_Universe

from MDAnalysis.core import levels

class TestGroupSlicing(object):
    """All Groups (Atom, Residue, Segment) should slice like a numpy array

    TODO
    ----
    TopologyGroup is technically called group, add this in too!
    """
    # test Universe is 5:1 mapping 3 times
    length = {'atom':125, 'residue':25, 'segment':5}

    def test_groups(self):
        u = make_Universe()

        for group in [u.atoms, u.residues, u.segments]:
            yield self._check_len, group
            yield self._check_boolean_slicing, group, list
            yield self._check_boolean_slicing, group, np.array
            yield self._check_integer_slicing, group
        for group, nparray in (
                [u.atoms, np.arange(self.length['atom'])],
                [u.residues, np.arange(self.length['residue'])],
                [u.segments, np.arange(self.length['segment'])],
            ):
            for sl in (slice(0, 10), slice(0, 2), slice(1, 3)):
                yield self._check_slice, group, nparray, sl

    def _check_len(self, group):
        assert_(len(group) == self.length[group.level])

    def _check_boolean_slicing(self, group, func):
        # func is the container type that will be used to slice
        group = group[:5]
        sli = func([True, False, False, True, True])
        assert_(len(group[sli]) == 3)

    def _check_integer_slicing(self, group):
        group = group[:5]
        sli = [0, 1]
        assert_(len(group[sli]) == 2)

    def _check_slice(self, group, other, sl):
        """Check that slicing a np array is identical"""
        g2 = group[sl]
        o2 = other[sl]

        assert_(len(g2) == len(o2))

    # boolean indexing
    # integer indexing
    # slicing
