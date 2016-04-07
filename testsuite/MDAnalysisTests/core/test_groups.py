import numpy as np
from numpy.testing import (
    assert_,
)

from MDAnalysisTests.core.groupbase import make_Universe

from MDAnalysis.core import levels

class TestGroupSlicing(object):
    length = {'atom':125, 'residue':25, 'segment':5}

    def test_groups(self):
        u = make_Universe()

        for group in [u.atoms, u.residues, u.segments]:
            yield self._check_len, group
            yield self._check_boolean_slicing, group, list
            yield self._check_boolean_slicing, group, np.array
            yield self._check_integer_slicing, group

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

    # boolean indexing
    # integer indexing
    # slicing
