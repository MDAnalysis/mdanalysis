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

    def _check_len(self, group):
        assert_(len(group) == self.length[group.level])
