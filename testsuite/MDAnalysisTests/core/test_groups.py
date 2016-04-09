import itertools
import numpy as np
from numpy.testing import (
    assert_,
    assert_raises,
)

from MDAnalysisTests.core.groupbase import make_Universe

import MDAnalysis as mda


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

        lvls = self.length.keys()

        nparrays = {
            lvl:np.arange(self.length[lvl])
            for lvl in lvls
        }
        groups = {
            'atom':u.atoms,
            'residue':u.residues,
            'segment':u.segments
        }

        for group in [u.atoms, u.residues, u.segments]:
            yield self._check_len, group
            for func in [list, np.array]:
                yield self._check_boolean_slicing, group, func
            yield self._check_indexerror, group

        for lvl in lvls:
            # Check slicing using slice objects
            for sl in (
                    slice(0, 10),
                    slice(0, 2),
                    slice(1, 3),
                    slice(0, 2, 2),
                    slice(0, -1),
                    slice(5, 1, -1),
                    slice(10, 0, -2),
            ):
                yield self._check_slice, groups[lvl], nparrays[lvl], sl

            # Check slicing using lists and arrays of integers
            for func in [list, lambda x: np.array(x, dtype=np.int64)]:
                for idx in (
                        [0, 1],
                        [0, 0, 0, 1],
                        [3, 2, 1],
                        [-1, -2, -3],
                        [],
                ):
                    yield self._check_slice, group, nparrays[lvl], func(idx)

        singulars = {
            'atom':mda.core.groups.Atom,
            'residue':mda.core.groups.Residue,
            'segment':mda.core.groups.Segment
        }
        for lvl in lvls:
            # Check integer getitem access
            for idx in [0, 1, -1, -2]:
                yield (self._check_integer_getitem, groups[lvl],
                       nparrays[lvl], idx, singulars[lvl])

    def _check_len(self, group):
        assert_(len(group) == self.length[group.level])

    def _check_boolean_slicing(self, group, func):
        # func is the container type that will be used to slice
        group = group[:5]
        sli = func([True, False, False, True, True])
        result = group[sli]
        assert_(len(result) == 3)
        for ref, val in zip(sli, group):
            if ref:
                assert_(val in result)
            else:
                assert_(not val in result)

    def _check_indexerror(self, group):
        assert_raises(IndexError, group.__getitem__, self.length[group.level])

    def _check_slice(self, group, other, sl):
        """Check that slicing a np array is identical"""
        g2 = group[sl]
        o2 = other[sl]

        assert_(len(g2) == len(o2))
        # Check identity of items in the sliced result
        for o, g in zip(o2, g2):
            if o in other:
                assert_(g in g2)
            else:
                assert_(not g in g2)

    def _check_integer_getitem(self, group, nparray, idx, singular):
        a = group[idx]
        ref = nparray[idx]

        assert_(a.ix == ref)
        assert_(isinstance(a, singular))


class TestGroupAddition(object):
    """Tests for combining Group objects

    A + A -> AG
    AG + A -> AG
    A + AG -> AG
    AG + AG -> AG

    SUMMATION
    """
    def test_addition(self):
        u = make_Universe()

        levels = ['atom', 'residue', 'segment']
        groups = {
            'atom':u.atoms[:5],
            'residue':u.residues[:5],
            'segment':u.segments[:5],
        }
        singles = {
            'atom':u.atoms[0],
            'residue':u.residues[0],
            'segment':u.segments[0],
        }

        groupclasses = {
            'atom':mda.core.groups.AtomGroup,
            'residue':mda.core.groups.ResidueGroup,
            'segment':mda.core.groups.SegmentGroup,
        }
        singleclasses = {
            'atom':mda.core.groups.Atom,
            'residue':mda.core.groups.Residue,
            'segment':mda.core.groups.Segment
        }

        for level in levels:
            group = groups[level]
            single = singles[level]
            # check that all combinations of group and singular work
            for x, y in itertools.product([group, single], repeat=2):
                yield self._check_addition, x, y, groupclasses[level]

        # Check that you can't add anything together cross-level
        for alevel, blevel in itertools.permutations(levels, 2):
            for typeA, typeB in itertools.product([singles, groups], repeat=2):
                yield self._check_crosslevel, typeA[alevel], typeB[blevel]
            ### A AG R RG
            # A R
            # A RG
            # AG R
            # AG RG

    @staticmethod
    def itr(x):
        # singular objects don't iterate
        try:
            x[0]
        except TypeError:
            return [x]
        else:
            return x

    def _check_addition(self, a, b, refclass):
        """Combine a and b, check length, returned type and ordering"""
        newgroup = a + b
        reflen = len(self.itr(a)) + len(self.itr(b))
        assert_(len(newgroup) == reflen)
        assert_(isinstance(newgroup, refclass))
        # Check ordering of created Group
        for x, y in zip(newgroup, itertools.chain(self.itr(a), self.itr(b))):
            assert_(x == y)

    def _check_crosslevel(self, a, b):
        def add(x, y):
            return x + y
        assert_raises(TypeError, add, a, b)


class TestGroupLevelTransition(object):
    """Test moving between different hierarchy levels

    AtomGroup
    ^
    v
    ResidueGroup
    ^
    v
    SegmentGroup

    *group_to_*group tests moves between levels
    _unique tests check that Upshifts only return unique higher level
    _listcomp tests check that Downshifts INCLUDE repeated elements
    _unique tests the unique method (performs set operation on self)
    """
    def setUp(self):
        self.u = make_Universe()

    def tearDown(self):
        del self.u

    def test_atomgroup_to_atomgroup(self):
        atm = self.u.atoms.atoms
        assert_(len(atm) == 125)
        assert_(isinstance(atm, mda.core.groups.AtomGroup))

    def test_atomgroup_to_residuegroup(self):
        res = self.u.atoms.residues
        assert_(len(res) == 25)
        assert_(isinstance(res, mda.core.groups.ResidueGroup))

    def test_atomgroup_to_segmentgroup(self):
        seg = self.u.atoms.segments
        assert_(len(seg) == 5)
        assert_(isinstance(seg, mda.core.groups.SegmentGroup))

    def test_residuegroup_to_atomgroup(self):
        atm = self.u.residues.atoms
        assert_(len(atm) == 125)
        assert_(isinstance(atm, mda.core.groups.AtomGroup))

    def test_residuegroup_to_residuegroup(self):
        res = self.u.residues.residues
        assert_(len(res) == 25)
        assert_(isinstance(res, mda.core.groups.ResidueGroup))

    def test_residuegroup_to_segmentgroup(self):
        seg = self.u.residues.segments
        assert_(len(seg) == 5)
        assert_(isinstance(seg, mda.core.groups.SegmentGroup))

    def test_segmentgroup_to_atomgroup(self):
        atm = self.u.segments.atoms
        assert_(len(atm) == 125)
        assert_(isinstance(atm, mda.core.groups.AtomGroup))

    def test_segmentgroup_to_residuegroup(self):
        res = self.u.segments.residues
        assert_(len(res) == 25)
        assert_(isinstance(res, mda.core.groups.ResidueGroup))

    def test_segmentgroup_to_segmentgroup(self):
        seg = self.u.segments.segments
        assert_(len(seg) == 5)
        assert_(isinstance(seg, mda.core.groups.SegmentGroup))

    def test_atom_to_residue(self):
        res = self.u.atoms[0].residue
        assert_(isinstance(res, mda.core.groups.Residue))

    def test_atom_to_segment(self):
        seg = self.u.atoms[0].segment
        assert_(isinstance(seg, mda.core.groups.Segment))

    def test_residue_to_atomgroup(self):
        ag = self.u.residues[0].atoms
        assert_(isinstance(ag, mda.core.groups.AtomGroup))

    def test_residue_to_segment(self):
        seg = self.u.residues[0].segment
        assert_(isinstance(seg, mda.core.groups.Segment))

    def test_segment_to_atomgroup(self):
        ag = self.u.segments[0].atoms
        assert_(isinstance(ag, mda.core.groups.AtomGroup))

    def test_segment_to_residuegroup(self):
        rg = self.u.segments[0].residues
        assert_(isinstance(rg, mda.core.groups.ResidueGroup))

    def test_atomgroup_to_residuegroup_unique(self):
        ag = self.u.atoms[:5] + self.u.atoms[10:15] + self.u.atoms[:5]

        assert_(len(ag.residues) == 2)

    def test_atomgroup_to_segmentgroup_unique(self):
        ag = self.u.atoms[0] + self.u.atoms[-1] + self.u.atoms[0]

        assert_(len(ag.segments) == 2)

    def test_residuegroup_to_segmentgroup_unique(self):
        rg = self.u.residues[0] + self.u.residues[6] + self.u.residues[1]

        assert_(len(rg.segments) == 2)

    def test_residuegroup_to_atomgroup_listcomp(self):
        rg = self.u.residues[0] + self.u.residues[0] + self.u.residues[4]

        assert_(len(rg.atoms) == 15)

    def test_segmentgroup_to_residuegroup_listcomp(self):
        sg = self.u.segments[0] + self.u.segments[0] + self.u.segments[1]

        assert_(len(sg.residues) == 15)

    def test_segmentgroup_to_atomgroup_listcomp(self):
        sg = self.u.segments[0] + self.u.segments[0] + self.u.segments[1]

        assert_(len(sg.atoms) == 75)

    def test_atomgroup_unique(self):
        ag = self.u.atoms[:10] + self.u.atoms[:10]

        assert_(len(ag) == 20)
        assert_(len(ag.unique) == 10)
        assert_(isinstance(ag.unique, mda.core.groups.AtomGroup))

    def test_residuegroup_unique(self):
        rg = self.u.residues[:5] + self.u.residues[:5]

        assert_(len(rg) == 10)
        assert_(len(rg.unique) == 5)
        assert_(isinstance(rg.unique, mda.core.groups.ResidueGroup))

    def test_segmentgroup_unique(self):
        sg = self.u.segments[0] + self.u.segments[1] + self.u.segments[0]

        assert_(len(sg) == 3)
        assert_(len(sg.unique) == 2)
        assert_(isinstance(sg.unique, mda.core.groups.SegmentGroup))
