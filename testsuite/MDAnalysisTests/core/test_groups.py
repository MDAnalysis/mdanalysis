import itertools
import numpy as np
from numpy.testing import (
    assert_,
    assert_array_equal,
    assert_equal,
    assert_raises,
)
import operator

from MDAnalysisTests import make_Universe
from MDAnalysis.core import groups


class TestGroupSlicing(object):
    """All Groups (Atom, Residue, Segment) should slice like a numpy array

    TODO
    ----
    TopologyGroup is technically called group, add this in too!
    """
    def test_groups(self):
        u = make_Universe()

        levels = ['atom', 'residue', 'segment']

        # test Universe is 5:1 mapping 3 times
        length = {'atom': 125,
                  'residue': 25,
                  'segment': 5}
        
        nparrays = {
            level: np.arange(length[level])
            for level in levels
        }
        group_dict = {
            'atom': u.atoms,
            'residue': u.residues,
            'segment': u.segments
        }
        singulars = {
            'atom': groups.Atom,
            'residue': groups.Residue,
            'segment': groups.Segment
        }

        for level in levels:
            group = group_dict[level]
            yield self._check_len, group, length[level]
            for func in [list, np.array]:
                yield self._check_boolean_slicing, group, func
            yield self._check_indexerror, group, length[level]
            yield self._check_n_atoms, group
            yield self._check_n_residues, group
            yield self._check_n_segments, group

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
                yield self._check_slice, group_dict[level], nparrays[level], sl

            # Check slicing using lists and arrays of integers
            for func in [list, lambda x: np.array(x, dtype=np.int64)]:
                for idx in (
                        [0, 1],
                        [0, 0, 0, 1],
                        [3, 2, 1],
                        [-1, -2, -3],
                        [],
                ):
                    yield self._check_slice, group, nparrays[level], func(idx)

            # Check integer getitem access
            for idx in [0, 1, -1, -2]:
                yield (self._check_integer_getitem, group_dict[level],
                       nparrays[level], idx, singulars[level])

    @staticmethod
    def _check_n_atoms(group):
        assert_(len(group.atoms) == group.n_atoms)

    @staticmethod
    def _check_n_residues(group):
        assert_(len(group.residues) == group.n_residues)

    @staticmethod
    def _check_n_segments(group):
        assert_(len(group.segments) == group.n_segments)

    def _check_len(self, group, ref):
        assert_(len(group) == ref)

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
                assert_(val not in result)

    def _check_indexerror(self, group, idx):
        assert_raises(IndexError, group.__getitem__, idx)

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
                assert_(g not in g2)

    def _check_integer_getitem(self, group, nparray, idx, singular):
        a = group[idx]
        ref = nparray[idx]

        assert_(a.ix == ref)
        assert_(isinstance(a, singular))


class TestGroupAddition(object):
    """Tests for combining Group objects

    Contents
    --------
    Addition of Groups should work like list addition
    Addition of Singular objects should make Group
      A + A -> AG
      AG + A -> AG
      A + AG -> AG
      AG + AG -> AG
    Cross level addition (eg AG + RG) raises TypeError
    Sum() should work on an iterable of many same level Components/Groups
    Groups contain items "x in y"
    """
    def test_addition(self):
        u = make_Universe()

        levels = ['atom', 'residue', 'segment']
        group_dict = {
            'atom': u.atoms[:5],
            'residue': u.residues[:5],
            'segment': u.segments[:5],
        }
        singles = {
            'atom': u.atoms[0],
            'residue': u.residues[0],
            'segment': u.segments[0],
        }

        groupclasses = {
            'atom': groups.AtomGroup,
            'residue': groups.ResidueGroup,
            'segment': groups.SegmentGroup,
        }
        # TODO: actually use this
        singleclasses = {
            'atom': groups.Atom,
            'residue': groups.Residue,
            'segment': groups.Segment
        }

        for level in levels:
            group = group_dict[level]
            single = singles[level]
            # check that all combinations of group and singular work
            for x, y in itertools.product([group, single], repeat=2):
                yield self._check_addition, x, y, groupclasses[level]

            for x, y, z in itertools.product([group, single], repeat=3):
                yield self._check_sum, x, y, z, groupclasses[level]
                yield self._check_bad_sum, x, y, z

            yield self._check_contains, group
            yield self._check_contains_false, group
            for olevel in levels:
                if level == olevel:
                    continue
                yield self._check_contains_wronglevel, group, group_dict[olevel]

        # Check that you can't add anything together cross-level
        for alevel, blevel in itertools.permutations(levels, 2):
            for typeA, typeB in itertools.product([singles, group_dict], repeat=2):
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

    def _check_sum(self, a, b, c, refclass):
        # weird hack in radd allows this
        summed = sum([a, b, c])

        assert_(isinstance(summed, refclass))
        assert_equal(len(summed),
                     len(self.itr(a)) + len(self.itr(b)) + len(self.itr(c)))
        for x, y in zip(summed,
                        itertools.chain(self.itr(a), self.itr(b), self.itr(c))):
            assert_(x == y)

    @staticmethod
    def _check_bad_sum(a, b, c):
        # sum with bad first argument
        assert_raises(TypeError, sum, [10, a, b, c])

    def _check_crosslevel(self, a, b):
        def add(x, y):
            return x + y
        assert_raises(TypeError, add, a, b)

    def _check_contains(self, group):
        assert_(group[2] in group)

    def _check_contains_false(self, group):
        assert_(not group[3] in group[:2])

    def _check_contains_wronglevel(self, group, group2):
        assert_(not group[2] in group2)


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
        assert_(isinstance(atm, groups.AtomGroup))

    def test_atomgroup_to_residuegroup(self):
        res = self.u.atoms.residues
        assert_(len(res) == 25)
        assert_(isinstance(res, groups.ResidueGroup))

    def test_atomgroup_to_segmentgroup(self):
        seg = self.u.atoms.segments
        assert_(len(seg) == 5)
        assert_(isinstance(seg, groups.SegmentGroup))

    def test_residuegroup_to_atomgroup(self):
        atm = self.u.residues.atoms
        assert_(len(atm) == 125)
        assert_(isinstance(atm, groups.AtomGroup))

    def test_residuegroup_to_residuegroup(self):
        res = self.u.residues.residues
        assert_(len(res) == 25)
        assert_(isinstance(res, groups.ResidueGroup))

    def test_residuegroup_to_segmentgroup(self):
        seg = self.u.residues.segments
        assert_(len(seg) == 5)
        assert_(isinstance(seg, groups.SegmentGroup))

    def test_segmentgroup_to_atomgroup(self):
        atm = self.u.segments.atoms
        assert_(len(atm) == 125)
        assert_(isinstance(atm, groups.AtomGroup))

    def test_segmentgroup_to_residuegroup(self):
        res = self.u.segments.residues
        assert_(len(res) == 25)
        assert_(isinstance(res, groups.ResidueGroup))

    def test_segmentgroup_to_segmentgroup(self):
        seg = self.u.segments.segments
        assert_(len(seg) == 5)
        assert_(isinstance(seg, groups.SegmentGroup))

    def test_atom_to_residue(self):
        res = self.u.atoms[0].residue
        assert_(isinstance(res, groups.Residue))

    def test_atom_to_segment(self):
        seg = self.u.atoms[0].segment
        assert_(isinstance(seg, groups.Segment))

    def test_residue_to_atomgroup(self):
        ag = self.u.residues[0].atoms
        assert_(isinstance(ag, groups.AtomGroup))
        assert_(len(ag) == 5)

    def test_residue_to_segment(self):
        seg = self.u.residues[0].segment
        assert_(isinstance(seg, groups.Segment))

    def test_segment_to_atomgroup(self):
        ag = self.u.segments[0].atoms
        assert_(isinstance(ag, groups.AtomGroup))
        assert_(len(ag) == 25)

    def test_segment_to_residuegroup(self):
        rg = self.u.segments[0].residues
        assert_(isinstance(rg, groups.ResidueGroup))
        assert_(len(rg) == 5)

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
        assert_(isinstance(ag.unique, groups.AtomGroup))

    def test_residuegroup_unique(self):
        rg = self.u.residues[:5] + self.u.residues[:5]

        assert_(len(rg) == 10)
        assert_(len(rg.unique) == 5)
        assert_(isinstance(rg.unique, groups.ResidueGroup))

    def test_segmentgroup_unique(self):
        sg = self.u.segments[0] + self.u.segments[1] + self.u.segments[0]

        assert_(len(sg) == 3)
        assert_(len(sg.unique) == 2)
        assert_(isinstance(sg.unique, groups.SegmentGroup))


class TestComponentComparisons(object):
    """Use of operators (< > == != <= >=) with Atom, Residue, and Segment"""
    @staticmethod
    def _check_lt(a, b, c):
        assert_(a < b)
        assert_(a < c)
        assert_(not b < a)
        assert_(not a < a)

    @staticmethod
    def _check_gt(a, b, c):
        assert_(b > a)
        assert_(c > a)
        assert_(not a > c)
        assert_(not a > a)

    @staticmethod
    def _check_ge(a, b, c):
        assert_(b >= a)
        assert_(c >= a)
        assert_(b >= b)
        assert_(not b >= c)

    @staticmethod
    def _check_le(a, b, c):
        assert_(b <= c)
        assert_(b <= b)
        assert_(not b <= a)

    @staticmethod
    def _check_neq(a, b, c):
        assert_(a != b)
        assert_(not a != a)

    @staticmethod
    def _check_eq(a, b, c):
        assert_(a == a)
        assert_(not a == b)

    @staticmethod
    def _check_sorting(a, b, c):
        assert_(sorted([b, a, c]) == [a, b, c])

    @staticmethod
    def _check_crosslevel_cmp(a, b):
        assert_raises(TypeError, operator.lt, a, b)
        assert_raises(TypeError, operator.le, a, b)
        assert_raises(TypeError, operator.gt, a, b)
        assert_raises(TypeError, operator.ge, a, b)

    @staticmethod
    def _check_crosslevel_eq(a, b):
        assert_raises(TypeError, operator.eq, a, b)
        assert_raises(TypeError, operator.ne, a, b)

    def test_comparions(self):
        u = make_Universe()
        for level in [u.atoms, u.residues, u.segments]:
            a, b, c = level[0], level[1], level[2]
            yield self._check_lt, a, b, c
            yield self._check_gt, a, b, c
            yield self._check_ge, a, b, c
            yield self._check_le, a, b, c
            yield self._check_neq, a, b, c
            yield self._check_eq, a, b, c
            yield self._check_sorting, a, b, c
        # Comparing Atoms To Residues etc
        for a, b in itertools.permutations((u.atoms[0], u.residues[0], u.segments[0]), 2):
            yield self._check_crosslevel_cmp, a, b
            yield self._check_crosslevel_eq, a, b

class TestMetaclassMagic(object):
    # tests for the weird voodoo we do with metaclasses
    @staticmethod
    def test_new_class():
        u = make_Universe(trajectory=True)

        # should be able to subclass AtomGroup as normal
        class NewGroup(groups.AtomGroup):
            pass

        ng = NewGroup(np.array([0, 1, 2]), u)

        assert_(isinstance(ng, NewGroup))

        ag = u.atoms[[0, 1, 2]]

        assert_array_equal(ng.positions, ag.positions)


class TestGroupBy(object):
    # tests for the method 'groupby'
    def setUp(self):
        self.u = make_Universe(('types', 'charges', 'resids'))

    def tearDown(self):
        del self.u

    def test_groupby_float(self):
        gb = self.u.atoms.groupby('charges')

        for ref in [-1.5, -0.5, 0.0, 0.5, 1.5]:
            assert_(ref in gb)
            g = gb[ref]
            assert_(all(g.charges == ref))
            assert_(len(g) == 25)
        
    def test_groupby_string(self):
        gb = self.u.atoms.groupby('types')

        assert_(len(gb) == 5)
        for ref in ['TypeA', 'TypeB', 'TypeC', 'TypeD', 'TypeE']:
            assert_(ref in gb)
            g = gb[ref]
            assert_(all(g.types == ref))
            assert_(len(g) == 25)

    def test_groupby_int(self):
        gb = self.u.atoms.groupby('resids')

        for g in gb.values():
            assert_(len(g) == 5)
