# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
from six.moves import range

import itertools
import numpy as np
from numpy.testing import (
    dec,
    assert_,
    assert_array_equal,
    assert_equal,
    assert_raises,
    assert_warns,
)
import operator
import six

import MDAnalysis as mda
from MDAnalysisTests import make_Universe, parser_not_found, assert_nowarns
from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysis.core import groups
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import Segids
from MDAnalysis.topology.base import change_squash


class TestGroupProperties(object):
    """ Test attributes of all groups
    """

    def setUp(self):
        self.u = make_Universe(trajectory=True)
        self.group_dict = {
            'atom': self.u.atoms,
            'residue': self.u.residues,
            'segment': self.u.segments
        }

    def test_dimensions(self):
        dimensions = np.arange(6)

        for group in six.itervalues(self.group_dict):
            group.dimensions = dimensions.copy()
            assert_array_equal(group.dimensions, dimensions)
            assert_equal(self.u.dimensions, group.dimensions)



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



class TestReprs(object):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    def test_atom_repr(self):
        at = self.u.atoms[0]
        assert_(repr(at) == '<Atom 1: N of type 56 of resname MET, resid 1 and segid 4AKE>')

    def test_residue_repr(self):
        res = self.u.residues[0]
        assert_(repr(res) == '<Residue MET, 1>')

    def test_segment_repr(self):
        seg = self.u.segments[0]
        assert_(repr(seg) == '<Segment 4AKE>')

    def test_atomgroup_repr(self):
        ag = self.u.atoms[:10]
        assert_(repr(ag) == '<AtomGroup with 10 atoms>')

    def test_atomgroup_str_short(self):
        ag = self.u.atoms[:2]
        assert_(str(ag) == '<AtomGroup [<Atom 1: N of type 56 of resname MET, resid 1 and segid 4AKE>, <Atom 2: HT1 of type 2 of resname MET, resid 1 and segid 4AKE>]>')

    def test_atomgroup_str_long(self):
        ag = self.u.atoms[:11]
        assert_(str(ag).startswith('<AtomGroup [<Atom 1: N of type 56 of resname MET,'))
        assert_('...' in str(ag))
        assert_(str(ag).endswith(', resid 1 and segid 4AKE>]>'))

    def test_residuegroup_repr(self):
        rg = self.u.residues[:10]
        assert_(repr(rg) == '<ResidueGroup with 10 residues>')

    def test_residuegroup_str_short(self):
        rg = self.u.residues[:2]
        assert_(str(rg) == '<ResidueGroup [<Residue MET, 1>, <Residue ARG, 2>]>')

    def test_residuegroup_str_long(self):
        rg = self.u.residues[:11]
        assert_(str(rg).startswith('<ResidueGroup [<Residue MET, 1>,'))
        assert_('...' in str(rg))
        assert_(str(rg).endswith(', <Residue ALA, 11>]>'))

    def test_segmentgroup_repr(self):
        sg = self.u.segments[:10]
        assert_(repr(sg) == '<SegmentGroup with 1 segment>')

    def test_segmentgroup_str(self):
        sg = self.u.segments[:10]
        assert_(str(sg) == '<SegmentGroup [<Segment 4AKE>]>')



class TestGroupBaseOperators(object):
    @staticmethod
    def _test_len(a, b, c, d, e):
        assert_equal(len(a), 4)
        assert_equal(len(b), 5)
        assert_equal(len(c), 2)
        assert_equal(len(d), 0)
        assert_equal(len(e), 3)

    @staticmethod
    def _test_len_duplicated_and_scrambled(a, b, c, d, e):
        assert_equal(len(a), 7)
        assert_equal(len(b), 8)
        assert_equal(len(c), 6)
        assert_equal(len(d), 0)
        assert_equal(len(e), 5)

    @staticmethod
    def _test_equal(a, b, c, d, e):
        assert_(a == a)
        assert_(a != b)
        assert_(not a == b)
        assert_(not a[0:1] == a[0],
                'Element should not equal single element group.')

    @staticmethod
    def _test_issubset(a, b, c, d, e):
        assert_(c.issubset(a))
        assert_(not c.issubset(e))
        assert_(not a.issubset(c))
        assert_(d.issubset(a))
        assert_(not a.issubset(d))

    @staticmethod
    def _test_is_strict_subset(a, b, c, d, e):
        assert_(c.is_strict_subset(a))
        assert_(not c.is_strict_subset(e))
        assert_(not a.is_strict_subset(a))

    @staticmethod
    def _test_issuperset(a, b, c, d, e):
        assert_(a.issuperset(c))
        assert_(not e.issuperset(c))
        assert_(not c.issuperset(a))
        assert_(a.issuperset(d))
        assert_(not d.issuperset(a))

    @staticmethod
    def _test_is_strict_superset(a, b, c, d, e):
        assert_(a.is_strict_superset(c))
        assert_(not c.is_strict_superset(e))
        assert_(not a.is_strict_superset(a))

    @staticmethod
    def _test_concatenate(a, b, c, d, e):
        cat_ab = a.concatenate(b)
        assert_(cat_ab[:len(a)] == a)
        assert_(cat_ab[len(a):] == b)

        cat_ba = b.concatenate(a)
        assert_(cat_ba[:len(b)] == b)
        assert_(cat_ba[len(b):] == a)

        cat_aa = a.concatenate(a)
        assert_(cat_aa[:len(a)] == a)
        assert_(cat_aa[len(a):] == a)

        cat_ad = a.concatenate(d)
        assert_(cat_ad == a)

        cat_da = d.concatenate(a)
        assert_(cat_da == a)

    @staticmethod
    def _test_union(a, b, c, d, e):
        union_ab = a.union(b)
        assert_(union_ab.ix.tolist() == sorted(union_ab.ix))
        assert_(list(sorted(set(union_ab.ix))) == list(sorted(union_ab.ix)))

        assert_(a.union(b) == b.union(a))
        assert_array_equal(a.union(a).ix, np.arange(1, 5))
        assert_(a.union(d), np.arange(1, 5))

    @staticmethod
    def _test_intersection(a, b, c, d, e):
        intersect_ab = a.intersection(b)
        assert_array_equal(intersect_ab.ix, np.arange(3, 5))
        assert_(a.intersection(b) == b.intersection(a))
        assert_equal(len(a.intersection(d)), 0)

    @staticmethod
    def _test_subtract(a, b, c, d, e):
        subtract_ab = a.subtract(b)
        assert_array_equal(subtract_ab.ix, np.array([1, 2, 1, 2]))
        subtract_ba = b.subtract(a)
        assert_array_equal(subtract_ba, np.array([7, 6, 5, 7]))
        subtract_ad = a.subtract(d)
        assert_equal(subtract_ad, a)
        subtract_ae = a.subtract(e)
        assert_equal(subtract_ae, a)

    @staticmethod
    def _test_difference(a, b, c, d, e):
        difference_ab = a.difference(b)
        assert_array_equal(difference_ab.ix, np.arange(1, 3))

        difference_ba = b.difference(a)
        assert_array_equal(difference_ba.ix, np.arange(5, 8))

        assert_array_equal(a.difference(d).ix, np.arange(1, 5))
        assert_array_equal(a.difference(e).ix, np.arange(1, 5))

    @staticmethod
    def _test_symmetric_difference(a, b, c, d, e):
        symdiff_ab = a.symmetric_difference(b)
        assert_array_equal(symdiff_ab.ix, np.array(list(range(1, 3)) +
                                                   list(range(5, 8))))
        assert_(a.symmetric_difference(b) == b.symmetric_difference(a))
        assert_array_equal(a.symmetric_difference(e).ix, np.arange(1, 8))

    @staticmethod
    def _test_isdisjoint(a, b, c, d, e):
        assert_(a.isdisjoint(e))
        assert_(e.isdisjoint(a))
        assert_(a.isdisjoint(d))
        assert_(d.isdisjoint(a))
        assert_(not a.isdisjoint(b))

    @staticmethod
    def make_groups(u, level):
        #   0123456789
        # a  ****
        # b    *****
        # c    **
        # e      ***
        # d empty
        #
        # None of the group start at 0, nor ends at the end. Each group
        # has a different size. The end of a slice is not the last element.
        # This increase the odds of catching errors.
        a = getattr(u, level)[1:5]
        b = getattr(u, level)[3:8]
        c = getattr(u, level)[3:5]
        d = getattr(u, level)[0:0]
        e = getattr(u, level)[5:8]
        return a, b, c, d, e

    @staticmethod
    def make_groups_duplicated_and_scrambled(u, level):
        # The content of the groups is the same as for make_groups, but the
        # elements can appear several times and their order is scrambled.
        a = getattr(u, level)[[1, 3, 2, 1, 2, 4, 4]]
        b = getattr(u, level)[[7, 4, 4, 6, 5, 3, 7, 6]]
        c = getattr(u, level)[[4, 4, 3, 4, 3, 3]]
        d = getattr(u, level)[0:0]
        e = getattr(u, level)[[6, 5, 7, 7, 6]]
        return a, b, c, d, e

    def test_groupbase_operators(self):
        n_segments = 10
        n_residues = n_segments * 5
        n_atoms = n_residues * 5
        u = make_Universe(size=(n_atoms, n_residues, n_segments))
        for level in ('atoms', 'residues', 'segments'):
            a, b, c, d, e = self.make_groups(u, level)
            yield self._test_len, a, b, c, d, e
            yield self._test_equal, a, b, c, d, e
            yield self._test_concatenate, a, b, c, d, e
            yield self._test_union, a, b, c, d, e
            yield self._test_intersection, a, b, c, d, e
            yield self._test_difference, a, b, c, d, e
            yield self._test_symmetric_difference, a, b, c, d, e
            yield self._test_issubset, a, b, c, d, e
            yield self._test_is_strict_subset, a, b, c, d, e
            yield self._test_issuperset, a, b, c, d, e
            yield self._test_is_strict_superset, a, b, c, d, e
            yield self._test_isdisjoint, a, b, c, d, e

    def test_groupbase_operators_duplicated_and_scrambled(self):
        n_segments = 10
        n_residues = n_segments * 5
        n_atoms = n_residues * 5
        u = make_Universe(size=(n_atoms, n_residues, n_segments))
        for level in ('atoms', 'residues', 'segments'):
            a, b, c, d, e = self.make_groups_duplicated_and_scrambled(u, level)
            yield self._test_len_duplicated_and_scrambled, a, b, c, d, e
            yield self._test_equal, a, b, c, d, e
            yield self._test_concatenate, a, b, c, d, e
            yield self._test_union, a, b, c, d, e
            yield self._test_intersection, a, b, c, d, e
            yield self._test_difference, a, b, c, d, e
            yield self._test_symmetric_difference, a, b, c, d, e
            yield self._test_issubset, a, b, c, d, e
            yield self._test_is_strict_subset, a, b, c, d, e
            yield self._test_issuperset, a, b, c, d, e
            yield self._test_is_strict_superset, a, b, c, d, e
            yield self._test_isdisjoint, a, b, c, d, e

    def test_only_same_level(self):
        def dummy(self, other):
            return True

        def failing_pairs(left, right):
            assert_raises(TypeError, _only_same_level(dummy), left, right)

        def succeeding_pairs(left, right):
            assert_(_only_same_level(dummy)(left, right))

        _only_same_level = mda.core.groups._only_same_level
        u = make_Universe()

        components = (u.atoms[0], u.residues[0], u.segments[0])
        groups = (u.atoms, u.residues, u.segments)

        # Do inter-levels pairs of groups fail as expected?
        for left, right in itertools.permutations(groups, 2):
            yield failing_pairs, left, right

        # Do inter-levels pairs of components
        for left, right in itertools.permutations(components, 2):
            yield failing_pairs, left, right

        # Do inter-levels pairs of components/groups fail as expected?
        indexes = range(len(groups))
        for idx_left, idx_right in itertools.permutations(indexes, 2):
            left = groups[idx_left]
            right = groups[idx_right]
            yield failing_pairs, left, right
            yield failing_pairs, right, left

        # Do succeeding pair actually succeed
        for level in ('atoms', 'residues', 'segments'):
            # Groups
            left = getattr(u, level)[0:2]
            right = getattr(u, level)[1:3]
            yield succeeding_pairs, left, right

            # Components
            left = getattr(u, level)[0]
            right = getattr(u, level)[1]
            yield succeeding_pairs, left, right

            # Mixed
            left = getattr(u, level)[0:2]
            right = getattr(u, level)[1]
            yield succeeding_pairs, left, right
            yield succeeding_pairs, right, left

        # Does the function fail with inputs that are not components or groups
        yield failing_pairs, u.atoms, 'invalid'

    def test_only_same_level_different_universes(self):
        def dummy(self, other):
            return True

        u = make_Universe()
        u2 = make_Universe()
        _only_same_level = mda.core.groups._only_same_level
        assert_raises(ValueError, _only_same_level(dummy), u.atoms, u2.atoms)

    def test_shortcut_overriding(self):
        def check_operator(op, method, level):
            left = getattr(u, level)[1:3]
            right = getattr(u, level)[2:4]
            assert_equal(op(left, right), getattr(left, method)(right))

        operators = (
            (operator.add, 'concatenate'),
            (operator.sub, 'difference'),
            (operator.and_, 'intersection'),
            (operator.or_, 'union'),
            (operator.xor, 'symmetric_difference'),
        )
        levels = ('atoms', 'residues', 'segments')

        n_segments = 5
        n_residues = n_segments * 3
        n_atoms = n_residues * 3
        u = make_Universe(size=(n_atoms, n_residues, n_segments))

        for op, method in operators:
            for level in levels:
                yield check_operator, op, method, level


class TestGroupHash(object):
    """
    Groups should be hashable.

    See issue #1397
    """
    def test_hash_exists(self):
        def _hash_type(group):
            assert_(isinstance(hash(group), int))

        u = make_Universe(size=(3, 3, 3))
        for level in ('atoms', 'residues', 'segments'):
            group = getattr(u, level)
            yield _hash_type, group

    def test_hash_equality(self):
        def _hash_equal(a, b):
            assert_equal(hash(a), hash(b))

        u = make_Universe(size=(3, 3, 3))
        for level in ('atoms', 'residues', 'segments'):
            a = getattr(u, level)[0:-1]
            b = getattr(u, level)[0:-1]
            yield _hash_equal, a, b

    def test_hash_difference(self):
        def _hash_not_equal(a, b):
            assert_(hash(a) != hash(b))

        u = make_Universe(size=(3, 3, 3))
        for level in ('atoms', 'residues', 'segments'):
            a = getattr(u, level)[:-1]
            b = getattr(u, level)[1:]
            yield _hash_not_equal, a, b

    def test_hash_difference_cross(self):
        def _hash_not_equal(a, b):
            assert_(hash(a) != hash(b))

        u = make_Universe(size=(3, 3, 3))
        levels = ('atoms', 'residues', 'segments')
        for level_a, level_b in itertools.permutations(levels, 2):
            a = getattr(u, level_a)[0:-1]
            b = getattr(u, level_b)[0:-1]
            yield _hash_not_equal, a, b

    def test_hash_diff_cross_universe(self):
        def _hash_not_equal(a, b):
            assert_(hash(a) != hash(b))

        u = make_Universe(size=(3, 3, 3))
        u2 = make_Universe(size=(3, 3, 3))
        for level in ('atoms', 'residues', 'segments'):
            a = getattr(u, level)
            b = getattr(u2, level)
            yield _hash_not_equal, a, b


class TestAtomGroup(object):

    @staticmethod
    def test_PDB_atom_repr():
        u = make_Universe(extras=('altLocs', 'names', 'types', 'resnames', 'resids', 'segids'))
        assert_equal("<Atom 1: AAA of type TypeA of resname RsA, resid 1 and segid SegA and altLoc A>", u.atoms[0].__repr__())


class TestInstantSelectorDeprecationWarnings(object):
    def setUp(self):
        self.u = make_Universe(("resids", "resnames", "segids", "names"))

    def test_AtomGroup_warn_getitem(self):
        name = self.u.atoms[0].name
        assert_warns(DeprecationWarning, lambda x: self.u.atoms[x], name)

    def test_AtomGroup_nowarn_getitem_index(self):
        assert_nowarns(DeprecationWarning, lambda x: self.u.atoms[x], 0)

    def test_AtomGroup_nowarn_segids_attribute(self):
        assert_nowarns(DeprecationWarning, lambda x: getattr(self.u.atoms, x), "segids")

    def test_AtomGroup_warn_getattr(self):
        name = self.u.atoms[0].name
        assert_warns(DeprecationWarning, lambda x: getattr(self.u.atoms, x), name)

    def test_ResidueGroup_warn_getattr_resname(self):
        name = self.u.residues[0].resname
        assert_warns(DeprecationWarning, lambda x: getattr(self.u.residues, x), name)

    def test_Segment_warn_getattr_resname(self):
        name = self.u.residues[0].resname
        assert_warns(DeprecationWarning, lambda x: getattr(self.u.segments[0], x), name)

    def test_Segment_warn_getattr_rRESNUM(self):
        assert_warns(DeprecationWarning, lambda x: getattr(self.u.segments[0], x), 'r1')

    def test_SegmentGroup_warn_getattr(self):
        name = self.u.segments[0].segid
        assert_warns(DeprecationWarning, lambda x: getattr(self.u.segments, x), name)

    def test_SegmentGroup_nowarn_getitem(self):
        assert_nowarns(DeprecationWarning, lambda x: self.u.segments[x], 0)

