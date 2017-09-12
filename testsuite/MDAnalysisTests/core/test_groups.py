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
    assert_,
    assert_array_equal,
    assert_equal,
    assert_warns,
)
import pytest
import operator
import six

import MDAnalysis as mda
from MDAnalysisTests import make_Universe, assert_nowarns
from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysis.core import groups


class TestGroupProperties(object):
    """ Test attributes of all groups
    """
    @pytest.fixture()
    def u(self):
        return make_Universe(trajectory=True)

    @pytest.fixture()
    def group_dict(self, u):
        return {
            'atom': u.atoms,
            'residue': u.residues,
            'segment': u.segments
        }

    def test_dimensions(self, u, group_dict):
        dimensions = np.arange(6)

        for group in six.itervalues(group_dict):
            group.dimensions = dimensions.copy()
            assert_array_equal(group.dimensions, dimensions)
            assert_equal(u.dimensions, group.dimensions)


class TestGroupSlicing(object):
    """All Groups (Atom, Residue, Segment) should slice like a numpy array

    TODO
    ----
    TopologyGroup is technically called group, add this in too!
    """
    u = make_Universe()

    # test universe is 5:1 mapping 3 times
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
    slices = (
        slice(0, 10),
        slice(0, 2),
        slice(1, 3),
        slice(0, 2, 2),
        slice(0, -1),
        slice(5, 1, -1),
        slice(10, 0, -2),
    )
    length = {'atom': 125,
              'residue': 25,
              'segment': 5}

    levels = ('atom', 'residue', 'segment')

    @pytest.fixture(params=levels)
    def level(self, request):
        return request.param

    @pytest.fixture
    def group(self, level):
        return self.group_dict[level]

    @pytest.fixture
    def nparray(self, level):
        return np.arange(self.length[level])

    @pytest.fixture
    def singular(self, level):
        return self.singulars[level]

    def test_n_atoms(self, group):
        assert_(len(group.atoms) == group.n_atoms)

    def test_n_residues(self, group):
        assert_(len(group.residues) == group.n_residues)

    def test_n_segments(self, group):
        assert_(len(group.segments) == group.n_segments)

    def test_len(self, group, level):
        ref = self.length[level]
        assert_(len(group) == ref)

    @pytest.mark.parametrize('func', [list, np.array])
    def test_boolean_slicing(self, group, func):
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

    def test_indexerror(self, group, level):
        idx = self.length[level]
        with pytest.raises(IndexError):
            group.__getitem__(idx)

    @pytest.mark.parametrize('sl,func', itertools.product((
        slice(0, 10),
        slice(0, 2),
        slice(1, 3),
        slice(0, 2, 2),
        slice(0, -1),
        slice(5, 1, -1),
        slice(10, 0, -2),
    ), [list, lambda x: np.array(x, dtype=np.int64)]))
    def test_slice(self, group, nparray, sl, func):
        """Check that slicing a np array is identical"""
        g2 = group[sl]
        o2 = nparray[sl]

        assert_(len(g2) == len(o2))
        # Check identity of items in the sliced result
        for o, g in zip(o2, g2):
            if o in nparray:
                assert_(g in g2)
            else:
                assert_(g not in g2)

    @pytest.mark.parametrize('idx', [0, 1, -1, -2])
    def test_integer_getitem(self, group, nparray, idx, singular):
        a = group[idx]
        ref = nparray[idx]

        assert_(a.ix == ref)
        assert_(isinstance(a, singular))


def _yield_groups(group_dict, singles, levels, groupclasses, repeat):
    for level in levels:
        for groups in itertools.product([group_dict[level], singles[level]],
                                        repeat=repeat):
            yield  list(groups) + [groupclasses[level]]

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

    @pytest.fixture(params=levels)
    def level(self, request):
        return request.param

    @pytest.fixture
    def group(self, level):
        return self.group_dict[level]

    @pytest.fixture
    def single(self, level):
        return self.singles[level]

    @pytest.fixture
    def two_groups(self, group, single):
        return itertools.product([group, single], repeat=2)

    @pytest.fixture
    def three_groups(self, group, single):
        return itertools.product([group, single], repeat=3)

    @staticmethod
    def itr(x):
        # singular objects don't iterate
        try:
            x[0]
        except TypeError:
            return [x]
        else:
            return x

    @pytest.mark.parametrize(
        'a, b, refclass',
        _yield_groups(group_dict, singles, levels, groupclasses, repeat=2)
    )
    def test_addition(self, a, b, refclass):
        """Combine a and b, check length, returned type and ordering"""
        newgroup = a + b
        reflen = len(self.itr(a)) + len(self.itr(b))
        assert_(len(newgroup) == reflen)
        assert_(isinstance(newgroup, refclass))
        # Check ordering of created Group
        for x, y in zip(newgroup, itertools.chain(self.itr(a), self.itr(b))):
            assert_(x == y)

    @pytest.mark.parametrize(
        'a, b, c, refclass',
        _yield_groups(group_dict, singles, levels, groupclasses, repeat=3)
    )
    def test_sum(self, a, b, c, refclass):
        # weird hack in radd allows this
        summed = sum([a, b, c])

        assert_(isinstance(summed, refclass))
        assert_equal(len(summed),
                     len(self.itr(a)) + len(self.itr(b)) + len(self.itr(c)))
        for x, y in zip(summed,
                        itertools.chain(self.itr(a), self.itr(b), self.itr(c))):
            assert_(x == y)

    @pytest.mark.parametrize(
        'a, b, c, refclass',
        _yield_groups(group_dict, singles, levels, groupclasses, repeat=3)
    )
    def test_bad_sum(self, a, b, c, refclass):
        # sum with bad first argument
        with pytest.raises(TypeError):
            sum([10, a, b, c])

    def test_contains(self, group):
        assert_(group[2] in group)

    def test_contains_false(self, group):
        assert_(not group[3] in group[:2])

    @pytest.mark.parametrize(
        'one_level, other_level',
        [
            (l1, l2)
            for l1, l2
            in itertools.product(levels, repeat=2)
            if l1 != l2
        ]
    )
    def test_contains_wronglevel(self, one_level, other_level):
        group = self.group_dict[one_level]
        group2 = self.group_dict[other_level]
        assert_(not group[2] in group2)

    @pytest.mark.parametrize(
        'a, b',
        [
            (typeA[alevel], typeB[blevel])
            for (typeA, typeB), (alevel, blevel)
            in itertools.product(
                itertools.product([singles, group_dict], repeat=2),
                itertools.permutations(levels, 2)
            )
        ]
    )
    def test_crosslevel(self, a, b):
        with pytest.raises(TypeError):
            a + b


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
        
    @pytest.fixture()
    def u(self):
        return make_Universe()

    def test_atomgroup_to_atomgroup(self, u):
        atm = u.atoms.atoms
        assert_(len(atm) == 125)
        assert_(isinstance(atm, groups.AtomGroup))

    def test_atomgroup_to_residuegroup(self, u):
        res = u.atoms.residues
        assert_(len(res) == 25)
        assert_(isinstance(res, groups.ResidueGroup))

    def test_atomgroup_to_segmentgroup(self, u):
        seg = u.atoms.segments
        assert_(len(seg) == 5)
        assert_(isinstance(seg, groups.SegmentGroup))

    def test_residuegroup_to_atomgroup(self, u):
        atm = u.residues.atoms
        assert_(len(atm) == 125)
        assert_(isinstance(atm, groups.AtomGroup))

    def test_residuegroup_to_residuegroup(self, u):
        res = u.residues.residues
        assert_(len(res) == 25)
        assert_(isinstance(res, groups.ResidueGroup))

    def test_residuegroup_to_segmentgroup(self, u):
        seg = u.residues.segments
        assert_(len(seg) == 5)
        assert_(isinstance(seg, groups.SegmentGroup))

    def test_segmentgroup_to_atomgroup(self, u):
        atm = u.segments.atoms
        assert_(len(atm) == 125)
        assert_(isinstance(atm, groups.AtomGroup))

    def test_segmentgroup_to_residuegroup(self, u):
        res = u.segments.residues
        assert_(len(res) == 25)
        assert_(isinstance(res, groups.ResidueGroup))

    def test_segmentgroup_to_segmentgroup(self, u):
        seg = u.segments.segments
        assert_(len(seg) == 5)
        assert_(isinstance(seg, groups.SegmentGroup))

    def test_atom_to_residue(self, u):
        res = u.atoms[0].residue
        assert_(isinstance(res, groups.Residue))

    def test_atom_to_segment(self, u):
        seg = u.atoms[0].segment
        assert_(isinstance(seg, groups.Segment))

    def test_residue_to_atomgroup(self, u):
        ag = u.residues[0].atoms
        assert_(isinstance(ag, groups.AtomGroup))
        assert_(len(ag) == 5)

    def test_residue_to_segment(self, u):
        seg = u.residues[0].segment
        assert_(isinstance(seg, groups.Segment))

    def test_segment_to_atomgroup(self, u):
        ag = u.segments[0].atoms
        assert_(isinstance(ag, groups.AtomGroup))
        assert_(len(ag) == 25)

    def test_segment_to_residuegroup(self, u):
        rg = u.segments[0].residues
        assert_(isinstance(rg, groups.ResidueGroup))
        assert_(len(rg) == 5)

    def test_atomgroup_to_residuegroup_unique(self, u):
        ag = u.atoms[:5] + u.atoms[10:15] + u.atoms[:5]

        assert_(len(ag.residues) == 2)

    def test_atomgroup_to_segmentgroup_unique(self, u):
        ag = u.atoms[0] + u.atoms[-1] + u.atoms[0]

        assert_(len(ag.segments) == 2)

    def test_residuegroup_to_segmentgroup_unique(self, u):
        rg = u.residues[0] + u.residues[6] + u.residues[1]

        assert_(len(rg.segments) == 2)

    def test_residuegroup_to_atomgroup_listcomp(self, u):
        rg = u.residues[0] + u.residues[0] + u.residues[4]

        assert_(len(rg.atoms) == 15)

    def test_segmentgroup_to_residuegroup_listcomp(self, u):
        sg = u.segments[0] + u.segments[0] + u.segments[1]

        assert_(len(sg.residues) == 15)

    def test_segmentgroup_to_atomgroup_listcomp(self, u):
        sg = u.segments[0] + u.segments[0] + u.segments[1]

        assert_(len(sg.atoms) == 75)

    def test_atomgroup_unique(self, u):
        ag = u.atoms[:10] + u.atoms[:10]

        assert_(len(ag) == 20)
        assert_(len(ag.unique) == 10)
        assert_(isinstance(ag.unique, groups.AtomGroup))

    def test_residuegroup_unique(self, u):
        rg = u.residues[:5] + u.residues[:5]

        assert_(len(rg) == 10)
        assert_(len(rg.unique) == 5)
        assert_(isinstance(rg.unique, groups.ResidueGroup))

    def test_segmentgroup_unique(self, u):
        sg = u.segments[0] + u.segments[1] + u.segments[0]

        assert_(len(sg) == 3)
        assert_(len(sg.unique) == 2)
        assert_(isinstance(sg.unique, groups.SegmentGroup))


class TestComponentComparisons(object):
    """Use of operators (< > == != <= >=) with Atom, Residue, and Segment"""
    u = make_Universe()
    levels = [u.atoms, u.residues, u.segments]

    @pytest.fixture(params=levels)
    def abc(self, request):
        level = request.param
        return level[0], level[1], level[2]

    @pytest.fixture
    def a(self, abc):
        return abc[0]

    @pytest.fixture
    def b (self, abc):
        return abc[1]

    @pytest.fixture
    def c(self, abc):
        return abc[2]

    def test_lt(self, a, b, c):
        assert_(a < b)
        assert_(a < c)
        assert_(not b < a)
        assert_(not a < a)

    def test_gt(self, a, b, c):
        assert_(b > a)
        assert_(c > a)
        assert_(not a > c)
        assert_(not a > a)

    def test_ge(self, a, b, c):
        assert_(b >= a)
        assert_(c >= a)
        assert_(b >= b)
        assert_(not b >= c)

    def test_le(self, a, b, c):
        assert_(b <= c)
        assert_(b <= b)
        assert_(not b <= a)

    def test_neq(self, a, b, c):
        assert_(a != b)
        assert_(not a != a)

    def test_eq(self, a, b, c):
        assert_(a == a)
        assert_(not a == b)

    def test_sorting(self, a, b, c):
        assert_(sorted([b, a, c]) == [a, b, c])

    @pytest.mark.parametrize(
        'x, y',
        itertools.permutations((u.atoms[0], u.residues[0], u.segments[0]), 2)
    )
    def test_crosslevel_cmp(self, x, y):
        with pytest.raises(TypeError):
            operator.lt(x, y)
        with pytest.raises(TypeError):
            operator.le(x, y)
        with pytest.raises(TypeError):
            operator.gt(x, y)
        with pytest.raises(TypeError):
            operator.ge(x, y)

    @pytest.mark.parametrize(
        'x, y',
        itertools.permutations((u.atoms[0], u.residues[0], u.segments[0]), 2)
    )
    def test_crosslevel_eq(self, x, y):
        with pytest.raises(TypeError):
            operator.eq(x, y)

        with pytest.raises(TypeError):
            operator.ne(x, y)


class TestMetaclassMagic(object):
    # tests for the weird voodoo we do with metaclasses
    def test_new_class(self):
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
    @pytest.fixture()
    def u(self):
        return make_Universe(('types', 'charges', 'resids'))

    def test_groupby_float(self, u):
        gb = u.atoms.groupby('charges')

        for ref in [-1.5, -0.5, 0.0, 0.5, 1.5]:
            assert_(ref in gb)
            g = gb[ref]
            assert_(all(g.charges == ref))
            assert_(len(g) == 25)

    def test_groupby_string(self, u):
        gb = u.atoms.groupby('types')

        assert_(len(gb) == 5)
        for ref in ['TypeA', 'TypeB', 'TypeC', 'TypeD', 'TypeE']:
            assert_(ref in gb)
            g = gb[ref]
            assert_(all(g.types == ref))
            assert_(len(g) == 25)

    def test_groupby_int(self, u):
        gb = u.atoms.groupby('resids')

        for g in gb.values():
            assert_(len(g) == 5)


class TestReprs(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(PSF, DCD)

    def test_atom_repr(self, u):
        at = u.atoms[0]
        assert_(repr(at) == '<Atom 1: N of type 56 of resname MET, resid 1 and segid 4AKE>')

    def test_residue_repr(self, u):
        res = u.residues[0]
        assert_(repr(res) == '<Residue MET, 1>')

    def test_segment_repr(self, u):
        seg = u.segments[0]
        assert_(repr(seg) == '<Segment 4AKE>')

    def test_atomgroup_repr(self, u):
        ag = u.atoms[:10]
        assert_(repr(ag) == '<AtomGroup with 10 atoms>')

    def test_atomgroup_str_short(self, u):
        ag = u.atoms[:2]
        assert_(str(ag) == '<AtomGroup [<Atom 1: N of type 56 of resname MET, resid 1 and segid 4AKE>, <Atom 2: HT1 of type 2 of resname MET, resid 1 and segid 4AKE>]>')

    def test_atomgroup_str_long(self, u):
        ag = u.atoms[:11]
        assert_(str(ag).startswith('<AtomGroup [<Atom 1: N of type 56 of resname MET,'))
        assert_('...' in str(ag))
        assert_(str(ag).endswith(', resid 1 and segid 4AKE>]>'))

    def test_residuegroup_repr(self, u):
        rg = u.residues[:10]
        assert_(repr(rg) == '<ResidueGroup with 10 residues>')

    def test_residuegroup_str_short(self, u):
        rg = u.residues[:2]
        assert_(str(rg) == '<ResidueGroup [<Residue MET, 1>, <Residue ARG, 2>]>')

    def test_residuegroup_str_long(self, u):
        rg = u.residues[:11]
        assert_(str(rg).startswith('<ResidueGroup [<Residue MET, 1>,'))
        assert_('...' in str(rg))
        assert_(str(rg).endswith(', <Residue ALA, 11>]>'))

    def test_segmentgroup_repr(self, u):
        sg = u.segments[:10]
        assert_(repr(sg) == '<SegmentGroup with 1 segment>')

    def test_segmentgroup_str(self, u):
        sg = u.segments[:10]
        assert_(str(sg) == '<SegmentGroup [<Segment 4AKE>]>')


def _yield_mix(groups, components):
    indices = list(range(len(components)))
    for left, right in itertools.permutations(indices, 2):
        yield (groups[left], components[right])
        yield (components[left], groups[right])

def _yield_sliced_groups(u, slice_left, slice_right):
    for level in ('atoms', 'residues', 'segments'):
        yield (getattr(u, level)[slice_left], getattr(u, level)[slice_right])

class TestGroupBaseOperators(object):
    u = make_Universe()

    components = (u.atoms[0], u.residues[0], u.segments[0])
    component_groups = (u.atoms, u.residues, u.segments)

    @pytest.fixture(params=('atoms', 'residues', 'segments'))
    def level(self, request):
        return request.param

    @pytest.fixture
    def groups_simple(self, level):
        n_segments = 10
        n_residues = n_segments * 5
        n_atoms = n_residues * 5
        u = make_Universe(size=(n_atoms, n_residues, n_segments))
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

    @pytest.fixture
    def groups_duplicated_and_scrambled(self, level):
        # The content of the groups is the same as for make_groups, but the
        # elements can appear several times and their order is scrambled.
        n_segments = 10
        n_residues = n_segments * 5
        n_atoms = n_residues * 5
        u = make_Universe(size=(n_atoms, n_residues, n_segments))
        a = getattr(u, level)[[1, 3, 2, 1, 2, 4, 4]]
        b = getattr(u, level)[[7, 4, 4, 6, 5, 3, 7, 6]]
        c = getattr(u, level)[[4, 4, 3, 4, 3, 3]]
        d = getattr(u, level)[0:0]
        e = getattr(u, level)[[6, 5, 7, 7, 6]]
        return a, b, c, d, e

    @pytest.fixture(params=('simple', 'scrambled'))
    def groups(self, request, groups_simple, groups_duplicated_and_scrambled):
        return {'simple': groups_simple,
                'scrambled': groups_duplicated_and_scrambled}[request.param]

    def test_len(self, groups_simple):
        a, b, c, d, e = groups_simple
        assert_equal(len(a), 4)
        assert_equal(len(b), 5)
        assert_equal(len(c), 2)
        assert_equal(len(d), 0)
        assert_equal(len(e), 3)

    def test_len_duplicated_and_scrambled(self, groups_duplicated_and_scrambled):
        a, b, c, d, e = groups_duplicated_and_scrambled
        assert_equal(len(a), 7)
        assert_equal(len(b), 8)
        assert_equal(len(c), 6)
        assert_equal(len(d), 0)
        assert_equal(len(e), 5)

    def test_equal(self, groups):
        a, b, c, d, e = groups
        assert_(a == a)
        assert_(a != b)
        assert_(not a == b)
        assert_(not a[0:1] == a[0],
                'Element should not equal single element group.')

    def test_issubset(self, groups):
        a, b, c, d, e = groups
        assert_(c.issubset(a))
        assert_(not c.issubset(e))
        assert_(not a.issubset(c))
        assert_(d.issubset(a))
        assert_(not a.issubset(d))

    def test_is_strict_subset(self, groups):
        a, b, c, d, e = groups
        assert_(c.is_strict_subset(a))
        assert_(not c.is_strict_subset(e))
        assert_(not a.is_strict_subset(a))

    def test_issuperset(self, groups):
        a, b, c, d, e = groups
        assert_(a.issuperset(c))
        assert_(not e.issuperset(c))
        assert_(not c.issuperset(a))
        assert_(a.issuperset(d))
        assert_(not d.issuperset(a))

    def test_is_strict_superset(self, groups):
        a, b, c, d, e = groups
        assert_(a.is_strict_superset(c))
        assert_(not c.is_strict_superset(e))
        assert_(not a.is_strict_superset(a))

    def test_concatenate(self, groups):
        a, b, c, d, e = groups
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

    def test_union(self, groups):
        a, b, c, d, e = groups
        union_ab = a.union(b)
        assert_(union_ab.ix.tolist() == sorted(union_ab.ix))
        assert_(list(sorted(set(union_ab.ix))) == list(sorted(union_ab.ix)))

        assert_(a.union(b) == b.union(a))
        assert_array_equal(a.union(a).ix, np.arange(1, 5))
        assert_(a.union(d), np.arange(1, 5))

    def test_intersection(self, groups):
        a, b, c, d, e = groups
        intersect_ab = a.intersection(b)
        assert_array_equal(intersect_ab.ix, np.arange(3, 5))
        assert_(a.intersection(b) == b.intersection(a))
        assert_equal(len(a.intersection(d)), 0)

    def test_subtract(self, groups):
        a, b, c, d, e = groups
        subtract_ab = a.subtract(b)
        reference = np.array([1, 2, 1, 2])
        # The groups can be "simple" or "scrambled", if they are simple they
        # do not contain repetitions and they are ordered. The reference for
        # the simple groups should follow the same patern and should be ordered
        # and without repetitions.
        if len(a) == len(np.unique(a)):
            reference = np.unique(reference)
        assert_array_equal(subtract_ab.ix, reference)
        subtract_ba = b.subtract(a)
        reference = np.array([7, 6, 5, 7, 6])
        if len(a) == len(np.unique(a)):
            reference = np.unique(reference)
        assert_array_equal(subtract_ba.ix, reference)
        subtract_ad = a.subtract(d)
        assert_equal(subtract_ad, a)
        subtract_ae = a.subtract(e)
        assert_equal(subtract_ae, a)

    def test_difference(self, groups):
        a, b, c, d, e = groups
        difference_ab = a.difference(b)
        assert_array_equal(difference_ab.ix, np.arange(1, 3))

        difference_ba = b.difference(a)
        assert_array_equal(difference_ba.ix, np.arange(5, 8))

        assert_array_equal(a.difference(d).ix, np.arange(1, 5))
        assert_array_equal(a.difference(e).ix, np.arange(1, 5))

    def test_symmetric_difference(self, groups):
        a, b, c, d, e = groups
        symdiff_ab = a.symmetric_difference(b)
        assert_array_equal(symdiff_ab.ix, np.array(list(range(1, 3)) +
                                                   list(range(5, 8))))
        assert_(a.symmetric_difference(b) == b.symmetric_difference(a))
        assert_array_equal(a.symmetric_difference(e).ix, np.arange(1, 8))

    def test_isdisjoint(self, groups):
        a, b, c, d, e = groups
        assert_(a.isdisjoint(e))
        assert_(e.isdisjoint(a))
        assert_(a.isdisjoint(d))
        assert_(d.isdisjoint(a))
        assert_(not a.isdisjoint(b))

    @pytest.mark.parametrize('left, right', itertools.chain(
        # Do inter-levels pairs of groups fail as expected?
        itertools.permutations(component_groups, 2),
        # Do inter-levels pairs of components
        itertools.permutations(components, 2),
        # Do inter-levels pairs of components/groups fail as expected?
        _yield_mix(component_groups, components),
        # Does the function fail with inputs that are not components or groups
        ((u.atoms, 'invalid'), ),
    ))
    def test_failing_pairs(self, left, right):
        def dummy(self, other):
            return True

        with pytest.raises(TypeError):
            mda.core.groups._only_same_level(dummy)(left, right)

    @pytest.mark.parametrize('left, right', itertools.chain(
        # Groups
        _yield_sliced_groups(u, slice(0, 2), slice(1, 3)),
        # Components
        _yield_sliced_groups(u, 0, 1),
        # Mixed
        _yield_sliced_groups(u, slice(0, 2), 1),
        _yield_sliced_groups(u, 1, slice(0, 2)),
    ))
    def test_succeeding_pairs(self, left, right):
        def dummy(self, other):
            return True

        assert_(mda.core.groups._only_same_level(dummy)(left, right))

    def test_only_same_level_different_universes(self):
        def dummy(self, other):
            return True

        u = make_Universe()
        u2 = make_Universe()
        _only_same_level = mda.core.groups._only_same_level
        with pytest.raises(ValueError):
            _only_same_level(dummy)(u.atoms, u2.atoms)

    @pytest.mark.parametrize('op, method', ((operator.add, 'concatenate'),
                                            (operator.sub, 'difference'),
                                            (operator.and_, 'intersection'),
                                            (operator.or_, 'union'),
                                            (operator.xor, 'symmetric_difference')))
    def test_shortcut_overriding(self, op, method, level):
        def check_operator(op, method, level):
            left = getattr(u, level)[1:3]
            right = getattr(u, level)[2:4]
            assert_equal(op(left, right), getattr(left, method)(right))

        n_segments = 5
        n_residues = n_segments * 3
        n_atoms = n_residues * 3
        u = make_Universe(size=(n_atoms, n_residues, n_segments))

        check_operator(op, method, level)


class TestGroupHash(object):
    """
    Groups should be hashable.

    See issue #1397
    """
    levels = ('atoms', 'residues', 'segments')

    @pytest.fixture(params=levels)
    def level(self, request):
        return request.param

    @pytest.fixture(scope='class')
    def u(self):
        return make_Universe(size=(3, 3, 3))

    def test_hash_exists(self, u, level):
        group = getattr(u, level)
        assert isinstance(hash(group), int)

    def test_hash_equality(self, u, level):
        a = getattr(u, level)[0:-1]
        b = getattr(u, level)[0:-1]
        assert hash(a) == hash(b)

    def test_hash_difference(self, u, level):
        a = getattr(u, level)[:-1]
        b = getattr(u, level)[1:]
        assert hash(a) != hash(b)

    @pytest.mark.parametrize('level_a, level_b',
                             itertools.permutations(levels, 2))
    def test_hash_difference_cross(self, u, level_a, level_b):
        a = getattr(u, level_a)[0:-1]
        b = getattr(u, level_b)[0:-1]
        assert hash(a) != hash(b)

    def test_hash_diff_cross_universe(self, level, u):
        u2 = make_Universe(size=(3, 3, 3))
        a = getattr(u, level)
        b = getattr(u2, level)
        assert hash(a) != hash(b)


class TestAtomGroup(object):

    def test_PDB_atom_repr(self):
        u = make_Universe(extras=('altLocs', 'names', 'types', 'resnames', 'resids', 'segids'))
        assert_equal("<Atom 1: AAA of type TypeA of resname RsA, resid 1 and segid SegA and altLoc A>", u.atoms[0].__repr__())


class TestInstantSelectorDeprecationWarnings(object):
    @pytest.fixture()
    def u(self):
        return make_Universe(("resids", "resnames", "segids", "names"))

    def test_AtomGroup_warn_getitem(self, u):
        name = u.atoms[0].name
        assert_warns(DeprecationWarning, lambda x: u.atoms[x], name)

    def test_AtomGroup_nowarn_getitem_index(self, u):
        assert_nowarns(DeprecationWarning, lambda x: u.atoms[x], 0)

    def test_AtomGroup_nowarn_segids_attribute(self, u):
        assert_nowarns(DeprecationWarning, lambda x: getattr(u.atoms, x), "segids")

    def test_AtomGroup_warn_getattr(self, u):
        name = u.atoms[0].name
        assert_warns(DeprecationWarning, lambda x: getattr(u.atoms, x), name)

    def test_ResidueGroup_warn_getattr_resname(self, u):
        name = u.residues[0].resname
        assert_warns(DeprecationWarning, lambda x: getattr(u.residues, x), name)

    def test_Segment_warn_getattr_resname(self, u):
        name = u.residues[0].resname
        assert_warns(DeprecationWarning, lambda x: getattr(u.segments[0], x), name)

    def test_Segment_warn_getattr_rRESNUM(self, u):
        assert_warns(DeprecationWarning, lambda x: getattr(u.segments[0], x), 'r1')

    def test_SegmentGroup_warn_getattr(self, u):
        name = u.segments[0].segid
        assert_warns(DeprecationWarning, lambda x: getattr(u.segments, x), name)

    def test_SegmentGroup_nowarn_getitem(self, u):
        assert_nowarns(DeprecationWarning, lambda x: u.segments[x], 0)

