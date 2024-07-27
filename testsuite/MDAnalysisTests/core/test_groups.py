# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import itertools
import re
import numpy as np
from numpy.testing import (
    assert_array_equal,
    assert_equal,
    assert_almost_equal
)
import pytest
import operator
import warnings

import MDAnalysis as mda
from MDAnalysis.exceptions import NoDataError
from MDAnalysisTests import make_Universe, no_deprecated_call
from MDAnalysisTests.datafiles import PSF, DCD, TPR
from MDAnalysis.core import groups


class TestGroupProperties(object):
    """ Test attributes common to all groups
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

    uni = make_Universe() # can't use fixtures in @pytest.mark.parametrize

    def test_dimensions(self, u, group_dict):
        dimensions = np.arange(6)

        for group in group_dict.values():
            group.dimensions = dimensions.copy()
            assert_array_equal(group.dimensions, dimensions)
            assert_equal(u.dimensions, group.dimensions)
    
    @pytest.mark.parametrize('group', (uni.atoms[:2], uni.residues[:2],
                                       uni.segments[:2]))
    def test_group_isunique(self, group):
        assert len(group) == 2
        # Initially, cache must be empty:
        with pytest.raises(KeyError):
            _ = group._cache['isunique']
        # Check for correct value and type:
        assert group.isunique is True
        # Check if cache is set correctly:
        assert group._cache['isunique'] is True

        # Add duplicate element to group:
        group += group[0]
        assert len(group) == 3
        # Cache must be reset since the group changed:
        with pytest.raises(KeyError):
            _ = group._cache['isunique']
        # Check for correct value and type:
        assert group.isunique is False
        # Check if cache is set correctly:
        assert group._cache['isunique'] is False

        #Check empty group:
        group = group[[]]
        assert len(group) == 0
        # Cache must be empty:
        with pytest.raises(KeyError):
            _ = group._cache['isunique']
        # Check for correct value and type:
        assert group.isunique is True
        # Check if cache is set correctly:
        assert group._cache['isunique'] is True

    @pytest.mark.parametrize('group', (uni.atoms[:2], uni.residues[:2],
                                       uni.segments[:2]))
    def test_group_unique_nocache(self, group):
        # check unique group:
        assert len(group) == 2
        # assert caches are empty:
        for attr in ('isunique', 'sorted_unique', 'unsorted_unique'):
            assert attr not in group._cache

    @pytest.mark.parametrize('group', (uni.atoms[:2], uni.residues[:2],
                                       uni.segments[:2]))
    def test_create_unique_group_from_unique(self, group):
        unique_group = group.asunique(sorted=True)
        # assert identity and caches
        assert unique_group is group
        assert group._cache['sorted_unique'] is unique_group
        assert group._cache['unsorted_unique'] is unique_group
        assert group._cache['isunique'] is True
        assert group._cache['issorted']  # numpy.bool != bool

        # assert .unique copies
        assert group.unique is not group
        assert group.unique == group

        # add duplicate element to group:
        group += group[0]
        assert len(group) == 3

        # assert caches are cleared since the group changed:
        for attr in ('isunique', 'sorted_unique', 'unsorted_unique'):
            assert attr not in group._cache

        # now not unique
        assert group.isunique is False
        assert not group.issorted
        # assert that group.unique of the non-unique group is not the old one
        assert group.asunique() is not unique_group
        assert group.asunique() == unique_group
        # check if caches have been set correctly:
        assert group._cache['unsorted_unique'] is group.asunique()
        assert group._cache['sorted_unique'] is group.asunique()
        assert group.unique is not group.asunique()

        # check length and type:
        assert len(group.unique) == 2
        assert type(group.unique) is type(group)

        # check if caches of group.sorted_unique have been set correctly:
        assert group.sorted_unique._cache['isunique'] is True
        assert group.sorted_unique._cache['sorted_unique'] is group.sorted_unique
        # assert that repeated access yields the same object (not a copy):
        unique_group = group.sorted_unique
        assert unique_group is group.sorted_unique

    @pytest.mark.parametrize('ugroup', [uni.atoms, uni.residues, uni.segments])
    @pytest.mark.parametrize('ix, unique_ix', [
        ([0, 1], [0, 1]),
        ([4, 3, 3, 1], [1, 3, 4])
    ])
    def test_group_unique_returns_sorted_copy(self, ugroup, ix, unique_ix):
        # is copy
        group = ugroup[ix]
        assert group.unique is not group
        # sorted
        assert_equal(group.unique.ix, unique_ix)

    @pytest.mark.parametrize('ugroup', [uni.atoms, uni.residues, uni.segments])
    @pytest.mark.parametrize('ix, value', [
        ([4, 3, 3, 1], False),
        ([1, 3, 4], True),
        ([2, 2, 2, 4], True),
    ])
    def test_group_issorted(self, ugroup, ix, value):
        assert ugroup[ix].issorted == value

    @pytest.mark.parametrize('ugroup', [uni.atoms, uni.residues, uni.segments])
    @pytest.mark.parametrize('ix, sort, unique_ix, is_same', [
        ([1, 3, 4], True, [1, 3, 4], True),
        ([1, 3, 4], False, [1, 3, 4], True),
        ([4, 3, 1], True, [1, 3, 4], False),
        ([4, 3, 1], False, [4, 3, 1], True),
        ([1, 3, 3, 4], True, [1, 3, 4], False),
        ([1, 3, 3, 4], False, [1, 3, 4], False),
        ([4, 3, 3, 1], True, [1, 3, 4], False),
        ([4, 3, 3, 1], False, [4, 3, 1], False),
    ])
    def test_group_asunique(self, ugroup, ix, sort, unique_ix, is_same):
        group = ugroup[ix]
        unique_group = group.asunique(sorted=sort)
        assert_equal(unique_group.ix, unique_ix)
        if is_same:
            assert unique_group is group

    @pytest.mark.parametrize('ugroup', [uni.atoms, uni.residues, uni.segments])
    def test_group_return_sorted_unsorted_unique(self, ugroup):
        unsorted_unique = ugroup[[1, 3, 4]].asunique(sorted=False)
        assert 'unsorted_unique' in unsorted_unique._cache
        assert 'sorted_unique' not in unsorted_unique._cache
        assert 'issorted' not in unsorted_unique._cache
        assert 'isunique' in unsorted_unique._cache

        sorted_unique = unsorted_unique.asunique(sorted=True)
        assert sorted_unique is unsorted_unique
        assert unsorted_unique._cache['issorted']
        assert unsorted_unique._cache['sorted_unique'] is unsorted_unique

    @pytest.mark.parametrize('ugroup', [uni.atoms, uni.residues, uni.segments])
    def test_group_return_unsorted_sorted_unique(self, ugroup):
        unique = ugroup[[1, 3, 3, 4]]
        sorted_unique = unique.asunique(sorted=True)
        assert unique._cache['sorted_unique'] is sorted_unique
        assert 'unsorted_unique' not in unique._cache

        unsorted_unique = unique.asunique(sorted=False)
        assert unsorted_unique is sorted_unique
        assert unique._cache['unsorted_unique'] is sorted_unique


class TestEmptyAtomGroup(object):
    """ Test empty atom groups
    """
    u = mda.Universe(PSF, DCD)

    @pytest.mark.parametrize('ag', [u.residues[:1]])
    def test_passive_decorator(self, ag):
        assert_almost_equal(ag.center_of_mass(), np.array([10.52567673,  9.49548312, -8.15335145]))
        assert_almost_equal(ag.total_mass(), 133.209)
        assert_almost_equal(ag.moment_of_inertia(), np.array([[ 657.514361 ,  104.9446833,  110.4782   ],
                                                              [ 104.9446833,  307.4360346, -199.1794289],
                                                              [ 110.4782   , -199.1794289,  570.2924896]]))
        assert_almost_equal(ag.radius_of_gyration(), 2.400527938286)
        assert_almost_equal(ag.shape_parameter(), 0.61460819)
        assert_almost_equal(ag.asphericity(), 0.4892751412)
        assert_almost_equal(ag.principal_axes(), np.array([[ 0.7574113, -0.113481 ,  0.643001 ],
                                                           [ 0.5896252,  0.5419056, -0.5988993],
                                                           [-0.2804821,  0.8327427,  0.4773566]]))
        assert_almost_equal(ag.center_of_charge(), np.array([11.0800112,  8.8885659, -8.9886632]))
        assert_almost_equal(ag.total_charge(), 1)

    @pytest.mark.parametrize('ag', [mda.AtomGroup([],u)])
    def test_error_empty_group(self, ag):
        with pytest.raises(ValueError, match ="AtomGroup is empty"):
            ag.center_of_mass()
        with pytest.raises(ValueError, match ="AtomGroup is empty"):
            ag.total_mass()
        with pytest.raises(ValueError, match ="AtomGroup is empty"):
            ag.moment_of_inertia()
        with pytest.raises(ValueError, match ="AtomGroup is empty"):
            ag.radius_of_gyration()
        with pytest.raises(ValueError, match ="AtomGroup is empty"):
            ag.shape_parameter()
        with pytest.raises(ValueError, match ="AtomGroup is empty"):
            ag.asphericity()
        with pytest.raises(ValueError, match ="AtomGroup is empty"):
            ag.principal_axes()
        with pytest.raises(ValueError, match ="AtomGroup is empty"):
            ag.center_of_charge()
        with pytest.raises(ValueError, match ="AtomGroup is empty"):
            ag.total_charge()


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
        assert len(group.atoms) == group.n_atoms

    def test_n_residues(self, group):
        assert len(group.residues) == group.n_residues

    def test_n_segments(self, group):
        assert len(group.segments) == group.n_segments

    def test_len(self, group, level):
        ref = self.length[level]
        assert len(group) == ref

    @pytest.mark.parametrize('func', [list, np.array])
    def test_boolean_slicing(self, group, func):
        # func is the container type that will be used to slice
        group = group[:5]
        sli = func([True, False, False, True, True])
        result = group[sli]
        assert len(result) == 3
        for ref, val in zip(sli, group):
            if ref:
                assert val in result
            else:
                assert val not in result

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

        assert len(g2) == len(o2)
        # Check identity of items in the sliced result
        for o, g in zip(o2, g2):
            if o in nparray:
                assert g in g2
            else:
                assert g not in g2

    @pytest.mark.parametrize('idx', [0, 1, -1, -2])
    def test_integer_getitem(self, group, nparray, idx, singular):
        a = group[idx]
        ref = nparray[idx]

        assert a.ix == ref
        assert isinstance(a, singular)
 
    def test_none_getitem(self, group):
        with pytest.raises(TypeError):
            group[None]


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
        assert len(newgroup) == reflen
        assert isinstance(newgroup, refclass)
        # Check ordering of created Group
        for x, y in zip(newgroup, itertools.chain(self.itr(a), self.itr(b))):
            assert x == y

    @pytest.mark.parametrize(
        'a, b, c, refclass',
        _yield_groups(group_dict, singles, levels, groupclasses, repeat=3)
    )
    def test_sum(self, a, b, c, refclass):
        # weird hack in radd allows this
        summed = sum([a, b, c])

        assert isinstance(summed, refclass)
        assert_equal(len(summed),
                     len(self.itr(a)) + len(self.itr(b)) + len(self.itr(c)))
        for x, y in zip(summed,
                        itertools.chain(self.itr(a), self.itr(b), self.itr(c))):
            assert x == y

    @pytest.mark.parametrize(
        'a, b, c, refclass',
        _yield_groups(group_dict, singles, levels, groupclasses, repeat=3)
    )
    def test_bad_sum(self, a, b, c, refclass):
        # sum with bad first argument
        with pytest.raises(TypeError):
            sum([10, a, b, c])

    def test_contains(self, group):
        assert group[2] in group

    def test_contains_false(self, group):
        assert not group[3] in group[:2]

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
        assert not group[2] in group2

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
    """

    @pytest.fixture()
    def u(self):
        return make_Universe()

    def test_atomgroup_to_atomgroup(self, u):
        atm = u.atoms.atoms
        assert len(atm) == 125
        assert isinstance(atm, groups.AtomGroup)
        assert atm is u.atoms

    def test_atomgroup_to_residuegroup(self, u):
        atm = u.atoms
        res = atm.residues
        assert len(res) == 25
        assert isinstance(res, groups.ResidueGroup)
        assert res == u.residues
        assert res is not u.residues
        assert res._cache['isunique'] is True
        assert res._cache['sorted_unique'] is res

    def test_atomgroup_to_segmentgroup(self, u):
        seg = u.atoms.segments
        assert len(seg) == 5
        assert isinstance(seg, groups.SegmentGroup)
        assert seg == u.segments
        assert seg is not u.segments
        assert seg._cache['isunique'] is True
        assert seg._cache['sorted_unique'] is seg

    def test_residuegroup_to_atomgroup(self, u):
        res = u.residues
        atm = res.atoms
        assert len(atm) == 125
        assert isinstance(atm, groups.AtomGroup)
        assert atm == u.atoms
        assert atm is not u.atoms
        # clear res' uniqueness caches:
        if 'sorted_unique' in res._cache.keys():
            del res._cache['sorted_unique']
        if 'isunique' in res._cache.keys():
            del res._cache['isunique']
        atm = res.atoms
        # assert uniqueness caches of atm are empty:
        with pytest.raises(KeyError):
            _ = atm._cache['isunique']
        with pytest.raises(KeyError):
            _ = atm._cache['sorted_unique']
        # populate uniqueness cache of res:
        assert res.isunique
        atm = res.atoms
        # assert uniqueness caches of atm are set:
        assert atm._cache['isunique'] is True
        assert atm._cache['unsorted_unique'] is atm

    def test_residuegroup_to_residuegroup(self, u):
        res = u.residues.residues
        assert len(res) == 25
        assert isinstance(res, groups.ResidueGroup)
        assert res is u.residues

    def test_residuegroup_to_segmentgroup(self, u):
        seg = u.residues.segments
        assert len(seg) == 5
        assert isinstance(seg, groups.SegmentGroup)
        assert seg == u.segments
        assert seg is not u.segments
        assert seg._cache['isunique'] is True
        assert seg._cache['sorted_unique'] is seg

    def test_segmentgroup_to_atomgroup(self, u):
        seg = u.segments
        atm = seg.atoms
        assert len(atm) == 125
        assert isinstance(atm, groups.AtomGroup)
        assert atm == u.atoms
        assert atm is not u.atoms
        # clear seg's uniqueness caches:
        if 'sorted_unique' in seg._cache.keys():
            del seg._cache['sorted_unique']
        if 'isunique' in seg._cache.keys():
            del seg._cache['isunique']
        atm = seg.atoms
        # assert uniqueness caches of atm are empty:
        with pytest.raises(KeyError):
            _ = atm._cache['isunique']
        with pytest.raises(KeyError):
            _ = atm._cache['sorted_unique']
        # populate uniqueness cache of seg:
        assert seg.isunique
        atm = seg.atoms
        # assert uniqueness caches of atm are set:
        assert atm._cache['isunique'] is True
        assert atm._cache['unsorted_unique'] is atm

    def test_segmentgroup_to_residuegroup(self, u):
        seg = u.segments
        res = seg.residues
        assert len(res) == 25
        assert isinstance(res, groups.ResidueGroup)
        assert res == u.residues
        assert res is not u.residues
        # clear seg's uniqueness caches:
        if 'sorted_unique' in seg._cache.keys():
            del seg._cache['sorted_unique']
        if 'isunique' in seg._cache.keys():
            del seg._cache['isunique']
        res = seg.residues
        # assert uniqueness caches of res are empty:
        with pytest.raises(KeyError):
            _ = res._cache['isunique']
        with pytest.raises(KeyError):
            _ = res._cache['sorted_unique']
        # populate uniqueness cache of seg:
        assert seg.isunique
        res = seg.residues
        # assert uniqueness caches of res are set:
        assert res._cache['isunique'] is True
        assert res._cache['unsorted_unique'] is res

    def test_segmentgroup_to_segmentgroup(self, u):
        seg = u.segments.segments
        assert len(seg) == 5
        assert isinstance(seg, groups.SegmentGroup)
        assert seg is u.segments

    def test_atom_to_residue(self, u):
        res = u.atoms[0].residue
        assert isinstance(res, groups.Residue)

    def test_atom_to_segment(self, u):
        seg = u.atoms[0].segment
        assert isinstance(seg, groups.Segment)

    def test_residue_to_atomgroup(self, u):
        ag = u.residues[0].atoms
        assert isinstance(ag, groups.AtomGroup)
        assert len(ag) == 5
        assert ag._cache['isunique'] is True
        assert ag._cache['sorted_unique'] is ag
        del ag._cache['sorted_unique']
        del ag._cache['isunique']
        assert ag.isunique

    def test_residue_to_segment(self, u):
        seg = u.residues[0].segment
        assert isinstance(seg, groups.Segment)

    def test_segment_to_atomgroup(self, u):
        ag = u.segments[0].atoms
        assert isinstance(ag, groups.AtomGroup)
        assert len(ag) == 25
        assert ag._cache['isunique'] is True
        assert ag._cache['sorted_unique'] is ag
        del ag._cache['sorted_unique']
        del ag._cache['isunique']
        assert ag.isunique

    def test_segment_to_residuegroup(self, u):
        rg = u.segments[0].residues
        assert isinstance(rg, groups.ResidueGroup)
        assert len(rg) == 5
        assert rg._cache['isunique'] is True
        assert rg._cache['sorted_unique'] is rg
        del rg._cache['sorted_unique']
        del rg._cache['isunique']
        assert rg.isunique

    def test_atomgroup_to_residuegroup_unique(self, u):
        ag = u.atoms[:5] + u.atoms[10:15] + u.atoms[:5]
        rg = ag.residues
        assert len(rg) == 2
        assert rg._cache['isunique'] is True
        assert rg._cache['sorted_unique'] is rg

    def test_atomgroup_to_segmentgroup_unique(self, u):
        ag = u.atoms[0] + u.atoms[-1] + u.atoms[0]
        sg = ag.segments
        assert len(sg) == 2
        assert sg._cache['isunique'] is True
        assert sg._cache['sorted_unique'] is sg

    def test_residuegroup_to_segmentgroup_unique(self, u):
        rg = u.residues[0] + u.residues[6] + u.residues[1]
        sg = rg.segments
        assert len(sg) == 2
        assert sg._cache['isunique'] is True
        assert sg._cache['sorted_unique'] is sg

    def test_residuegroup_to_atomgroup_listcomp(self, u):
        rg = u.residues[0] + u.residues[0] + u.residues[4]
        ag = rg.atoms
        assert len(ag) == 15
        # assert uniqueness caches of ag are empty:
        with pytest.raises(KeyError):
            _ = ag._cache['isunique']
        with pytest.raises(KeyError):
            _ = ag._cache['sorted_unique']
        # populate uniqueness cache of rg:
        assert not rg.isunique
        ag = rg.atoms
        # ag uniqueness caches are now from residue
        assert not ag._cache['isunique']
        with pytest.raises(KeyError):
            _ = ag._cache['sorted_unique']

    def test_segmentgroup_to_residuegroup_listcomp(self, u):
        sg = u.segments[0] + u.segments[0] + u.segments[1]
        rg = sg.residues
        assert len(rg) == 15
        # assert uniqueness caches of rg are empty:
        with pytest.raises(KeyError):
            _ = rg._cache['isunique']
        with pytest.raises(KeyError):
            _ = rg._cache['sorted_unique']
        # populate uniqueness cache of sg:
        assert not sg.isunique
        rg = sg.residues
        # assert uniqueness caches of rg are now populated
        assert not rg._cache['isunique']
        with pytest.raises(KeyError):
            _ = rg._cache['sorted_unique']

    def test_segmentgroup_to_atomgroup_listcomp(self, u):
        sg = u.segments[0] + u.segments[0] + u.segments[1]
        ag = sg.atoms
        assert len(ag) == 75
        # assert uniqueness caches of ag are empty:
        with pytest.raises(KeyError):
            _ = ag._cache['isunique']
        with pytest.raises(KeyError):
            _ = ag._cache['sorted_unique']
        # populate uniqueness cache of sg:
        assert not sg.isunique
        ag = sg.atoms
        # ag uniqueness caches are now from segment
        assert not ag._cache['isunique']
        with pytest.raises(KeyError):
            _ = ag._cache['sorted_unique']


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
        assert a < b
        assert a < c
        assert not b < a
        assert not a < a

    def test_gt(self, a, b, c):
        assert b > a
        assert c > a
        assert not a > c
        assert not a > a

    def test_ge(self, a, b, c):
        assert b >= a
        assert c >= a
        assert b >= b
        assert not b >= c

    def test_le(self, a, b, c):
        assert b <= c
        assert b <= b
        assert not b <= a

    def test_neq(self, a, b, c):
        assert a != b
        assert not a != a

    def test_eq(self, a, b, c):
        assert a == a
        assert not a == b

    def test_sorting(self, a, b, c):
        assert sorted([b, a, c]) == [a, b, c]

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

        assert isinstance(ng, NewGroup)

        ag = u.atoms[[0, 1, 2]]

        assert_array_equal(ng.positions, ag.positions)


class TestGroupBy(object):
    # tests for the method 'groupby'
    @pytest.fixture()
    def u(self):
        return make_Universe(('segids', 'charges', 'resids'))

    def test_groupby_float(self, u):
        gb = u.atoms.groupby('charges')

        for ref in [-1.5, -0.5, 0.0, 0.5, 1.5]:
            assert ref in gb
            g = gb[ref]
            assert all(g.charges == ref)
            assert len(g) == 25

    @pytest.mark.parametrize('string', ['segids', b'segids', u'segids'])
    def test_groupby_string(self, u, string):
        gb = u.atoms.groupby(string)

        assert len(gb) == 5
        for ref in ['SegA', 'SegB', 'SegC', 'SegD', 'SegE']:
            assert ref in gb
            g = gb[ref]
            assert all(g.segids == ref)
            assert len(g) == 25

    def test_groupby_int(self, u):
        gb = u.atoms.groupby('resids')

        for g in gb.values():
            assert len(g) == 5

    # tests for multiple attributes as arguments

    def test_groupby_float_string(self, u):
        gb = u.atoms.groupby(['charges', 'segids'])

        for ref in [-1.5, -0.5, 0.0, 0.5, 1.5]:
            for subref in ['SegA','SegB','SegC','SegD','SegE']:
                assert (ref, subref) in gb.keys()
                a = gb[(ref, subref)]
                assert len(a) == 5
                assert all(a.charges == ref)
                assert all(a.segids  == subref)

    def test_groupby_int_float(self, u):
        gb = u.atoms.groupby(['resids', 'charges'])

        uplim=int(len(gb)/5+1)
        for ref in range(1, uplim):
            for subref in [-1.5, -0.5, 0.0, 0.5, 1.5]:
                assert (ref, subref) in gb.keys()
                a = gb[(ref, subref)]
                assert len(a) == 1
                assert all(a.resids == ref)
                assert all(a.charges == subref)

    def test_groupby_string_int(self, u):
        gb = u.atoms.groupby(['segids', 'resids'])

        assert len(gb) == 25
        res = 1
        for ref in ['SegA','SegB','SegC','SegD','SegE']:
            for subref in range(0, 5):
                assert (ref, res) in gb.keys()
                a = gb[(ref, res)]
                assert all(a.segids == ref)
                assert all(a.resids == res)
                res += 1


class TestReprs(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(PSF, DCD)

    def test_atom_repr(self, u):
        at = u.atoms[0]
        assert repr(at) == '<Atom 1: N of type 56 of resname MET, resid 1 and segid 4AKE>'

    def test_residue_repr(self, u):
        res = u.residues[0]
        assert repr(res) == '<Residue MET, 1>'

    def test_segment_repr(self, u):
        seg = u.segments[0]
        assert repr(seg) == '<Segment 4AKE>'

    def test_atomgroup_repr(self, u):
        ag = u.atoms[:10]
        assert repr(ag) == '<AtomGroup with 10 atoms>'

    def test_atomgroup_str_short(self, u):
        ag = u.atoms[:2]
        assert str(ag) == '<AtomGroup [<Atom 1: N of type 56 of resname MET, resid 1 and segid 4AKE>, <Atom 2: HT1 of type 2 of resname MET, resid 1 and segid 4AKE>]>'

    def test_atomgroup_str_long(self, u):
        ag = u.atoms[:11]
        assert str(ag).startswith('<AtomGroup [<Atom 1: N of type 56 of resname MET,')
        assert '...' in str(ag)
        assert str(ag).endswith(', resid 1 and segid 4AKE>]>')

    def test_residuegroup_repr(self, u):
        rg = u.residues[:10]
        assert repr(rg) == '<ResidueGroup with 10 residues>'

    def test_residuegroup_str_short(self, u):
        rg = u.residues[:2]
        assert str(rg) == '<ResidueGroup [<Residue MET, 1>, <Residue ARG, 2>]>'

    def test_residuegroup_str_long(self, u):
        rg = u.residues[:11]
        assert str(rg).startswith('<ResidueGroup [<Residue MET, 1>,')
        assert '...' in str(rg)
        assert str(rg).endswith(', <Residue ALA, 11>]>')

    def test_segmentgroup_repr(self, u):
        sg = u.segments[:10]
        assert repr(sg) == '<SegmentGroup with 1 segment>'

    def test_segmentgroup_str(self, u):
        sg = u.segments[:10]
        assert str(sg) == '<SegmentGroup [<Segment 4AKE>]>'


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
        assert a == a
        assert a != b
        assert not a == b
        assert not a[0:1] == a[0], \
            'Element should not equal single element group.'

    @pytest.mark.parametrize('group', (u.atoms[:2], u.residues[:2],
                                       u.segments[:2]))
    def test_copy(self, group):
        # make sure uniqueness caches of group are empty:
        with pytest.raises(KeyError):
            _ = group._cache['isunique']
        with pytest.raises(KeyError):
            _ = group._cache['sorted_unique']
        # make a copy:
        cgroup = group.copy()
        # check if cgroup is an identical copy of group:
        assert type(cgroup) is type(group)
        assert cgroup is not group
        assert cgroup == group
        # check if the copied group's uniqueness caches are empty:
        with pytest.raises(KeyError):
            _ = cgroup._cache['isunique']
        with pytest.raises(KeyError):
            _ = cgroup._cache['sorted_unique']
        # populate group's uniqueness caches:
        assert group.isunique
        # make a copy:
        cgroup = group.copy()
        # check if the copied group's uniqueness caches are set correctly:
        assert cgroup._cache['isunique'] is True
        # assert sorted_unique still doesn't exist
        assert 'sorted_unique' not in cgroup._cache
        # add duplicate element to group:
        group += group[0]
        # populate group's uniqueness caches:
        assert not group.isunique
        # make a copy:
        cgroup = group.copy()
        # check if the copied group's uniqueness caches are set correctly:
        assert cgroup._cache['isunique'] is False
        with pytest.raises(KeyError):
            _ = cgroup._cache['sorted_unique']
        # assert that duplicates are preserved:
        assert cgroup == group

    def test_issubset(self, groups):
        a, b, c, d, e = groups
        assert c.issubset(a)
        assert not c.issubset(e)
        assert not a.issubset(c)
        assert d.issubset(a)
        assert not a.issubset(d)

    def test_is_strict_subset(self, groups):
        a, b, c, d, e = groups
        assert c.is_strict_subset(a)
        assert not c.is_strict_subset(e)
        assert not a.is_strict_subset(a)

    def test_issuperset(self, groups):
        a, b, c, d, e = groups
        assert a.issuperset(c)
        assert not e.issuperset(c)
        assert not c.issuperset(a)
        assert a.issuperset(d)
        assert not d.issuperset(a)

    def test_is_strict_superset(self, groups):
        a, b, c, d, e = groups
        assert a.is_strict_superset(c)
        assert not c.is_strict_superset(e)
        assert not a.is_strict_superset(a)

    def test_concatenate(self, groups):
        a, b, c, d, e = groups
        cat_ab = a.concatenate(b)
        assert cat_ab[:len(a)] == a
        assert cat_ab[len(a):] == b

        cat_ba = b.concatenate(a)
        assert cat_ba[:len(b)] == b
        assert cat_ba[len(b):] == a

        cat_aa = a.concatenate(a)
        assert cat_aa[:len(a)] == a
        assert cat_aa[len(a):] == a

        cat_ad = a.concatenate(d)
        assert cat_ad == a

        cat_da = d.concatenate(a)
        assert cat_da == a

    def test_union(self, groups):
        a, b, c, d, e = groups
        union_ab = a.union(b)
        assert union_ab.ix.tolist() == sorted(union_ab.ix)
        assert list(sorted(set(union_ab.ix))) == list(sorted(union_ab.ix))

        assert a.union(b) == b.union(a)
        assert_array_equal(a.union(a).ix, np.arange(1, 5))
        assert a.union(d), np.arange(1, 5)

    def test_intersection(self, groups):
        a, b, c, d, e = groups
        intersect_ab = a.intersection(b)
        assert_array_equal(intersect_ab.ix, np.arange(3, 5))
        assert a.intersection(b) == b.intersection(a)
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
        assert a.symmetric_difference(b) == b.symmetric_difference(a)
        assert_array_equal(a.symmetric_difference(e).ix, np.arange(1, 8))

    def test_isdisjoint(self, groups):
        a, b, c, d, e = groups
        assert a.isdisjoint(e)
        assert e.isdisjoint(a)
        assert a.isdisjoint(d)
        assert d.isdisjoint(a)
        assert not a.isdisjoint(b)

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

        assert mda.core.groups._only_same_level(dummy)(left, right)

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


@pytest.fixture()
def attr_universe():
    return make_Universe(('names', 'resids', 'segids'))

class TestAttributeSetting(object):
    @pytest.mark.parametrize('groupname', ['atoms', 'residues', 'segments'])
    def test_setting_group_fail(self, attr_universe, groupname):
        group = getattr(attr_universe, groupname)

        with pytest.raises(AttributeError):
            group.this = 'that'

    @pytest.mark.parametrize('groupname', ['atoms', 'residues', 'segments'])
    def test_setting_component_fails(self, attr_universe, groupname):
        component = getattr(attr_universe, groupname)[0]

        with pytest.raises(AttributeError):
            component.this = 'that'

    @pytest.mark.parametrize('attr', ['name', 'resid', 'segid'])
    @pytest.mark.parametrize('groupname', ['atoms', 'residues', 'segments'])
    def test_group_set_singular(self, attr_universe, attr, groupname):
        # this should fail as you can't set the 'name' of a 'ResidueGroup'
        group = getattr(attr_universe, groupname)
        with pytest.raises(AttributeError):
            setattr(group, attr, 24)

    def test_atom_set_name(self, attr_universe):
        attr_universe.atoms[0].name = 'this'
        assert attr_universe.atoms[0].name == 'this'

    def test_atom_set_resid(self, attr_universe):
        with pytest.raises(NotImplementedError):
            attr_universe.atoms[0].resid = 24

    def test_atom_set_segid(self, attr_universe):
        with pytest.raises(NotImplementedError):
            attr_universe.atoms[0].segid = 'this'

    def test_residue_set_name(self, attr_universe):
        with pytest.raises(AttributeError):
            attr_universe.residues[0].name = 'this'

    def test_residue_set_resid(self, attr_universe):
        attr_universe.residues[0].resid = 24
        assert attr_universe.residues[0].resid == 24

    def test_residue_set_segid(self, attr_universe):
        with pytest.raises(NotImplementedError):
            attr_universe.residues[0].segid = 'this'

    def test_segment_set_name(self, attr_universe):
        with pytest.raises(AttributeError):
            attr_universe.segments[0].name = 'this'

    def test_segment_set_resid(self, attr_universe):
        with pytest.raises(AttributeError):
            attr_universe.segments[0].resid = 24

    def test_segment_set_segid(self, attr_universe):
        attr_universe.segments[0].segid = 'this'
        assert attr_universe.segments[0].segid == 'this'

    @pytest.mark.parametrize('attr', ['names', 'resids', 'segids'])
    @pytest.mark.parametrize('groupname', ['atoms', 'residues', 'segments'])
    def test_component_set_plural(self, attr, groupname):
        # this should fail as you can't set the 'Names' of an 'Atom'
        u = make_Universe(('names', 'resids', 'segids'))
        group = getattr(u, groupname)
        comp = group[0]
        with pytest.raises(AttributeError):
            setattr(comp, attr, 24)

class TestAttributeGetting(object):

    @staticmethod
    @pytest.fixture()
    def universe():
        return make_Universe(extras=('masses', 'altLocs'))

    @staticmethod
    @pytest.fixture()
    def atoms():
        u = make_Universe(extras=("masses",), size=(3, 1, 1))
        return u.atoms

    @pytest.mark.parametrize('attr', ['masses', 'altLocs'])
    def test_get_present_topattr_group(self, universe, attr):
        values = getattr(universe.atoms, attr)
        assert values is not None

    @pytest.mark.parametrize('attr', ['mass', 'altLoc'])
    def test_get_present_topattr_component(self, universe, attr):
        value = getattr(universe.atoms[0], attr)
        assert value is not None

    @pytest.mark.parametrize('attr,singular', [
        ('masses', 'mass'),
        ('altLocs', 'altLoc')])
    def test_get_plural_topattr_from_component(self, universe, attr, singular):
        with pytest.raises(AttributeError) as exc:
            getattr(universe.atoms[0], attr)
        assert ('Do you mean ' + singular) in str(exc.value)

    @pytest.mark.parametrize('attr,singular', [
        ('masses', 'mass'),
        ('altLocs', 'altLoc')])
    def test_get_sing_topattr_from_group(self, universe, attr, singular):
        with pytest.raises(AttributeError) as exc:
            getattr(universe.atoms, singular)
        assert ('Do you mean ' + attr) in str(exc.value)

    @pytest.mark.parametrize('attr,singular', [
        ('elements', 'element'),
        ('tempfactors', 'tempfactor'),
        ('bonds', 'bonds')])
    def test_get_absent_topattr_group(self, universe, attr, singular):
        with pytest.raises(NoDataError) as exc:
            getattr(universe.atoms, attr)
        assert 'does not contain ' + singular in str(exc.value)

    def test_get_non_topattr(self, universe):
        with pytest.raises(AttributeError) as exc:
            universe.atoms.jabberwocky
        assert 'has no attribute' in str(exc.value)
        
    def test_unwrap_without_bonds(self, universe):
        with pytest.raises(NoDataError) as exc:
            universe.atoms.unwrap()
        
        expected_message = (
            "AtomGroup.unwrap() not available; this AtomGroup lacks defined bonds. "
            "To resolve this, you can either:\n"
            "1. Guess the bonds at universe creation using `guess_bonds = True`, or\n"
            "2. Create a universe using a topology format where bonds are pre-defined."
        )
        expected_message_pattern = re.escape(expected_message)
        assert re.fullmatch(expected_message_pattern, str(exc.value))

    def test_get_absent_attr_method(self, universe):
        with pytest.raises(NoDataError) as exc:
            universe.atoms.total_charge()
        err = ('AtomGroup.total_charge() not available; '
               'this requires charges')
        assert str(exc.value) == err

    def test_get_absent_attrprop(self, universe):
        with pytest.raises(NoDataError) as exc:
            universe.atoms.fragindices
        err = ('AtomGroup.fragindices not available; '
               'this requires bonds')
        assert str(exc.value) == err

    def test_attrprop_wrong_group(self, universe):
        with pytest.raises(AttributeError) as exc:
            universe.atoms[0].fragindices
        err = ('fragindices is a property of AtomGroup, not Atom')
        assert str(exc.value) == err

    def test_attrmethod_wrong_group(self, universe):
        with pytest.raises(AttributeError) as exc:
            universe.atoms[0].center_of_mass()
        err = ('center_of_mass() is a method of AtomGroup, not Atom')
        assert str(exc.value) == err

    @pytest.mark.parametrize('attr', ['altlocs', 'alt_Locs'])
    def test_wrong_name(self, universe, attr):
        with pytest.raises(AttributeError) as exc:
            getattr(universe.atoms, attr)
        err = ('AtomGroup has no attribute {}. '
               'Did you mean altLocs?').format(attr)
        assert str(exc.value) == err

class TestInitGroup(object):
    @staticmethod
    @pytest.fixture(
        params=['atoms', 'residues', 'segments']
    )
    def components(request):
        # return list of Component and container class for all three levels
        u = make_Universe()

        group = getattr(u, request.param)

        cls = {
            'atoms': mda.AtomGroup,
            'residues': mda.ResidueGroup,
            'segments': mda.SegmentGroup,
        }[request.param]

        yield (u, [group[0], group[2], group[4]], cls)

    def test_object_init(self, components):
        u, objects, cls = components

        group = cls(objects)

        assert len(group) == len(objects)
        for obj in objects:
            assert obj in group

    def test_number_init(self, components):
        u, objects, cls = components

        group = cls([0, 2, 4], u)
        for obj in objects:
            assert obj in group

    def test_VE_no_uni(self, components):
        u, objects, cls = components

        with pytest.raises(TypeError):
            cls([0, 2, 4])  # missing Universe

    def test_VE_no_uni_2(self, components):
        u, objects, cls = components

        with pytest.raises(TypeError):
            cls(0, 2, 4)  # missing Universe


class TestDecorator(object):
    @groups._pbc_to_wrap
    @groups.check_wrap_and_unwrap
    def dummy_funtion(cls, compound="group", wrap=True, unwrap=True):
        return 0

    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('pbc', (True, False))
    @pytest.mark.parametrize('unwrap', (True, False))
    def test_wrap_and_unwrap_deprecation(self, compound, pbc, unwrap):

        if pbc and unwrap:
            with pytest.raises(ValueError):
                # We call a deprecated argument that does not appear in the
                # function's signature. This is done on purpose to test the
                # deprecation. We need to tell the linter.
                # pylint: disable-next=unexpected-keyword-arg
                self.dummy_funtion(compound=compound, pbc=pbc, unwrap=unwrap)
        else:
            with pytest.warns(DeprecationWarning):
                # We call a deprecated argument that does not appear in the
                # function's signature. This is done on purpose to test the
                # deprecation. We need to tell the linter.
                # pylint: disable-next=unexpected-keyword-arg
                assert self.dummy_funtion(compound=compound, pbc=pbc, unwrap=unwrap) == 0

    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('wrap', (True, False))
    @pytest.mark.parametrize('unwrap', (True, False))
    def test_wrap_and_unwrap(self, compound, wrap, unwrap):

        if wrap and unwrap:
            with pytest.raises(ValueError):
                self.dummy_funtion(compound=compound, wrap=wrap, unwrap=unwrap)
        else:
            assert self.dummy_funtion(compound=compound, wrap=wrap, unwrap=unwrap) == 0


@pytest.fixture()
def tpr():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                message="No coordinate reader found")
        return mda.Universe(TPR)

class TestGetConnectionsAtoms(object):
    """Test Atom and AtomGroup.get_connections"""

    @pytest.mark.parametrize("typename",
                             ["bonds", "angles", "dihedrals", "impropers"])
    def test_connection_from_atom_not_outside(self, tpr, typename):
        cxns = tpr.atoms[1].get_connections(typename, outside=False)
        assert len(cxns) == 0

    @pytest.mark.parametrize("typename, n_atoms", [
        ("bonds", 1),
        ("angles", 3),
        ("dihedrals", 4),
    ])
    def test_connection_from_atom_outside(self, tpr, typename, n_atoms):
        cxns = tpr.atoms[10].get_connections(typename, outside=True)
        assert len(cxns) == n_atoms

    @pytest.mark.parametrize("typename, n_atoms", [
        ("bonds", 9),
        ("angles", 15),
        ("dihedrals", 12),
    ])
    def test_connection_from_atoms_not_outside(self, tpr, typename,
                                               n_atoms):
        ag = tpr.atoms[:10]
        cxns = ag.get_connections(typename, outside=False)
        assert len(cxns) == n_atoms
        indices = np.ravel(cxns.to_indices())
        assert np.all(np.isin(indices, ag.indices))

    @pytest.mark.parametrize("typename, n_atoms", [
        ("bonds", 13),
        ("angles", 27),
        ("dihedrals", 38),
    ])
    def test_connection_from_atoms_outside(self, tpr, typename, n_atoms):
        ag = tpr.atoms[:10]
        cxns = ag.get_connections(typename, outside=True)
        assert len(cxns) == n_atoms
        indices = np.ravel(cxns.to_indices())
        assert not np.all(np.isin(indices, ag.indices))

    def test_invalid_connection_error(self, tpr):
        with pytest.raises(AttributeError, match="does not contain"):
            ag = tpr.atoms[:10]
            ag.get_connections("ureybradleys")

    @pytest.mark.parametrize("outside", [True, False])
    def test_get_empty_group(self, tpr, outside):
        imp = tpr.impropers
        ag = tpr.atoms[:10]
        cxns = ag.get_connections("impropers", outside=outside)
        assert len(imp) == 0
        assert len(cxns) == 0


class TestGetConnectionsResidues(object):
    """Test Residue and ResidueGroup.get_connections"""

    @pytest.mark.parametrize("typename, n_atoms", [
        ("bonds", 9),
        ("angles", 14),
        ("dihedrals", 9),
        ("impropers", 0),
    ])
    def test_connection_from_res_not_outside(self, tpr, typename, n_atoms):
        cxns = tpr.residues[10].get_connections(typename, outside=False)
        assert len(cxns) == n_atoms

    @pytest.mark.parametrize("typename, n_atoms", [
        ("bonds", 11),
        ("angles", 22),
        ("dihedrals", 27),
        ("impropers", 0),
    ])
    def test_connection_from_res_outside(self, tpr, typename, n_atoms):
        cxns = tpr.residues[10].get_connections(typename, outside=True)
        assert len(cxns) == n_atoms

    @pytest.mark.parametrize("typename, n_atoms", [
        ("bonds", 157),
        ("angles", 290),
        ("dihedrals", 351),
    ])
    def test_connection_from_residues_not_outside(self, tpr, typename,
                                                  n_atoms):
        ag = tpr.residues[:10]
        cxns = ag.get_connections(typename, outside=False)
        assert len(cxns) == n_atoms
        indices = np.ravel(cxns.to_indices())
        assert np.all(np.isin(indices, ag.atoms.indices))

    @pytest.mark.parametrize("typename, n_atoms", [
        ("bonds", 158),
        ("angles", 294),
        ("dihedrals", 360),
    ])
    def test_connection_from_residues_outside(self, tpr, typename, n_atoms):
        ag = tpr.residues[:10]
        cxns = ag.get_connections(typename, outside=True)
        assert len(cxns) == n_atoms
        indices = np.ravel(cxns.to_indices())
        assert not np.all(np.isin(indices, ag.atoms.indices))

    def test_invalid_connection_error(self, tpr):
        with pytest.raises(AttributeError, match="does not contain"):
            ag = tpr.residues[:10]
            ag.get_connections("ureybradleys")

    @pytest.mark.parametrize("outside", [True, False])
    def test_get_empty_group(self, tpr, outside):
        imp = tpr.impropers
        ag = tpr.residues[:10]
        cxns = ag.get_connections("impropers", outside=outside)
        assert len(imp) == 0
        assert len(cxns) == 0


@pytest.mark.parametrize("typename, n_inside", [
    ("intra_bonds", 9),
    ("intra_angles", 15),
    ("intra_dihedrals", 12),
])
def test_topologygroup_gets_connections_inside(tpr, typename, n_inside):
    ag = tpr.atoms[:10]
    cxns = getattr(ag, typename)
    assert len(cxns) == n_inside
    indices = np.ravel(cxns.to_indices())
    assert np.all(np.isin(indices, ag.indices))


@pytest.mark.parametrize("typename, n_outside", [
    ("bonds", 13),
    ("angles", 27),
    ("dihedrals", 38),
])
def test_topologygroup_gets_connections_outside(tpr, typename, n_outside):
    ag = tpr.atoms[:10]
    cxns = getattr(ag, typename)
    assert len(cxns) == n_outside
