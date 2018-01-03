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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""Tests for MDAnalysis.core.topologyattrs objects.

"""
from __future__ import division, absolute_import

import numpy as np

from numpy.testing import (
    assert_equal,
    assert_almost_equal,
)
import pytest
from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysisTests import make_Universe, no_deprecated_call

import MDAnalysis as mda
import MDAnalysis.core.topologyattrs as tpattrs
from MDAnalysis.core import groups
from MDAnalysis.core.topology import Topology
from MDAnalysis.exceptions import NoDataError


class DummyGroup(object):
    """Designed to mock an Group

    initiate with indices, these are then available as ._ix
    """
    def __init__(self, vals):
        self._ix = vals

    def __len__(self):
        return len(self._ix)

    @property
    def ix(self):
        return self._ix


class TopologyAttrMixin(object):
    """Mixin to test the common elements to all TopologyAttrs.

    10 atoms
    4 residues
    2 segments

    """
    # Reference data


    @pytest.fixture()
    def top(self):
        Ridx = np.array([0, 0, 2, 2, 1, 1, 3, 3, 1, 2])
        Sidx = np.array([0, 1, 1, 0])
        return Topology(10, 4, 2,
                        attrs=[self.attrclass(self.values.copy())],
                        atom_resindex=Ridx,
                        residue_segindex=Sidx)

    @pytest.fixture()
    def attr(self, top):
        return getattr(top, self.attrclass.attrname)

    def test_len(self, attr):
        assert len(attr) == len(attr.values)


class TestAtomAttr(TopologyAttrMixin):
    """Test atom-level TopologyAttrs.

    """
    values = np.array([7, 3, 69, 9993, 84, 194, 263, 501, 109, 5873])
    attrclass = tpattrs.AtomAttr

    def test_set_atom_VE(self):
        u = make_Universe(('names',))
        at = u.atoms[0]

        with pytest.raises(ValueError):
            setattr(at, 'name', ['oopsy', 'daisies'])

    def test_get_atoms(self, attr):
        result = attr.get_atoms(DummyGroup([2, 1]))

        assert len(result) == 2
        assert_equal(result,
                           self.values[[2, 1]])

    def test_set_atoms_singular(self, attr):
        # set len 2 Group to len 1 value
        dg = DummyGroup([3, 7])
        attr.set_atoms(dg, 567)
        assert_equal(attr.get_atoms(dg), np.array([567, 567]))

    def test_set_atoms_plural(self, attr):
        # set len 2 Group to len 2 values
        dg = DummyGroup([3, 7])
        attr.set_atoms(dg, np.array([23, 504]))
        assert_equal(attr.get_atoms(dg), np.array([23, 504]))

    def test_set_atoms_VE(self, attr):
        # set len 2 Group to wrong length values
        dg = DummyGroup([3, 7])
        with pytest.raises(ValueError):
            attr.set_atoms(dg, np.array([6, 7, 8, 9]))

    def test_get_residues(self, attr):
        """Unless overriden by child class, this should yield values for all
        atoms in residues.

        """
        result = attr.get_residues(DummyGroup([2, 1]))

        assert len(result) == 2
        assert_equal(result,
                           [self.values[[2, 3, 9]], self.values[[4, 5, 8]]])

    def test_get_segments(self, attr):
        """Unless overriden by child class, this should yield values for all
        atoms in segments.

        """
        result = attr.get_segments(DummyGroup([1]))

        assert len(result) == 1
        assert_equal(result,
                           [self.values[[4, 5, 8, 2, 3, 9]]])


class TestAtomids(TestAtomAttr):
    attrclass = tpattrs.Atomids


class TestIndicesClasses(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(PSF, DCD)

    def test_cant_set_atom_indices(self, u):
        with pytest.raises(AttributeError):
            u.atoms.indices = 1

    def test_cant_set_residue_indices(self, u):
        with pytest.raises(AttributeError):
            u.atoms.residues.resindices = 1

    def test_cant_set_segment_indices(self, u):
        with pytest.raises(AttributeError):
            u.atoms.segments.segindices = 1


class TestAtomnames(TestAtomAttr):
    values = np.array(['O', 'C', 'CA', 'N', 'CB', 'CG', 'CD', 'NA', 'CL', 'OW'],
                      dtype=np.object)
    attrclass = tpattrs.Atomnames


class AggregationMixin(TestAtomAttr):
    def test_get_residues(self, attr):
        assert_equal(attr.get_residues(DummyGroup([2, 1])),
                           np.array([self.values[[2, 3, 9]].sum(),
                                     self.values[[4, 5, 8]].sum()]))

    def test_get_segments(self, attr):
        assert_equal(attr.get_segments(DummyGroup([1])),
                           np.array([self.values[[4, 5, 8, 2, 3, 9]].sum()]))

    def test_get_segment(self, attr):
        assert_equal(attr.get_segments(DummyGroup(1)),
                           np.sum(self.values[[4, 5, 8, 2, 3, 9]]))


class TestMasses(AggregationMixin):
    attrclass = tpattrs.Masses


class TestCharges(AggregationMixin):
    values = np.array([+2, -1, 0, -1, +1, +2, 0, 0, 0, -1])
    attrclass = tpattrs.Charges


class TestResidueAttr(TopologyAttrMixin):
    """Test residue-level TopologyAttrs.

    """
    values = np.array([15.2, 395.6, 0.1, 9.8])
    attrclass = tpattrs.ResidueAttr

    def test_set_residue_VE(self):
        u = make_Universe(('resnames',))
        res = u.residues[0]
        with pytest.raises(ValueError):
            setattr(res, 'resname', ['wrong', 'length'])

    def test_get_atoms(self, attr):
        assert_equal(attr.get_atoms(DummyGroup([7, 3, 9])),
                           self.values[[3, 2, 2]])

    def test_get_residues(self, attr):
        assert_equal(attr.get_residues(DummyGroup([1, 2, 1, 3])),
                           self.values[[1, 2, 1, 3]])

    def test_set_residues_singular(self, attr):
        dg = DummyGroup([3, 0, 1])
        attr.set_residues(dg, 2)

        assert_almost_equal(attr.get_residues(dg),
                                  np.array([2, 2, 2]))

    def test_set_residues_plural(self, attr):
        attr.set_residues(DummyGroup([3, 0, 1]),
                               np.array([23, 504, 2]))
        assert_almost_equal(attr.get_residues(DummyGroup([3, 0, 1])),
                                  np.array([23, 504, 2]))

    def test_set_residues_VE(self, attr):
        dg = DummyGroup([3, 0, 1])

        with pytest.raises(ValueError):
            attr.set_residues(dg, np.array([4.5, 5.2]))

    def test_get_segments(self, attr):
        """Unless overriden by child class, this should yield values for all
        atoms in segments.

        """
        assert_equal(attr.get_segments(DummyGroup([0, 1, 1])),
                           [self.values[[0, 3]], self.values[[1, 2]], self.values[[1, 2]]])


class TestResids(TestResidueAttr):
    values = np.array([10, 11, 18, 20])
    attrclass = tpattrs.Resids

    @pytest.mark.xfail
    def test_set_atoms(self, attr):
        """Setting the resids of atoms changes their residue membership.

        """
        # moving resids doesn't currently work!
        assert 1 == 2

        # set with array
        attr.set_atoms(DummyGroup([3, 7]), np.array([11, 20]))
        assert_equal(attr.get_atoms(DummyGroup([3, 7])), np.array([11, 20]))

        # set to resid that no residue has (should raise exception)
        with pytest.raises(NoDataError):
            attr.set_atoms(DummyGroup([3, 7]), np.array([11, 21]))

    def test_set_residues(self, attr):
        attr.set_residues(DummyGroup([3, 0, 1]),
                               np.array([23, 504, 27]))
        assert_almost_equal(attr.get_residues(DummyGroup([3, 0, 1])),
                                  np.array([23, 504, 27]))


class TestResnames(TestResidueAttr):
    values = np.array(['VAL', 'LYS', 'VAL', 'POPG'], dtype=np.object)
    attrclass = tpattrs.Resnames

    def test_residuegroup_getattr_single(self):
        u = make_Universe(('resnames',))

        res = u.residues.RsB

        assert isinstance(res, groups.Residue)
        assert res == u.residues[1]

    def test_residuegroup_getattr_multiple(self):
        u = make_Universe(('resnames',))
        u.residues[:10].resnames = 'ABC'

        rg = u.residues.ABC

        assert isinstance(rg, groups.ResidueGroup)
        assert len(rg) == 10

    def test_residuegroup_getattr_AE(self):
        u = make_Universe(('resnames',))
        with pytest.raises(AttributeError):
            getattr(u.residues, 'foo')


    def test_segment_getattr_singular(self):
        u = make_Universe(('resnames',))

        res = u.segments[0].RsB

        assert isinstance(res, groups.Residue)
        assert res == u.residues[1]

    def test_segment_getattr_multiple(self):
        u = make_Universe(('resnames',))
        u.residues[:3].resnames = 'bar'

        rg = u.segments[0].bar

        assert isinstance(rg, groups.ResidueGroup)
        assert len(rg) == 3

    def test_segment_getattr_AE(self):
        u = make_Universe(('resnames',))
        with pytest.raises(AttributeError):
            getattr(u.segments[0], 'foo')


class TestSegmentAttr(TopologyAttrMixin):
    """Test segment-level TopologyAttrs.

    """
    values = np.array([-0.19, 500])
    attrclass = tpattrs.SegmentAttr

    def test_set_segment_VE(self):
        u = make_Universe(('segids',))
        seg = u.segments[0]
        with pytest.raises(ValueError):
            setattr(seg, 'segid', [1, 2, 3])

    def test_get_atoms(self, attr):
        assert_equal(attr.get_atoms(DummyGroup([2, 4, 1])),
                           self.values[[1, 1, 0]])

    def test_get_residues(self, attr):
        assert_equal(attr.get_residues(DummyGroup([1, 2, 1, 3])),
                           self.values[[1, 1, 1, 0]])

    def test_get_segments(self, attr):
        """Unless overriden by child class, this should yield values for all
        atoms in segments.

        """
        assert_equal(attr.get_segments(DummyGroup([1, 0, 0])),
                           self.values[[1, 0, 0]])

    def test_set_segments_singular(self, attr):
        dg = DummyGroup([0, 1])
        attr.set_segments(dg, 0.45)
        assert_equal(attr.get_segments(dg), np.array([0.45, 0.45]))

    def test_set_segments_plural(self, attr):
        dg = DummyGroup([0, 1])
        attr.set_segments(dg, np.array([23, -0.0002]))
        assert_equal(attr.get_segments(dg), np.array([23, -0.0002]))

    def test_set_segments_VE(self, attr):
        dg = DummyGroup([0, 1])
        with pytest.raises(ValueError):
            attr.set_segments(dg, np.array([4, 5, 6, 7]))


class TestAttr(object):
    @pytest.fixture()
    def ag(self):
        universe = mda.Universe(PSF, DCD)
        return universe.atoms  # prototypical AtomGroup

    def test_principal_axes(self, ag):
        assert_almost_equal(
            ag.principal_axes(),
            np.array([
                      [1.53389276e-03, 4.41386224e-02, 9.99024239e-01],
                      [1.20986911e-02, 9.98951474e-01, -4.41539838e-02],
                      [-9.99925632e-01, 1.21546132e-02, 9.98264877e-04]]))

    def test_align_principal_axes_with_self(self, ag):
        pa = ag.principal_axes()
        ag.align_principal_axis(0, pa[0])

        assert_almost_equal(ag.principal_axes(), pa)

    def test_align_principal_axes_with_x(self, ag):
        ag.align_principal_axis(0, [1, 0, 0])
        # This is a very loose check that the difference is not more then 0.5.
        # This is OK here because the rounding error in the calculation really
        # is that big.
        assert_almost_equal(np.abs(ag.principal_axes()), np.eye(3), decimal=1)


class TestCrossLevelAttributeSetting(object):
    """

    Can only get attributes belonging to higher level objects

    Atom.resid works!
    ResidueGroup.names = ['a', 'b', 'c'] doesn't work, Atom is below Residue

    Setting any attribute we can get should only work if they are the same level.

    Atom.resid = 4  should fail because resid belongs to Residue not Atom
    """

    u = make_Universe(('names', 'resids', 'segids'))

    # component and group in each level
    atomlevel = (u.atoms[0], u.atoms[:10])
    residuelevel = (u.residues[0], u.residues[:5])
    segmentlevel = (u.segments[0], u.segments[:2])
    levels = {0: atomlevel, 1: residuelevel, 2: segmentlevel}

    atomattr = 'names'
    residueattr = 'resids'
    segmentattr = 'segids'
    attrs = {0: atomattr, 1: residueattr, 2: segmentattr}

    @pytest.mark.parametrize('level_idx, level', levels.items())
    @pytest.mark.parametrize('attr_idx, attr', attrs.items())
    def test_set_crosslevel(self, level_idx, level, attr_idx, attr):
        if level_idx == attr_idx:
            # if we're on the same level, then this should work
            # ie Atom.mass = 12.0 is OK!
            return
        component, group = level
        # eg 'name', 'names' = 'names', 'names'[:-1]
        singular_attr, plural_attr = attr[:-1], attr

        # eg check ResidueGroup.names = 'newvalue' raises NIE
        # or ResidueGroup.segids = 'newvalue' raises NIE
        self._check_crosslevel_fail(group, plural_attr)

        if attr_idx < level_idx:
            # Segment.resid doesn't even exist as an attribute
            # so we don't have to check that setting fails
            # Atom.segid does exist as attribute,
            # but will fail to be set
            return
        self._check_crosslevel_fail(component, singular_attr)

    @staticmethod
    def _check_crosslevel_fail(item, attr):
        with pytest.raises(NotImplementedError):
            setattr(item, attr, 1.0)


class TestInstantSelectorDeprecation(object):
    """Test the deprecation warnings for instant selectors

    Instant selectors are deprecated since version 0.16.2. PR #1403 introduced
    deprecation warnings for these selectors.
    """
    @staticmethod
    @pytest.fixture()
    def universe():
        return mda.Universe(PSF, DCD)

    @pytest.mark.parametrize('instruction', (
        'universe.atoms.CA',
        'universe.residues.LYS',
        'universe.segments.s4AKE',
        'universe.s4AKE',
    ))
    def test_deprecation(self, universe, instruction):
        """Test that the warnings are issued when required.
        """
        with pytest.deprecated_call():
            exec(instruction)  #pylint: disable=W0122

    @pytest.mark.parametrize('instruction', (
        'universe.atoms',
        'universe.residues',
        'universe.segments',
    ))
    def test_no_deprecation(self, universe, instruction):
        """Test that the warnings are not issued when they should not.

        See issue #1476.
        """
        with no_deprecated_call():
            exec(instruction)  #pylint: disable=W0122
