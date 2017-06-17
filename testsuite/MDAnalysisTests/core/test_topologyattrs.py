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
"""Tests for MDAnalysis.core.topologyattrs objects.

"""
from __future__ import division, absolute_import
import numpy as np

from numpy.testing import (
    assert_,
    assert_array_equal,
    assert_array_almost_equal,
    assert_raises,
)
from nose.tools import raises
from MDAnalysisTests.plugins.knownfailure import knownfailure
from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysisTests import make_Universe

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
    Ridx = np.array([0, 0, 2, 2, 1, 1, 3, 3, 1, 2])
    Sidx = np.array([0, 1, 1, 0])

    def setUp(self):
        self.top = Topology(10, 4, 2,
                            attrs=[self.attrclass(self.values.copy())],
                            atom_resindex=self.Ridx,
                            residue_segindex=self.Sidx)
        self.attr = getattr(self.top, self.attrclass.attrname)

    def tearDown(self):
        del self.top

    def test_len(self):
        assert len(self.attr) == len(self.attr.values)


class TestAtomAttr(TopologyAttrMixin):
    """Test atom-level TopologyAttrs.

    """
    values = np.array([7, 3, 69, 9993, 84, 194, 263, 501, 109, 5873])
    attrclass = tpattrs.AtomAttr

    @staticmethod
    def test_set_atom_VE():
        u = make_Universe(('names',))
        at = u.atoms[0]

        assert_raises(ValueError, setattr, at, 'name', ['oopsy', 'daisies'])

    def test_get_atoms(self):
        result = self.attr.get_atoms(DummyGroup([2, 1]))

        assert_(len(result) == 2)
        assert_array_equal(result,
                           self.values[[2, 1]])

    def test_set_atoms_singular(self):
        # set len 2 Group to len 1 value
        dg = DummyGroup([3, 7])
        self.attr.set_atoms(dg, 567)
        assert_array_equal(self.attr.get_atoms(dg), np.array([567, 567]))

    def test_set_atoms_plural(self):
        # set len 2 Group to len 2 values
        dg = DummyGroup([3, 7])
        self.attr.set_atoms(dg, np.array([23, 504]))
        assert_array_equal(self.attr.get_atoms(dg), np.array([23, 504]))

    def test_set_atoms_VE(self):
        # set len 2 Group to wrong length values
        dg = DummyGroup([3, 7])
        assert_raises(ValueError, self.attr.set_atoms, dg, np.array([6, 7, 8, 9]))

    def test_get_residues(self):
        """Unless overriden by child class, this should yield values for all
        atoms in residues.

        """
        result = self.attr.get_residues(DummyGroup([2, 1]))

        assert_(len(result) == 2)
        assert_array_equal(result,
                           [self.values[[2, 3, 9]], self.values[[4, 5, 8]]])

    def test_get_segments(self):
        """Unless overriden by child class, this should yield values for all
        atoms in segments.

        """
        result = self.attr.get_segments(DummyGroup([1]))

        assert_(len(result) == 1)
        assert_array_equal(result,
                           [self.values[[4, 5, 8, 2, 3, 9]]])


class TestAtomids(TestAtomAttr):
    attrclass = tpattrs.Atomids


class TestIndicesClasses(object):
    def setUp(self):
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    @raises(AttributeError)
    def test_cant_set_atom_indices(self):
        self.u.atoms.indices = 1

    @raises(AttributeError)
    def test_cant_set_residue_indices(self):
        self.u.atoms.residues.resindices = 1

    @raises(AttributeError)
    def test_cant_set_segment_indices(self):
        self.u.atoms.segments.segindices = 1


class TestAtomnames(TestAtomAttr):
    values = np.array(['O', 'C', 'CA', 'N', 'CB', 'CG', 'CD', 'NA', 'CL', 'OW'],
                      dtype=np.object)
    attrclass = tpattrs.Atomnames


class AggregationMixin(TestAtomAttr):
    def test_get_residues(self):
        assert_array_equal(self.attr.get_residues(DummyGroup([2, 1])),
                           np.array([self.values[[2, 3, 9]].sum(),
                                     self.values[[4, 5, 8]].sum()]))

    def test_get_segments(self):
        assert_array_equal(self.attr.get_segments(DummyGroup([1])),
                           np.array([self.values[[4, 5, 8, 2, 3, 9]].sum()]))

    def test_get_segment(self):
        assert_array_equal(self.attr.get_segments(DummyGroup(1)),
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

    @staticmethod
    def test_set_residue_VE():
        u = make_Universe(('resnames',))
        res = u.residues[0]

        assert_raises(ValueError, setattr, res, 'resname', ['wrong', 'length'])

    def test_get_atoms(self):
        assert_array_equal(self.attr.get_atoms(DummyGroup([7, 3, 9])),
                           self.values[[3, 2, 2]])

    def test_get_residues(self):
        assert_array_equal(self.attr.get_residues(DummyGroup([1, 2, 1, 3])),
                           self.values[[1, 2, 1, 3]])

    def test_set_residues_singular(self):
        dg = DummyGroup([3, 0, 1])
        self.attr.set_residues(dg, 2)

        assert_array_almost_equal(self.attr.get_residues(dg),
                                  np.array([2, 2, 2]))

    def test_set_residues_plural(self):
        self.attr.set_residues(DummyGroup([3, 0, 1]),
                               np.array([23, 504, 2]))
        assert_array_almost_equal(self.attr.get_residues(DummyGroup([3, 0, 1])),
                                  np.array([23, 504, 2]))

    def test_set_residues_VE(self):
        dg = DummyGroup([3, 0, 1])

        assert_raises(ValueError, self.attr.set_residues, dg, np.array([4.5, 5.2]))

    def test_get_segments(self):
        """Unless overriden by child class, this should yield values for all
        atoms in segments.

        """
        assert_array_equal(self.attr.get_segments(DummyGroup([0, 1, 1])),
                           [self.values[[0, 3]], self.values[[1, 2]], self.values[[1, 2]]])


class TestResids(TestResidueAttr):
    values = np.array([10, 11, 18, 20])
    attrclass = tpattrs.Resids

    @knownfailure
    def test_set_atoms(self):
        """Setting the resids of atoms changes their residue membership.

        """
        # moving resids doesn't currently work!
        assert_(1 == 2)

        # set with array
        self.attr.set_atoms(DummyGroup([3, 7]), np.array([11, 20]))
        assert_array_equal(self.attr.get_atoms(DummyGroup([3, 7])), np.array([11, 20]))

        # set to resid that no residue has (should raise exception)
        assert_raises(NoDataError, self.attr.set_atoms, DummyGroup([3, 7]), np.array([11, 21]))

    def test_set_residues(self):
        self.attr.set_residues(DummyGroup([3, 0, 1]),
                               np.array([23, 504, 27]))
        assert_array_almost_equal(self.attr.get_residues(DummyGroup([3, 0, 1])),
                                  np.array([23, 504, 27]))


class TestResnames(TestResidueAttr):
    values = np.array(['VAL', 'LYS', 'VAL', 'POPG'], dtype=np.object)
    attrclass = tpattrs.Resnames

    def test_residuegroup_getattr_single(self):
        u = make_Universe(('resnames',))

        res = u.residues.RsB

        assert_(isinstance(res, groups.Residue))
        assert_(res == u.residues[1])

    def test_residuegroup_getattr_multiple(self):
        u = make_Universe(('resnames',))
        u.residues[:10].resnames = 'ABC'

        rg = u.residues.ABC

        assert_(isinstance(rg, groups.ResidueGroup))
        assert_(len(rg) == 10)

    def test_residuegroup_getattr_AE(self):
        u = make_Universe(('resnames',))

        assert_raises(AttributeError, getattr, u.residues, 'foo')

    def test_segment_getattr_singular(self):
        u = make_Universe(('resnames',))

        res = u.segments[0].RsB

        assert_(isinstance(res, groups.Residue))
        assert_(res == u.residues[1])

    def test_segment_getattr_multiple(self):
        u = make_Universe(('resnames',))
        u.residues[:3].resnames = 'bar'

        rg = u.segments[0].bar

        assert_(isinstance(rg, groups.ResidueGroup))
        assert_(len(rg) == 3)

    def test_segment_getattr_AE(self):
        u = make_Universe(('resnames',))

        assert_raises(AttributeError, getattr, u.segments[0], 'foo')


class TestSegmentAttr(TopologyAttrMixin):
    """Test segment-level TopologyAttrs.

    """
    values = np.array([-0.19, 500])
    attrclass = tpattrs.SegmentAttr

    @staticmethod
    def test_set_segment_VE():
        u = make_Universe(('segids',))
        seg = u.segments[0]

        assert_raises(ValueError, setattr, seg, 'segid', [1, 2, 3])

    def test_get_atoms(self):
        assert_array_equal(self.attr.get_atoms(DummyGroup([2, 4, 1])),
                           self.values[[1, 1, 0]])

    def test_get_residues(self):
        assert_array_equal(self.attr.get_residues(DummyGroup([1, 2, 1, 3])),
                           self.values[[1, 1, 1, 0]])

    def test_get_segments(self):
        """Unless overriden by child class, this should yield values for all
        atoms in segments.

        """
        assert_array_equal(self.attr.get_segments(DummyGroup([1, 0, 0])),
                           self.values[[1, 0, 0]])

    def test_set_segments_singular(self):
        dg = DummyGroup([0, 1])
        self.attr.set_segments(dg, 0.45)
        assert_array_equal(self.attr.get_segments(dg), np.array([0.45, 0.45]))

    def test_set_segments_plural(self):
        dg = DummyGroup([0, 1])
        self.attr.set_segments(dg, np.array([23, -0.0002]))
        assert_array_equal(self.attr.get_segments(dg), np.array([23, -0.0002]))

    def test_set_segments_VE(self):
        dg = DummyGroup([0, 1])
        assert_raises(ValueError, self.attr.set_segments, dg, np.array([4, 5, 6, 7]))


class TestAttr(object):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.ag = self.universe.atoms  # prototypical AtomGroup

    def tearDown(self):
        del self.universe
        del self.ag

    def test_principal_axes(self):
        assert_array_almost_equal(
            self.ag.principal_axes(),
            np.array([
                      [1.53389276e-03, 4.41386224e-02, 9.99024239e-01],
                      [1.20986911e-02, 9.98951474e-01, -4.41539838e-02],
                      [-9.99925632e-01, 1.21546132e-02, 9.98264877e-04]]))

    def test_align_principal_axes_with_self(self):
        pa = self.ag.principal_axes()
        self.ag.align_principal_axis(0, pa[0])


        assert_array_almost_equal(self.ag.principal_axes(), pa)

    def test_align_principal_axes_with_x(self):
        self.ag.align_principal_axis(0, [1, 0, 0])
        # This is a very loose check that the difference is not more then 0.5.
        # This is OK here because the rounding error in the calculation really
        # is that big.
        assert_(np.allclose(np.abs(self.ag.principal_axes()), np.eye(3),
                            rtol=0, atol=0.1))


class TestCrossLevelAttributeSetting(object):
    """

    Can only get attributes belonging to higher level objects

    Atom.resid works!
    ResidueGroup.names = ['a', 'b', 'c'] doesn't work, Atom is below Residue

    Setting any attribute we can get should only work if they are the same level.

    Atom.resid = 4  should fail because resid belongs to Residue not Atom
    """
    @staticmethod
    def _check_crosslevel_fail(item, attr):
        assert_raises(NotImplementedError, setattr, item, attr, 1.0)

    def test_set_crosslevel(self):
        u = make_Universe(('names', 'resids', 'segids'))

        # component and group in each level
        atomlevel = (u.atoms[0], u.atoms[:10])
        residuelevel = (u.residues[0], u.residues[:5])
        segmentlevel = (u.segments[0], u.segments[:2])
        levels = {0:atomlevel, 1:residuelevel, 2:segmentlevel}

        atomattr = 'names'
        residueattr = 'resids'
        segmentattr = 'segids'
        attrs = {0:atomattr, 1:residueattr, 2:segmentattr}

        # loop over Atom, Residue, Segment level
        for level_idx, level in levels.items():
            # loop over an Attribute native to each level
            for attr_idx, attr in attrs.items():
                if level_idx == attr_idx:
                    # if we're on the same level, then this should work
                    # ie Atom.mass = 12.0 is OK!
                    continue
                component, group = level
                # eg 'name', 'names' = 'names', 'names'[:-1]
                singular_attr, plural_attr = attr[:-1], attr

                # eg check ResidueGroup.names = 'newvalue' raises NIE
                # or ResidueGroup.segids = 'newvalue' raises NIE
                yield self._check_crosslevel_fail, group, plural_attr

                if attr_idx < level_idx:
                    # Segment.resid doesn't even exist as an attribute
                    # so we don't have to check that setting fails
                    # Atom.segid does exist as attribute,
                    # but will fail to be set
                    continue
                yield self._check_crosslevel_fail, component, singular_attr
