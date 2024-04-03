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
"""Tests for MDAnalysis.core.topologyattrs objects.

"""
import numpy as np

from numpy.testing import (
    assert_equal,
    assert_almost_equal,
)
import pytest
from MDAnalysisTests.datafiles import PSF, DCD, PDB_CHECK_RIGHTHAND_PA, MMTF
from MDAnalysisTests import make_Universe

import MDAnalysis as mda
import MDAnalysis.core.topologyattrs as tpattrs
from MDAnalysis.core._get_readers import get_reader_for
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

    @pytest.fixture()
    def universe(self, top):
        return mda.Universe(top)

    def test_len(self, attr):
        assert len(attr) == len(attr.values)


class TestAtomAttr(TopologyAttrMixin):
    """Test atom-level TopologyAttrs.

    """
    values = np.array([7, 3, 69, 9993, 84, 194, 263, 501, 109, 5873])
    single_value = 567
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
        attr.set_atoms(dg, self.single_value)
        assert_equal(attr.get_atoms(dg),
                     np.array([self.single_value, self.single_value]))

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
                      dtype=object)
    single_value = 'Ca2'
    attrclass = tpattrs.Atomnames

    @pytest.fixture()
    def u(self):
        return mda.Universe(PSF, DCD)

    def test_prev_emptyresidue(self, u):
        # Checking group size rather than doing
        # array comparison on zero-sized arrays (Issue #3535)
        groupsize = len(u.residues[[]]._get_prev_residues_by_resid().atoms)
        assert groupsize == 0
        assert groupsize == len(u.residues[[]].atoms)

    def test_next_emptyresidue(self, u):
        # See above re: checking size for zero-sized arrays
        groupsize = len(u.residues[[]]._get_next_residues_by_resid().atoms)
        assert groupsize == 0
        assert groupsize == len(u.residues[[]].atoms)


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
    single_value = 2
    values = np.array([15.2, 395.6, 0.1, 9.8])
    attrclass = tpattrs.ResidueAttr

    def test_set_residue_VE(self, universe):
        # setting e.g. resname to 2 values should fail with VE
        res = universe.residues[0]
        with pytest.raises(ValueError):
            setattr(res, self.attrclass.singular, self.values[:2])

    def test_get_atoms(self, attr):
        assert_equal(attr.get_atoms(DummyGroup([7, 3, 9])),
                     self.values[[3, 2, 2]])

    def test_get_atom(self, universe):
        attr = getattr(universe.atoms[0], self.attrclass.singular)
        assert_equal(attr, self.values[0])

    def test_get_residues(self, attr):
        assert_equal(attr.get_residues(DummyGroup([1, 2, 1, 3])),
                     self.values[[1, 2, 1, 3]])

    def test_set_residues_singular(self, attr):
        dg = DummyGroup([3, 0, 1])
        attr.set_residues(dg, self.single_value)

        assert_equal(attr.get_residues(dg),
                     np.array([self.single_value]*3, dtype=self.values.dtype))

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


class TestResnames(TestResidueAttr):
    attrclass = tpattrs.Resnames
    single_value = 'xyz'
    values = np.array(['a', 'b', '', 'd'], dtype=object)


class TestICodes(TestResnames):
    attrclass = tpattrs.ICodes


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

    @pytest.fixture()
    def universe_pa(self):
        return mda.Universe(PDB_CHECK_RIGHTHAND_PA)

    def test_principal_axes_handedness(self, universe_pa):
        e_vec = universe_pa.atoms.principal_axes()
        assert_almost_equal(np.dot(np.cross(e_vec[0], e_vec[1]), e_vec[2]), 1.0)

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


class TestRecordTypes(object):
    def test_record_types_default(self):
        u = make_Universe()

        u.add_TopologyAttr('record_type')

        assert u.atoms[0].record_type == 'ATOM'
        assert_equal(u.atoms[:10].record_types, 'ATOM')

    @pytest.fixture()
    def rectype_uni(self):
        # standard 125/25/5 universe
        u = make_Universe()
        u.add_TopologyAttr('record_type')
        # first 25 atoms are ATOM (first 5 residues, first segment)
        # 25 to 50th are HETATM (res 5:10, second segment)
        # all after are ATOM
        u.atoms[:25].record_types = 'ATOM'
        u.atoms[25:50].record_types = 'HETATM'
        u.atoms[50:].record_types = 'ATOM'

        return u

    def test_encoding(self, rectype_uni):
        ag = rectype_uni.atoms[:10]

        ag[0].record_type = 'ATOM'
        ag[1:4].record_types = 'HETATM'

        assert ag[0].record_type == 'ATOM'
        assert ag[1].record_type == 'HETATM'

    def test_residue_record_types(self, rectype_uni):
        rt = rectype_uni.residues.record_types

        assert isinstance(rt, list)
        assert len(rt) == 25

        # check return type explicitly
        # some versions of numpy allow bool to str comparison
        assert not rt[0].dtype == bool
        assert (rt[0] == 'ATOM').all()
        assert (rt[5] == 'HETATM').all()

    def test_segment_record_types(self, rectype_uni):
        rt = rectype_uni.segments.record_types

        assert isinstance(rt, list)
        assert len(rt) == 5

        assert not rt[0].dtype == bool
        assert (rt[0] == 'ATOM').all()
        assert (rt[1] == 'HETATM').all()


def test_static_typing():
    ta = tpattrs.Charges(['1.0', '2.0', '3.0'])

    assert isinstance(ta.values, np.ndarray)
    assert ta.values.dtype == float


def test_static_typing_from_empty():
    u = mda.Universe.empty(3)

    u.add_TopologyAttr('masses', values=['1.0', '2.0', '3.0'])

    assert isinstance(u._topology.masses.values, np.ndarray)
    assert isinstance(u.atoms[0].mass, float)


@pytest.mark.parametrize('level, transplant_name', (
    ('atoms', 'center_of_mass'),
    ('atoms', 'center_of_charge'),
    ('atoms', 'total_charge'),
    ('residues', 'total_charge'),
))
def test_stub_transplant_methods(level, transplant_name):
    u = mda.Universe.empty(n_atoms=2)
    group = getattr(u, level)
    with pytest.raises(NoDataError):
        getattr(group, transplant_name)()


@pytest.mark.parametrize('level, transplant_name', (
    ('universe', 'models'),
    ('atoms', 'n_fragments'),
))
def test_stub_transplant_property(level, transplant_name):
    u = mda.Universe.empty(n_atoms=2)
    group = getattr(u, level)
    with pytest.raises(NoDataError):
        getattr(group, transplant_name)


def test_warn_selection_for_strange_dtype():
    err = "A selection keyword could not be automatically generated"

    with pytest.warns(UserWarning, match=err):
        class Star(tpattrs.TopologyAttr):
            singular = "star"  # turns out test_imports doesn't like emoji
            attrname = "stars"  # :(
            per_object = "atom"
            dtype = dict


class TestDeprecateBFactor:

    MATCH = "use the tempfactor attribute"

    @pytest.fixture()
    def universe(self):
        return mda.Universe(MMTF)

    def test_deprecate_bfactors_get_group(self, universe):
        with pytest.warns(DeprecationWarning, match=self.MATCH):
            universe.atoms.bfactors

    def test_deprecate_bfactors_get_atom(self, universe):
        with pytest.warns(DeprecationWarning, match=self.MATCH):
            assert universe.atoms[0].bfactor == universe.atoms[0].tempfactor

    def test_deprecate_bfactors_set_group(self, universe):
        with pytest.warns(DeprecationWarning, match=self.MATCH):
            universe.atoms[:2].bfactors = [3.14, 10]
        assert universe.atoms.tempfactors[0] == 3.14
        assert universe.atoms.tempfactors[1] == 10

        with pytest.warns(DeprecationWarning, match=self.MATCH):
            assert universe.atoms.bfactors[0] == 3.14
            assert universe.atoms.bfactors[1] == 10

    def test_deprecate_bfactors_set_atom(self, universe):
        with pytest.warns(DeprecationWarning, match=self.MATCH):
            universe.atoms[0].bfactor = 3.14
        assert universe.atoms[0].tempfactor == 3.14
        with pytest.warns(DeprecationWarning, match=self.MATCH):
            assert universe.atoms[0].bfactor == 3.14

    def test_deprecate_bfactor_sel(self, universe):
        with pytest.warns(DeprecationWarning, match=self.MATCH):
            universe.select_atoms("bfactor 3")


class TestStringInterning:
    # try and trip up the string interning we use for string attributes
    @pytest.fixture
    def universe(self):
        u = mda.Universe.empty(n_atoms=10, n_residues=2,
                               atom_resindex=[0]*5 + [1] * 5)
        u.add_TopologyAttr('names', values=['A'] * 10)
        u.add_TopologyAttr('resnames', values=['ResA', 'ResB'])
        u.add_TopologyAttr('segids', values=['SegA'])

        return u

    @pytest.mark.parametrize('newname', ['ResA', 'ResB'])
    def test_add_residue(self, universe, newname):
        newres = universe.add_Residue(resname=newname)

        assert newres.resname == newname

        ag = universe.atoms[2]
        ag.residue = newres

        assert ag.resname == newname

    @pytest.mark.parametrize('newname', ['SegA', 'SegB'])
    def test_add_segment(self, universe, newname):
        newseg = universe.add_Segment(segid=newname)

        assert newseg.segid == newname

        rg = universe.residues[0]
        rg.segment = newseg

        assert rg.atoms[0].segid == newname

    def test_issue3437(self, universe):
        newseg = universe.add_Segment(segid='B')

        ag = universe.residues[0].atoms

        ag.residues.segments = newseg

        assert 'B' in universe.segments.segids

        ag2 = universe.select_atoms('segid B')

        assert len(ag2) == 5
        assert (ag2.ix == ag.ix).all()


class Testcenter_of_charge():

    compounds = ['group', 'segments', 'residues', 'molecules', 'fragments']

    @pytest.fixture
    def u(self):
        """A universe containing two dimers with a finite dipole moment."""
        universe = mda.Universe.empty(n_atoms=4,
                                      n_residues=2,
                                      n_segments=2,
                                      atom_resindex=[0, 0, 1, 1],
                                      residue_segindex=[0, 1])

        universe.add_TopologyAttr("masses", [1, 0, 0, 1])
        universe.add_TopologyAttr("charges", [1, -1, -1, 1])
        universe.add_TopologyAttr("bonds", ((0, 1), (2, 3)))
        universe.add_TopologyAttr("resids", [0, 1])
        universe.add_TopologyAttr("molnums", [0, 1])

        positions = np.array([[0, 0, 0], [0, 1, 0], [2, 1, 0], [2, 2, 0]])

        universe.trajectory = get_reader_for(positions)(positions,
                                                        order='fac',
                                                        n_atoms=4)

        for ts in universe.trajectory:
            ts.dimensions = np.array([1, 2, 3, 90, 90, 90])

        return universe

    @pytest.mark.parametrize('compound', compounds)
    def test_coc(self, u, compound):
        coc = u.atoms.center_of_charge(compound=compound)
        if compound == "group":
            coc_ref = [1, 1, 0]
        else:
            coc_ref = [[0, 0.5, 0], [2, 1.5, 0]]
        assert_equal(coc, coc_ref)

    @pytest.mark.parametrize('compound', compounds)
    def test_coc_wrap(self, u, compound):
        coc = u.atoms[:2].center_of_charge(compound=compound, wrap=True)
        assert_equal(coc.flatten(), [0, 0.5, 0])

    @pytest.mark.parametrize('compound', compounds)
    def test_coc_unwrap(self, u, compound):
        u.atoms.wrap
        coc = u.atoms[:2].center_of_charge(compound=compound, unwrap=True)
        assert_equal(coc.flatten(), [0, -0.5, 0])
