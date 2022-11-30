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
import numpy as np
from numpy.testing import assert_equal, assert_almost_equal

import MDAnalysis as mda
from MDAnalysis.exceptions import DuplicateWarning, NoDataError
from MDAnalysisTests.datafiles import PSF, DCD, GRO
from MDAnalysisTests.core.util import UnWrapUniverse
import pytest

levels = ('atoms', 'residues', 'segments')

class TestAccumulate(object):
    """Tests the functionality of *Group.accumulate()."""
    @pytest.fixture(params=levels)
    def group(self, request):
        u = mda.Universe(PSF, DCD)
        return getattr(u, request.param)

    def test_accumulate_str_attribute(self, group):
        assert_almost_equal(group.accumulate("masses"), np.sum(group.atoms.masses))

    def test_accumulate_different_func(self, group):
        assert_almost_equal(
            group.accumulate("masses", function=np.prod),
            np.prod(group.atoms.masses))

    @pytest.mark.parametrize('name, compound', (('resindices', 'residues'),
                                                ('segindices', 'segments'),
                                                ('molnums', 'molecules'),
                                                ('fragindices', 'fragments')))
    @pytest.mark.parametrize('level', levels)
    def test_accumulate_str_attribute_compounds(self, name, compound, level):
        u = UnWrapUniverse()
        group = getattr(u, level)
        ref = [sum(a.masses) for a in group.atoms.groupby(name).values()]
        vals = group.accumulate("masses", compound=compound)
        assert_almost_equal(vals, ref, decimal=5)

    def test_accumulate_wrongname(self, group):
        with pytest.raises(AttributeError):
            group.accumulate("foo")

    def test_accumulate_wrongcomponent(self, group):
        with pytest.raises(ValueError):
            group.accumulate("masses", compound="foo")

    @pytest.mark.parametrize('level', levels)
    def test_accumulate_nobonds(self, level):
        group = getattr(mda.Universe(GRO), level)
        with pytest.raises(NoDataError):
            group.accumulate("masses", compound="fragments")

    @pytest.mark.parametrize('level', levels)
    def test_accumulate_nomolnums(self, level):
        group = getattr(mda.Universe(GRO), level)
        with pytest.raises(NoDataError):
            group.accumulate("masses", compound="molecules")

    def test_accumulate_array_attribute(self, group):
        a = np.ones((len(group.atoms), 2, 5))
        assert_equal(group.accumulate(a), np.sum(a, axis=0))

    def test_accumulate_array_attribute_wrongshape(self, group):
        with pytest.raises(ValueError):
            group.accumulate(np.ones(len(group.atoms) - 1))

    @pytest.mark.parametrize('name, compound', (('resindices', 'residues'),
                                                ('segindices', 'segments'),
                                                ('molnums', 'molecules'),
                                                ('fragindices', 'fragments')))
    @pytest.mark.parametrize('level', levels)
    def test_accumulate_array_attribute_compounds(self, name, compound, level):
        u = UnWrapUniverse()
        group = getattr(u, level)
        ref = [np.ones((len(a), 2, 5)).sum(axis=0) for a in group.atoms.groupby(name).values()]
        assert_equal(group.accumulate(np.ones((len(group.atoms), 2, 5)), compound=compound), ref)

class TestTotals(object):
    """Tests the functionality of *Group.total*() like total_mass
    and total_charge.
    """
    @pytest.fixture(params=levels)
    def group(self, request):
        u = mda.Universe(PSF, DCD)
        return getattr(u, request.param)

    def test_total_charge(self, group):
        assert_almost_equal(group.total_charge(), -4.0, decimal=4)

    @pytest.mark.parametrize('name, compound',
                             (('resids', 'residues'), ('segids', 'segments'),
                              ('fragindices', 'fragments')))
    def test_total_charge_compounds(self, group, name, compound):
        ref = [sum(a.charges) for a in group.atoms.groupby(name).values()]
        assert_almost_equal(group.total_charge(compound=compound), ref)

    @pytest.mark.filterwarnings(  # Prevents regression of issue #2990
        "error:"
        "Using a non-tuple sequence for multidimensional indexing is deprecated:"
        "FutureWarning"
    )
    def test_total_charge_duplicates(self, group):
        group2 = group + group[0]
        ref = group.total_charge() + group[0].charge
        with pytest.warns(DuplicateWarning) as w:
            assert_almost_equal(group2.total_charge(), ref)
            assert len(w) == 1

    def test_total_mass(self, group):
        assert_almost_equal(group.total_mass(), 23582.043)

    @pytest.mark.parametrize('name, compound',
                             (('resids', 'residues'), ('segids', 'segments'),
                              ('fragindices', 'fragments')))
    def test_total_mass_compounds(self, group, name, compound):
        ref = [sum(a.masses) for a in group.atoms.groupby(name).values()]
        assert_almost_equal(group.total_mass(compound=compound), ref)

    @pytest.mark.filterwarnings(  # Prevents regression of issue #2990
        "error:"
        "Using a non-tuple sequence for multidimensional indexing is deprecated:"
        "FutureWarning"
    )
    def test_total_mass_duplicates(self, group):
        group2 = group + group[0]
        ref = group.total_mass() + group2[0].mass
        with pytest.warns(DuplicateWarning) as w:
            assert_almost_equal(group2.total_mass(), ref)
            assert len(w) == 1
