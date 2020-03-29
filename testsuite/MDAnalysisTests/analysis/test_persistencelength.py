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
from __future__ import print_function, division, absolute_import

import pytest

import MDAnalysis as mda
from MDAnalysis.analysis import polymer
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.core.topologyattrs import Bonds

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from numpy.testing import assert_almost_equal, assert_equal

from MDAnalysisTests.datafiles import Plength, TRZ_psf, TRZ


class TestPersistenceLength(object):
    @staticmethod
    @pytest.fixture()
    def u():
        return mda.Universe(Plength)

    @staticmethod
    @pytest.fixture()
    def p(u):
        ags = [r.atoms.select_atoms('name C* N*')
               for r in u.residues]

        p = polymer.PersistenceLength(ags)
        return p

    @staticmethod
    @pytest.fixture()
    def p_run(p):
        return p.run()

    def test_ag_ValueError(self, u):
        ags = [u.atoms[:10], u.atoms[10:110]]
        with pytest.raises(ValueError):
            polymer.PersistenceLength(ags)

    def test_run(self, p_run):
        assert len(p_run.results) == 280

    def test_lb(self, p_run):
        assert_almost_equal(p_run.lb, 1.485, 3)

    def test_fit(self, p_run):
        assert_almost_equal(p_run.lp, 6.504, 3)
        assert len(p_run.fit) == len(p_run.results)

    def test_raise_NoDataError(self, p):
        #Ensure that a NoDataError is raised if perform_fit()
        # is called before the run() method of AnalysisBase
        with pytest.raises(NoDataError):
            p._perform_fit()

    def test_plot_ax_return(self, p_run):
        '''Ensure that a matplotlib axis object is
        returned when plot() is called.'''
        actual = p_run.plot()
        expected = matplotlib.axes.Axes
        assert isinstance(actual, expected)

    def test_plot_with_ax(self, p_run):
        fig, ax = plt.subplots()

        ax2 = p_run.plot(ax=ax)

        assert ax2 is ax
    
    def test_current_axes(self, p_run):
        fig, ax = plt.subplots()
        ax2 = p_run.plot(ax=None)
        assert ax2 is not ax


class TestFitExponential(object):
    x = np.linspace(0, 250, 251)
    a_ref = 20.0
    y = np.exp(-x / a_ref)

    def test_fit_simple(self):
        a = polymer.fit_exponential_decay(self.x, self.y)
        assert a == self.a_ref

    def test_fit_noisy(self):
        noise = np.sin(self.x) * 0.01
        y2 = noise + self.y

        a = polymer.fit_exponential_decay(self.x, y2)

        assert_almost_equal(a, self.a_ref, decimal=3)
        # assert np.rint(a) == self.a_ref


class TestSortBackbone(object):
    @staticmethod
    @pytest.fixture(scope='class')
    def u():
        return mda.Universe(TRZ_psf, TRZ)

    def test_missing_bonds(self):
        u = mda.Universe(Plength)

        with pytest.raises(NoDataError):
            polymer.sort_backbone(u.atoms[:10])

    def test_sortbb(self, u):
        # grab backbone atoms out of order
        # 0 1 4 6 8 - correct
        ag = u.atoms[[4, 1, 0, 8, 6]]

        s_ag = polymer.sort_backbone(ag)

        assert_equal(s_ag.ids, [0, 1, 4, 6, 8])

    def test_not_fragment(self, u):
        # two fragments don't work
        bad_ag = u.residues[0].atoms[:2] + u.residues[1].atoms[:2]
        with pytest.raises(ValueError):
            polymer.sort_backbone(bad_ag)

    def test_branches(self, u):
        # includes side branches, can't sort
        bad_ag = u.atoms[:10]  # include -H etc

        with pytest.raises(ValueError):
            polymer.sort_backbone(bad_ag)

    def test_circular(self):
        u = mda.Universe.empty(6, trajectory=True)
        # circular structure
        bondlist = [(0, 1), (1, 2), (2, 3),
                    (3, 4), (4, 5), (5, 0)]
        u.add_TopologyAttr(Bonds(bondlist))

        with pytest.raises(ValueError) as ex:
            polymer.sort_backbone(u.atoms)
        assert 'cyclical' in str(ex.value)
