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
from __future__ import absolute_import

import pytest
import matplotlib
import numpy as np

from numpy.testing import assert_equal, assert_almost_equal

import MDAnalysis as mda
from MDAnalysis.analysis import secondary_structure as ss
from MDAnalysisTests.datafiles import TPR, XTC, ADK_DSSP, ADK_DSSP_SIMPLE

XTC_modes = ['C', 'E', 'E', 'E', 'E', 'E', 'C', 'C', 'T', 'T', 'S', 'C', 'H',
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C',
             'C', 'E', 'E', 'E', 'T', 'T', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
             'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
             'H', 'H', 'H', 'T', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H',
             'H', 'H', 'H', 'H', 'H', 'H', 'T', 'T', 'S', 'G', 'G', 'G', 'S',
             'S', 'C', 'E', 'E', 'E', 'E', 'S', 'C', 'C', 'C', 'S', 'H', 'H',
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'T', 'T', 'C', 'C', 'C', 'S',
             'E', 'E', 'E', 'E', 'E', 'E', 'C', 'C', 'H', 'H', 'H', 'H', 'H',
             'H', 'H', 'H', 'H', 'T', 'E', 'E', 'E', 'S', 'T', 'T', 'T', 'S',
             'S', 'E', 'E', 'E', 'T', 'T', 'T', 'B', 'C', 'S', 'S', 'S', 'T',
             'T', 'B', 'C', 'T', 'T', 'T', 'C', 'C', 'B', 'C', 'S', 'S', 'S',
             'G', 'G', 'G', 'S', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
             'H', 'H', 'H', 'H', 'H', 'H', 'T', 'S', 'C', 'C', 'E', 'E', 'E',
             'E', 'E', 'C', 'S', 'S', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
             'H', 'H', 'H', 'H', 'C', 'C']

XTC_modes_simple = ['Coil', 'Strand', 'Strand', 'Strand', 'Strand', 'Strand', 'Coil',
                    'Coil', 'Coil', 'Coil', 'Coil', 'Coil', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Coil', 'Coil', 'Coil', 'Strand', 'Strand',
                    'Strand', 'Coil', 'Coil', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Coil', 'Coil',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Coil', 'Coil',
                    'Coil', 'Coil', 'Coil', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Coil', 'Coil', 'Coil', 'Helix', 'Helix', 'Helix', 'Coil', 'Coil',
                    'Coil', 'Strand', 'Strand', 'Strand', 'Strand', 'Coil', 'Coil',
                    'Coil', 'Coil', 'Coil', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Coil', 'Coil',
                    'Coil', 'Coil', 'Coil', 'Coil', 'Strand', 'Strand', 'Strand',
                    'Strand', 'Strand', 'Strand', 'Coil', 'Coil', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Coil', 'Strand', 'Strand', 'Strand', 'Coil', 'Coil', 'Coil',
                    'Coil', 'Coil', 'Coil', 'Strand', 'Strand', 'Strand', 'Coil',
                    'Coil', 'Coil', 'Strand', 'Coil', 'Coil', 'Coil', 'Coil', 'Coil',
                    'Coil', 'Strand', 'Coil', 'Coil', 'Coil', 'Coil', 'Coil', 'Coil',
                    'Strand', 'Coil', 'Coil', 'Coil', 'Coil', 'Helix', 'Helix',
                    'Helix', 'Coil', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Coil', 'Coil', 'Coil', 'Coil', 'Strand',
                    'Strand', 'Strand', 'Strand', 'Strand', 'Coil', 'Coil', 'Coil',
                    'Coil', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Helix',
                    'Helix', 'Helix', 'Helix', 'Helix', 'Helix', 'Coil', 'Coil']

XTC_code_counts = {'G': np.array([7, 6, 3, 3, 3, 8, 3, 6, 0, 3]),
                   'H': np.array([100,  96,  94,  97,  89, 101,  97,  96,  97,  96]),
                   'I': np.array([0, 0, 0, 0, 5, 0, 0, 0, 0, 0]),
                   'E': np.array([24, 28, 28, 27, 28, 29, 28, 31, 30, 30]),
                   'B': np.array([5, 4, 3, 5, 4, 1, 3, 3, 1, 2]),
                   'T': np.array([23, 29, 34, 28, 29, 19, 32, 27, 25, 27]),
                   'S': np.array([20, 20, 22, 22, 23, 27, 22, 21, 31, 21]),
                   'C': np.array([35, 31, 30, 32, 33, 29, 29, 30, 30, 35])}

XTC_simple_counts = {'Helix': np.array([107, 102,  97, 100,  97, 109, 100, 102,  97,  99]),
                     'Coil': np.array([78, 80, 86, 82, 85, 75, 83, 78, 86, 83]),
                     'Strand': np.array([29, 32, 31, 32, 32, 30, 31, 34, 31, 32])}


class TestDSSP(object):

    prot_err_msg = 'Should not have empty secondary structure in protein residues'
    other_err_msg = 'Should not have secondary structure code in non-protein residues'
    plot_err_msg = "DSSP.plot_content() did not produce an Axes instance"

    @pytest.fixture()
    def universe(self):
        return mda.Universe(TPR, XTC)

    @pytest.fixture()
    def dssp(self, universe):
        pytest.importorskip('mdtraj')
        dssp = ss.DSSP(universe, select='backbone', o_name='O O1')
        dssp.run()
        return dssp

    @pytest.fixture()
    def mdtraj_codes(self):
        return np.load(ADK_DSSP)

    @pytest.fixture()
    def mdtraj_simple(self, dssp):
        mdtraj_simple = np.load(ADK_DSSP_SIMPLE).astype(object)
        for code in ['H', 'E', 'C']:
            mdtraj_simple[mdtraj_simple ==
                          code] = dssp.ss_codes_to_simple[code]
        return mdtraj_simple

    def test_correct_assignments(self, dssp, mdtraj_codes, mdtraj_simple):
        assert_equal(dssp.ss_codes, mdtraj_codes)
        assert_equal(dssp.ss_simple, mdtraj_simple)
        assert_equal(dssp.ss_mode, XTC_modes)
        assert_equal(dssp.simple_mode, XTC_modes_simple)
        for k in dssp.ss_counts.keys():
            assert_almost_equal(dssp.ss_counts[k], XTC_code_counts[k])
        for k in dssp.simple_counts.keys():
            assert_almost_equal(dssp.simple_counts[k], XTC_simple_counts[k])

    def test_add_topology_attr(self, universe):
        assert not hasattr(universe._topology, 'secondary_structures')
        dssp = ss.DSSP(universe, select='backbone', o_name='O O1',
                       add_topology_attr=True)
        dssp.run()
        assert hasattr(universe._topology, 'secondary_structures')

        secstruct = universe.residues.secondary_structures
        assert not np.any(secstruct[:214] == ''), self.prot_err_msg
        assert np.all(secstruct[214:] == ''), self.other_err_msg

    def test_non_protein(self, universe, mdtraj_codes, mdtraj_simple):
        n_prot = 214
        dssp = ss.DSSP(universe, select='all', o_name='O O1')
        dssp.run()

        assert_equal(dssp.ss_codes[:, :n_prot], mdtraj_codes)
        assert_equal(dssp.ss_simple[:, :n_prot], mdtraj_simple)
        assert_equal(dssp.ss_mode[:n_prot], XTC_modes)
        assert_equal(dssp.simple_mode[:n_prot], XTC_modes_simple)

        for k in dssp.ss_counts.keys():
            assert_almost_equal(dssp.ss_counts[k], XTC_code_counts[k])
        for k in dssp.simple_counts.keys():
            assert_almost_equal(dssp.simple_counts[k], XTC_simple_counts[k])

        assert np.all(dssp.ss_codes[:, n_prot:] == ''), self.other_err_msg
        assert np.all(dssp.ss_simple[:, n_prot:] == ''), self.other_err_msg
        assert np.all(dssp.ss_mode[n_prot:] == ''), self.other_err_msg
        assert np.all(dssp.simple_mode[n_prot:] == ''), self.other_err_msg

    def test_plot_all_lines(self, dssp):
        ax = dssp.plot_content(kind='line', simple=False)
        assert isinstance(ax, matplotlib.axes.Axes), self.plot_err_msg
        lines = ax.get_lines()[:]
        assert len(lines) == 8

    def test_plot_simple_lines(self, dssp):
        ax = dssp.plot_content(kind='line', simple=True)
        assert isinstance(ax, matplotlib.axes.Axes), self.plot_err_msg
        lines = ax.get_lines()[:]
        assert len(lines) == 3

    def test_plot_all_bars(self, dssp):
        ax = dssp.plot_content(kind='bar', simple=False)
        assert isinstance(ax, matplotlib.axes.Axes), self.plot_err_msg
        rects = [r for r in ax.get_children() if isinstance(
            r, matplotlib.patches.Rectangle)]
        assert len(rects) == (8 * dssp.n_frames) + 1

    def test_plot_n_xticks(self, dssp):
        ax = dssp.plot_content(n_xticks=3)
        locs = ax.get_xticks()
        assert len(locs) == 3

        ax2 = dssp.plot_content(n_xticks=100)
        locs2 = ax2.get_xticks()
        assert len(locs2) == 10
