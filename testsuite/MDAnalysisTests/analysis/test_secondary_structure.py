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
from MDAnalysisTests import executable_not_found
from MDAnalysisTests.datafiles import (TPR, XTC,
                                       ADK_DSSP, ADK_DSSP_SIMPLE,
                                       ADK_MKDSSP, MKDSSP_phi_psi_sasa,
                                       ADK_STRIDE, STRIDE_phi_psi_sasa,
                                       )


def mdtraj_found():
    try:
        import mdtraj
    except ImportError:
        return False
    else:
        return True


class BaseTestSecondaryStructure(object):
    prot_err_msg = 'Protein residues do not match secondary structure found in {}'
    other_err_msg = 'Should not have secondary structure code in non-protein residues'
    plot_err_msg = 'plot_content() did not produce an Axes instance'

    n_prot = 214

    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(TPR, XTC)

    def test_correct_assignments(self, dssp, all_codes):
        assert_equal(dssp.ss_codes, all_codes)
        assert_equal(dssp.ss_mode[::3], self.modes_3)
        assert_equal(dssp.simple_mode[:10], self.modes_simple_10)
        for k in dssp.ss_counts.keys():
            assert_almost_equal(dssp.ss_counts[k], self.code_counts[k])
        for k in dssp.simple_counts.keys():
            assert_almost_equal(dssp.simple_counts[k], self.simple_counts[k])

    def test_add_topology_attr(self, universe):
        assert not hasattr(universe._topology, 'secondary_structures')
        dssp = self.analysis_class(universe, select='backbone',
                                   add_topology_attr=True)
        dssp.run()
        assert hasattr(universe._topology, 'secondary_structures')

        secstruct = universe.residues.secondary_structures
        modes = dssp.ss_mode[:self.n_prot]
        prot_msg = self.prot_err_msg.format(self.analysis_class.__name__)
        assert np.all(secstruct[:self.n_prot] == modes), prot_msg
        assert np.all(secstruct[self.n_prot:] == ''), self.other_err_msg

    def test_non_protein(self, universe, all_codes):
        dssp = self.analysis_class(universe, select='all')
        dssp.run()

        assert_equal(dssp.ss_codes[:, :self.n_prot], all_codes)
        assert_equal(dssp.ss_mode[:self.n_prot:3], self.modes_3)
        assert_equal(dssp.simple_mode[:10], self.modes_simple_10)

        for k in dssp.ss_counts.keys():
            assert_almost_equal(dssp.ss_counts[k], self.code_counts[k])
        for k in dssp.simple_counts.keys():
            assert_almost_equal(dssp.simple_counts[k], self.simple_counts[k])

        assert np.all(dssp.ss_codes[:, self.n_prot:] == ''), self.other_err_msg
        assert np.all(dssp.ss_simple[:, self.n_prot:]
                      == ''), self.other_err_msg
        assert np.all(dssp.ss_mode[self.n_prot:] == ''), self.other_err_msg
        assert np.all(dssp.simple_mode[self.n_prot:] == ''), self.other_err_msg

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

@pytest.mark.skipif(not mdtraj_found(),
                    reason="Tests skipped because mdtraj not installed")
class TestDSSP(BaseTestSecondaryStructure):

    modes_3 = ['C', 'E', 'C', 'T', 'H', 'H', 'H', 'H', 'C', 'E', 'T', 'H', 'H',
               'H', 'C', 'H', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'H', 'T', 'G',
               'S', 'E', 'S', 'C', 'H', 'H', 'H', 'T', 'C', 'E', 'E', 'C', 'H',
               'H', 'H', 'E', 'T', 'S', 'E', 'T', 'C', 'S', 'B', 'T', 'C', 'S',
               'G', 'S', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'S', 'E',
               'E', 'S', 'H', 'H', 'H', 'H', 'C']

    modes_simple_10 = ['Coil', 'Strand', 'Strand', 'Strand', 'Strand', 'Strand', 'Coil',
                       'Coil', 'Coil', 'Coil']

    code_counts = {'G': np.array([7, 6, 3, 3, 3, 8, 3, 6, 0, 3]),
                   'H': np.array([100,  96,  94,  97,  89, 101,  97,  96,  97,  96]),
                   'I': np.array([0, 0, 0, 0, 5, 0, 0, 0, 0, 0]),
                   'E': np.array([24, 28, 28, 27, 28, 29, 28, 31, 30, 30]),
                   'B': np.array([5, 4, 3, 5, 4, 1, 3, 3, 1, 2]),
                   'T': np.array([23, 29, 34, 28, 29, 19, 32, 27, 25, 27]),
                   'S': np.array([20, 20, 22, 22, 23, 27, 22, 21, 31, 21]),
                   'C': np.array([35, 31, 30, 32, 33, 29, 29, 30, 30, 35])}

    simple_counts = {'Helix': np.array([107, 102,  97, 100,  97, 109, 100, 102,  97,  99]),
                     'Coil': np.array([78, 80, 86, 82, 85, 75, 83, 78, 86, 83]),
                     'Strand': np.array([29, 32, 31, 32, 32, 30, 31, 34, 31, 32])}

    analysis_class = ss.DSSP

    @pytest.fixture()
    def dssp(self, universe):
        pytest.importorskip('mdtraj')
        dssp = ss.DSSP(universe, select='backbone', o_name='O O1')
        dssp.run()
        return dssp

    @pytest.fixture()
    def all_codes(self):
        return np.load(ADK_DSSP)

    @pytest.fixture()
    def all_simple(self, dssp):
        mdtraj_simple = np.load(ADK_DSSP_SIMPLE).astype(object)
        for code in ['H', 'E', 'C']:
            matches = mdtraj_simple == code
            mdtraj_simple[matches] = dssp.ss_codes_to_simple[code]
        return mdtraj_simple

    def test_correct_assignments(self, dssp, all_codes, all_simple):
        super(TestDSSP, self).test_correct_assignments(dssp, all_codes)
        assert_equal(dssp.ss_simple, all_simple)

    def test_non_protein(self, universe, all_codes, all_simple):
        dssp = self.analysis_class(universe, select='all')
        dssp.run()

        assert_equal(dssp.ss_codes[:, :self.n_prot], all_codes)
        assert_equal(dssp.ss_simple[:, :self.n_prot], all_simple)
        assert_equal(dssp.ss_mode[:self.n_prot:3], self.modes_3)
        assert_equal(dssp.simple_mode[:10], self.modes_simple_10)

        for k in dssp.ss_counts.keys():
            assert_almost_equal(dssp.ss_counts[k], self.code_counts[k])
        for k in dssp.simple_counts.keys():
            assert_almost_equal(dssp.simple_counts[k], self.simple_counts[k])

        assert np.all(dssp.ss_codes[:, self.n_prot:] == ''), self.other_err_msg
        assert np.all(dssp.ss_simple[:, self.n_prot:]
                      == ''), self.other_err_msg
        assert np.all(dssp.ss_mode[self.n_prot:] == ''), self.other_err_msg
        assert np.all(dssp.simple_mode[self.n_prot:] == ''), self.other_err_msg


class BaseTestSecondaryStructureWrapper(BaseTestSecondaryStructure):

    phi_psi_sasa_file = None

    def test_executable_not_found(self, universe):
        with pytest.raises(OSError) as exc:
            dssp = self.analysis_class(universe, executable='foo')
        assert 'executable not found' in str(exc.value)

    def test_correct_assignments(self, dssp, all_codes):
        super(BaseTestSecondaryStructureWrapper, self).test_correct_assignments(dssp, all_codes)
        phi, psi, sasa = np.load(self.phi_psi_sasa_file)
        assert_almost_equal(dssp.phi, phi)
        assert_almost_equal(dssp.psi, psi)
        assert_almost_equal(dssp.sasa, sasa)

    

@pytest.mark.skipif(executable_not_found('mkdssp'),
                    reason="Tests skipped because mkdssp (DSSP) not found")
class TestDSSPWrapper(BaseTestSecondaryStructureWrapper):

    phi_psi_sasa_file = MKDSSP_phi_psi_sasa

    modes_3 = ['C', 'E', 'C', 'T', 'H', 'H', 'H', 'H', 'C', 'E', 'T', 'H', 'H',
               'H', 'C', 'H', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'H', 'T', 'G',
               'S', 'E', 'S', 'C', 'H', 'H', 'H', 'T', 'C', 'E', 'E', 'C', 'H',
               'H', 'H', 'E', 'T', 'C', 'E', 'T', 'C', 'S', 'B', 'T', 'C', 'C',
               'G', 'S', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'S', 'E',
               'E', 'S', 'H', 'H', 'H', 'H', '']

    modes_simple_10 = ['Coil', 'Strand', 'Strand', 'Strand', 'Strand', 'Strand', 'Coil',
                       'Coil', 'Coil', 'Coil']

    code_counts = {'G': np.array([7, 6, 3, 3, 3, 8, 3, 6, 0, 3]),
                   'H': np.array([99,  96,  94,  97,  89, 101,  97,  96,  97,  96]),
                   'I': np.array([0, 0, 0, 0, 5, 0, 0, 0, 0, 0]),
                   'E': np.array([20, 28, 28, 23, 28, 29, 28, 31, 30, 30]),
                   'B': np.array([5, 4, 3, 5, 4, 1, 3, 3, 1, 2]),
                   'T': np.array([23, 28, 33, 28, 28, 19, 25, 27, 22, 26]),
                   'S': np.array([13, 17, 19, 17, 17, 23, 19, 10, 25, 18]),
                   'C': np.array([46, 34, 33, 40, 39, 32, 38, 40, 38, 38])}

    simple_counts = {'Helix': np.array([106, 102,  97, 100,  97, 109, 100, 102,  97,  99]),
                     'Coil': np.array([82, 79, 85, 85, 84, 74, 82, 77, 85, 82]),
                     'Strand': np.array([25, 32, 31, 28, 32, 30, 31, 34, 31, 32])}

    analysis_class = ss.DSSPWrapper

    @pytest.fixture(scope='class')
    def dssp(self, universe):
        dssp = ss.DSSPWrapper(universe, executable='mkdssp',
                              select='backbone')
        dssp.run()
        return dssp

    @pytest.fixture()
    def all_codes(self):
        return np.load(ADK_MKDSSP)


@pytest.mark.skipif(executable_not_found('stride'),
                    reason="Tests skipped because STRIDE not found")
class TestStrideWrapper(BaseTestSecondaryStructureWrapper):

    phi_psi_sasa_file = STRIDE_phi_psi_sasa

    modes_3 = ['C', 'E', 'C', 'T', 'H', 'H', 'H', 'H', 'C', 'E', 'C', 'H', 'H',
               'H', 'C', 'H', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'H', 'C', 'T',
               'T', 'E', 'T', 'T', 'H', 'H', 'H', 'C', 'C', 'E', 'E', 'C', 'H',
               'H', 'H', 'E', 'T', 'C', 'E', 'T', 'T', 'T', 'B', 'T', 'C', 'C',
               'G', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'E',
               'E', 'T', 'H', 'H', 'H', 'H', 'C']

    modes_simple_10 = ['Coil', 'Strand', 'Strand', 'Strand', 'Strand', 'Strand', 'Coil',
                       'Coil', 'Coil', 'Coil']

    code_counts = {'G': np.array([3, 3, 0, 3, 3, 0, 3, 3, 0, 3]),
                   'H': np.array([101,  99, 100,  93, 102, 106, 102, 105, 101, 102]),
                   'I': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                   'E': np.array([26, 29, 28, 28, 28, 32, 27, 30, 30, 30]),
                   'B': np.array([3, 4, 3, 3, 6, 1, 3, 3, 1, 2]),
                   'T': np.array([38, 41, 43, 48, 39, 38, 34, 34, 41, 38]),
                   'S': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                   'C': np.array([43, 38, 40, 39, 36, 37, 45, 39, 41, 39])}

    simple_counts = {'Helix': np.array([104, 102, 100,  96, 105, 106, 105, 108, 101, 105]),
                     'Coil': np.array([81, 79, 83, 87, 75, 75, 79, 73, 82, 77]),
                     'Strand': np.array([29, 33, 31, 31, 34, 33, 30, 33, 31, 32])}

    analysis_class = ss.STRIDEWrapper

    @pytest.fixture(scope='class')
    def dssp(self, universe):
        dssp = ss.STRIDEWrapper(universe, executable='stride',
                                select='backbone')
        dssp.run()
        return dssp

    @pytest.fixture()
    def all_codes(self):
        return np.load(ADK_STRIDE)
