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
import MDAnalysis as mda
import pytest

from numpy.testing import assert_equal, assert_almost_equal
import numpy as np
import networkx as NX

from MDAnalysis.analysis.leaflet import LeafletFinder, optimize_cutoff, LipidEnrichment
from MDAnalysisTests.datafiles import (Martini_membrane_gro, 
                                       GRO_MEMPROT, 
                                       XTC_MEMPROT,
                                       )

LIPID_HEAD_STRING = "name PO4"

def lines2one(lines):
    """Join lines and squash all whitespace"""
    return " ".join(" ".join(lines).split())

@pytest.fixture()
def universe():
    return mda.Universe(Martini_membrane_gro)


@pytest.fixture()
def lipid_heads(universe):
    return universe.select_atoms(LIPID_HEAD_STRING)


def test_leaflet_finder(universe, lipid_heads):
    lfls = LeafletFinder(universe, lipid_heads, pbc=True)
    top_heads, bottom_heads = lfls.groups()
    # Make top be... on top.
    if top_heads.center_of_geometry()[2] < bottom_heads.center_of_geometry()[2]:
        top_heads,bottom_heads = (bottom_heads,top_heads)
    assert_equal(top_heads.indices, np.arange(1,2150,12),
                 err_msg="Found wrong leaflet lipids")
    assert_equal(bottom_heads.indices, np.arange(2521,4670,12),
                 err_msg="Found wrong leaflet lipids")


def test_string_vs_atomgroup_proper(universe, lipid_heads):
    lfls_ag = LeafletFinder(universe, lipid_heads, pbc=True)
    lfls_string = LeafletFinder(universe, LIPID_HEAD_STRING, pbc=True)
    groups_ag = lfls_ag.groups()
    groups_string = lfls_string.groups()
    assert_equal(groups_string[0].indices, groups_ag[0].indices)
    assert_equal(groups_string[1].indices, groups_ag[1].indices)


def test_optimize_cutoff(universe, lipid_heads):
    cutoff, N = optimize_cutoff(universe, lipid_heads, pbc=True)
    assert N == 2
    assert_almost_equal(cutoff, 10.5, decimal=4)


def test_pbc_on_off(universe, lipid_heads):
    lfls_pbc_on = LeafletFinder(universe, lipid_heads, pbc=True)
    lfls_pbc_off = LeafletFinder(universe, lipid_heads, pbc=False)
    assert lfls_pbc_on.graph.size() > lfls_pbc_off.graph.size()


def test_pbc_on_off_difference(universe, lipid_heads):
    lfls_pbc_on = LeafletFinder(universe, lipid_heads, cutoff=7, pbc=True)
    lfls_pbc_off = LeafletFinder(universe, lipid_heads, cutoff=7, pbc=False)
    pbc_on_graph = lfls_pbc_on.graph
    pbc_off_graph = lfls_pbc_off.graph
    diff_graph = NX.difference(pbc_on_graph, pbc_off_graph)
    assert_equal(set(diff_graph.edges), {(69, 153), (73, 79),
                                        (206, 317), (313, 319)})


@pytest.mark.parametrize("sparse", [True, False, None])
def test_sparse_on_off_none(universe, lipid_heads, sparse):
    lfls_ag = LeafletFinder(universe, lipid_heads, cutoff=15.0, pbc=True,
                                sparse=sparse)
    assert_almost_equal(len(lfls_ag.graph.edges), 1903, decimal=4)


def test_cutoff_update(universe, lipid_heads):
    lfls_ag = LeafletFinder(universe, lipid_heads, cutoff=15.0, pbc=True)
    lfls_ag.update(cutoff=1.0)
    assert_almost_equal(lfls_ag.cutoff, 1.0, decimal=4)
    assert_almost_equal(len(lfls_ag.groups()), 360, decimal=4)


def test_cutoff_update_default(universe, lipid_heads):
    lfls_ag = LeafletFinder(universe, lipid_heads, cutoff=15.0, pbc=True)
    lfls_ag.update()
    assert_almost_equal(lfls_ag.cutoff, 15.0, decimal=4)
    assert_almost_equal(len(lfls_ag.groups()), 2, decimal=4)


def test_write_selection(universe, lipid_heads, tmpdir):
    lfls_ag = LeafletFinder(universe, lipid_heads, cutoff=15.0, pbc=True)
    with tmpdir.as_cwd():
        filename = lfls_ag.write_selection('leaflet.vmd')
        expected_output = lines2one([
            """# leaflets based on select=<AtomGroup with 360 atoms> cutoff=15.000000
        # MDAnalysis VMD selection
        atomselect macro leaflet_1 {index 1 13 25 37 49 61 73 85 \\
        97 109 121 133 145 157 169 181 \\
        193 205 217 229 241 253 265 277 \\
        289 301 313 325 337 349 361 373 \\
        385 397 409 421 433 445 457 469 \\
        481 493 505 517 529 541 553 565 \\
        577 589 601 613 625 637 649 661 \\
        673 685 697 709 721 733 745 757 \\
        769 781 793 805 817 829 841 853 \\
        865 877 889 901 913 925 937 949 \\
        961 973 985 997 1009 1021 1033 1045 \\
        1057 1069 1081 1093 1105 1117 1129 1141 \\
        1153 1165 1177 1189 1201 1213 1225 1237 \\
        1249 1261 1273 1285 1297 1309 1321 1333 \\
        1345 1357 1369 1381 1393 1405 1417 1429 \\
        1441 1453 1465 1477 1489 1501 1513 1525 \\
        1537 1549 1561 1573 1585 1597 1609 1621 \\
        1633 1645 1657 1669 1681 1693 1705 1717 \\
        1729 1741 1753 1765 1777 1789 1801 1813 \\
        1825 1837 1849 1861 1873 1885 1897 1909 \\
        1921 1933 1945 1957 1969 1981 1993 2005 \\
        2017 2029 2041 2053 2065 2077 2089 2101 \\
        2113 2125 2137 2149 }
        # MDAnalysis VMD selection
        atomselect macro leaflet_2 {index 2521 2533 2545 2557 2569 2581 2593 2605 \\
        2617 2629 2641 2653 2665 2677 2689 2701 \\
        2713 2725 2737 2749 2761 2773 2785 2797 \\
        2809 2821 2833 2845 2857 2869 2881 2893 \\
        2905 2917 2929 2941 2953 2965 2977 2989 \\
        3001 3013 3025 3037 3049 3061 3073 3085 \\
        3097 3109 3121 3133 3145 3157 3169 3181 \\
        3193 3205 3217 3229 3241 3253 3265 3277 \\
        3289 3301 3313 3325 3337 3349 3361 3373 \\
        3385 3397 3409 3421 3433 3445 3457 3469 \\
        3481 3493 3505 3517 3529 3541 3553 3565 \\
        3577 3589 3601 3613 3625 3637 3649 3661 \\
        3673 3685 3697 3709 3721 3733 3745 3757 \\
        3769 3781 3793 3805 3817 3829 3841 3853 \\
        3865 3877 3889 3901 3913 3925 3937 3949 \\
        3961 3973 3985 3997 4009 4021 4033 4045 \\
        4057 4069 4081 4093 4105 4117 4129 4141 \\
        4153 4165 4177 4189 4201 4213 4225 4237 \\
        4249 4261 4273 4285 4297 4309 4321 4333 \\
        4345 4357 4369 4381 4393 4405 4417 4429 \\
        4441 4453 4465 4477 4489 4501 4513 4525 \\
        4537 4549 4561 4573 4585 4597 4609 4621 \\
        4633 4645 4657 4669 }

"""])

        assert lines2one(open('leaflet.vmd').readlines()) == expected_output


def test_component_index_is_not_none(universe, lipid_heads):
    lfls_ag = LeafletFinder(universe, lipid_heads, cutoff=15.0, pbc=True)
    assert_almost_equal(len(lfls_ag.groups(component_index=0)), 180, decimal=4)

class BaseTestLipidEnrichment:
    pytest.importorskip('sklearn')
    protein_sel = 'protein'

    @pytest.fixture()
    def universe(self):
        return mda.Universe(*self.files)

    @pytest.fixture()
    def lipen(self, universe):
        return LipidEnrichment(universe, select_protein=self.protein_sel,
                               select_headgroup=self.headgroup_sel,
                               select_residues=self.lipid_sel,
                               enrichment_cutoff=self.cutoff).run()


class TestLipidEnrichmentMembrane(BaseTestLipidEnrichment):
    files = [Martini_membrane_gro]
    lipid_sel = 'resname DPPC'
    headgroup_sel = 'name PO4'
    cutoff = 6
    n_lipids = 360

    def test_leaflet_finder(self, lipen):
        top_heads, bottom_heads = lipen.leaflet_headgroups
        n_residues = len(top_heads) + len(bottom_heads)
        assert_equal(n_residues, self.n_lipids)
        assert_equal(top_heads.indices, np.arange(1,2150,12),
                    err_msg="Found wrong leaflet lipids")
        assert_equal(bottom_heads.indices, np.arange(2521,4670,12),
                    err_msg="Found wrong leaflet lipids")
    
    def test_empty_results(self, lipen):
        assert len(lipen.leaflets) == 2
        top, bottom = lipen.leaflets
        assert len(top) == 2
        assert len(bottom) == 2
        assert_equal(top['DPPC']['Near protein'], 0)
        assert_equal(top['all']['Near protein'], 0)

class TestLipidEnrichmentMemProtAA(BaseTestLipidEnrichment):
    files = [GRO_MEMPROT, XTC_MEMPROT]
    lipid_sel = 'resname POPE POPG'
    headgroup_sel = 'name P*'
    cutoff = 4

    n_u_POPE = 113
    n_l_POPE = 108
    n_u_POPG = 28
    n_l_POPG = 27
    n_upper = 113+28
    n_lower = 108+27
    n_lipids = 276

    def test_leaflet_finder(self, lipen):
        # by visual inspection
        tops, bottoms = lipen.leaflet_residues
        n_residues = len(tops) + len(bottoms)
        assert_equal(n_residues, self.n_lipids)
        top_rids = list(range(297, 410)) + list(range(518, 546))
        assert_equal(tops.resids, top_rids,
                    err_msg="Found wrong leaflet lipids")
        bottom_rids = list(range(410, 518)) + list(range(546, 573))
        assert_equal(bottoms.resids, bottom_rids,
                     err_msg="Found wrong leaflet lipids")

    def test_results(self, lipen):
        upper, lower = lipen.leaflets
        assert_equal(set(upper.keys()), {'POPE', 'POPG', 'all'})
        assert_equal(set(lower.keys()), {'POPE', 'POPG', 'all'})

        # counts near protein by selection language ('around 4 protein')
        n_u_pope = np.array([4, 7, 6, 2, 3])
        n_l_pope = np.array([8, 11, 10, 10, 7])
        n_u_popg = np.zeros(5)
        n_l_popg = np.array([5, 5, 3, 6, 7])
        n_u = n_u_pope+n_u_popg
        n_l = n_l_pope+n_l_popg
        assert_equal(upper['POPE']['Near protein'], n_u_pope)
        assert_equal(lower['POPE']['Near protein'], n_l_pope)
        assert_equal(upper['POPG']['Near protein'], n_u_popg)
        assert_equal(lower['POPG']['Near protein'], n_l_popg)
        assert_equal(upper['all']['Near protein'], n_u)
        assert_equal(lower['all']['Near protein'], n_l)

        # Fraction
        f_l_pope = n_l_pope/n_l
        f_l_popg = 1-f_l_pope
        assert_almost_equal(upper['POPE']['Fraction near protein'], np.ones(5))
        assert_almost_equal(lower['POPE']['Fraction near protein'], f_l_pope)
        assert_almost_equal(upper['POPG']['Fraction near protein'], np.zeros(5))
        assert_almost_equal(lower['POPG']['Fraction near protein'], f_l_popg)

        # enrichment
        assert_almost_equal(upper['POPE']['Enrichment'], self.n_upper/self.n_u_POPE)
        assert_almost_equal(lower['POPE']['Enrichment'], f_l_pope/0.8)
        assert_almost_equal(upper['POPG']['Enrichment'], 0)
        assert_almost_equal(lower['POPG']['Enrichment'], f_l_popg/0.2)

    def test_results_summary(self, lipen):
        upper, lower = lipen.leaflets_summary
        assert_equal(set(upper.keys()), {'POPE', 'POPG', 'all'})
        assert_equal(set(lower.keys()), {'POPE', 'POPG', 'all'})

        # counts near protein
        assert_almost_equal(upper['POPE']['Average near protein'], 4.4)
        assert_almost_equal(lower['POPE']['Average near protein'], 9.2)
        assert_almost_equal(upper['POPG']['Average near protein'], 0)
        assert_almost_equal(lower['POPG']['Average near protein'], 5.2)
        assert_almost_equal(upper['all']['Average near protein'], 4.4)
        assert_almost_equal(lower['all']['Average near protein'], 14.4)
        assert_almost_equal(upper['POPE']['SD near protein'], 1.8547, decimal=4)
        assert_almost_equal(lower['POPE']['SD near protein'], 1.4697, decimal=4)
        assert_almost_equal(upper['POPG']['SD near protein'], 0)
        assert_almost_equal(lower['POPG']['SD near protein'], 1.3266, decimal=4)

        # Fraction
        assert_almost_equal(upper['POPE']['Average fraction near protein'], 1)
        assert_almost_equal(lower['POPE']['Average fraction near protein'], 0.6389, decimal=4)
        assert_almost_equal(upper['POPG']['Average fraction near protein'], 0)
        assert_almost_equal(lower['POPG']['Average fraction near protein'], 0.3611, decimal=4)
        assert_almost_equal(upper['POPE']['SD fraction near protein'], 0)
        assert_almost_equal(lower['POPE']['SD fraction near protein'], 0.0888, decimal=4)
        assert_almost_equal(upper['POPG']['SD fraction near protein'], 0)
        assert_almost_equal(lower['POPG']['SD fraction near protein'], 0.0888, decimal=4)

        # enrichment
        assert_almost_equal(upper['POPE']['Average enrichment'], 1.2478, decimal=4)
        assert_almost_equal(lower['POPE']['Average enrichment'], 0.7986, decimal=4)
        assert_almost_equal(upper['POPG']['Average enrichment'], 0)
        assert_almost_equal(lower['POPG']['Average enrichment'], 1.8056, decimal=4)
        assert_almost_equal(upper['POPE']['SD enrichment'], 0)
        assert_almost_equal(lower['POPE']['SD enrichment'], 0.1109, decimal=4)
        assert_almost_equal(upper['POPG']['SD enrichment'], 0)
        assert_almost_equal(lower['POPG']['SD enrichment'], 0.4438, decimal=4)

    def test_results_summary_df(self, lipen):
        df = lipen.summary_as_dataframe()
        upper = df[df.Leaflet == 1]
        lower = df[df.Leaflet == 2]

        # counts near protein
        assert_almost_equal(upper.loc['POPE', 'Average near protein'], 4.4)
        assert_almost_equal(lower.loc['POPE', 'Average near protein'], 9.2)
        assert_almost_equal(upper.loc['POPG', 'Average near protein'], 0)
        assert_almost_equal(lower.loc['POPG', 'Average near protein'], 5.2)
        assert_almost_equal(upper.loc['all', 'Average near protein'], 4.4)
        assert_almost_equal(lower.loc['all', 'Average near protein'], 14.4)
        assert_almost_equal(upper.loc['POPE', 'SD near protein'], 1.8547, decimal=4)
        assert_almost_equal(lower.loc['POPE', 'SD near protein'], 1.4697, decimal=4)
        assert_almost_equal(upper.loc['POPG', 'SD near protein'], 0)
        assert_almost_equal(lower.loc['POPG', 'SD near protein'], 1.3266, decimal=4)

        # Fraction
        assert_almost_equal(upper.loc['POPE', 'Average fraction near protein'], 1)
        assert_almost_equal(lower.loc['POPE', 'Average fraction near protein'], 0.6389, decimal=4)
        assert_almost_equal(upper.loc['POPG', 'Average fraction near protein'], 0)
        assert_almost_equal(lower.loc['POPG', 'Average fraction near protein'], 0.3611, decimal=4)
        assert_almost_equal(upper.loc['POPE', 'SD fraction near protein'], 0)
        assert_almost_equal(lower.loc['POPE', 'SD fraction near protein'], 0.0888, decimal=4)
        assert_almost_equal(upper.loc['POPG', 'SD fraction near protein'], 0)
        assert_almost_equal(lower.loc['POPG', 'SD fraction near protein'], 0.0888, decimal=4)

        # enrichment
        assert_almost_equal(upper.loc['POPE', 'Average enrichment'], 1.2478, decimal=4)
        assert_almost_equal(lower.loc['POPE', 'Average enrichment'], 0.7986, decimal=4)
        assert_almost_equal(upper.loc['POPG', 'Average enrichment'], 0)
        assert_almost_equal(lower.loc['POPG', 'Average enrichment'], 1.8056, decimal=4)
        assert_almost_equal(upper.loc['POPE', 'SD enrichment'], 0)
        assert_almost_equal(lower.loc['POPE', 'SD enrichment'], 0.1109, decimal=4)
        assert_almost_equal(upper.loc['POPG', 'SD enrichment'], 0)
        assert_almost_equal(lower.loc['POPG', 'SD enrichment'], 0.4438, decimal=4)