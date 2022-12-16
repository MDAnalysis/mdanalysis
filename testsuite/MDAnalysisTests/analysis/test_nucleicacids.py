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
from MDAnalysis.analysis.nucleicacids import WatsonCrickDist
from MDAnalysisTests.datafiles import RNA_PSF, RNA_PDB
from numpy.testing import assert_allclose


@pytest.fixture(scope='module')
def u():
    return mda.Universe(RNA_PSF, RNA_PDB)


@pytest.fixture(scope='module')
def wc_rna(u):
    strand: mda.AtomGroup = u.select_atoms("segid RNAA")
    strand1 = [strand.residues[0], strand.residues[21]]
    strand2 = [strand.residues[1], strand.residues[22]]

    WC = WatsonCrickDist(strand1, strand2)
    WC.run()
    return WC


def test_wc_dist(wc_rna):
    assert_allclose(wc_rna.results.pair_distances[0, 0], 4.3874702, atol=1e-3)
    assert_allclose(wc_rna.results.pair_distances[0, 1], 4.1716404, atol=1e-3)


def test_wc_dist_sel_indices_deprecated(wc_rna):
    wmsg = ("Accessing results via selection indices is "
            "deprecated and will be removed in MDAnalysis 2.5.0")

    with pytest.deprecated_call(match=wmsg):
        for i in range(wc_rna._n_sel):
            assert_allclose(
                    wc_rna.results.pair_distances[:, i],
                    wc_rna.results[i][0]
            )


def test_wc_dist_times_deprecated(wc_rna):
    wmsg = ("The `times` results attribute is deprecated and will "
            "be removed in MDAnalysis 2.5.0.")

    with pytest.deprecated_call(match=wmsg):
        assert_allclose(wc_rna.times, wc_rna.results.times)


@pytest.mark.parametrize('key', [2, 'parsnips', 'time', -1])
def test_wc_dis_results_keyerrs(wc_rna, key):

    with pytest.raises(KeyError, match=f"{key}"):
        wc_rna.results[key]
