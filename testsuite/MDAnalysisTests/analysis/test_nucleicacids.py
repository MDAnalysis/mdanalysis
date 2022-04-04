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
from MDAnalysis.analysis import nuclinfo, nucleicacids
from MDAnalysisTests.datafiles import RNA_PSF, RNA_PDB
from numpy.testing import (
    assert_almost_equal,
    assert_allclose,
)


@pytest.fixture(scope='module')
def u():
    return mda.Universe(RNA_PSF, RNA_PDB)

@pytest.mark.parametrize('i, bp, seg1, seg2, expected_value', (
    ( 1,  2, 'RNAA', 'RNAA', 4.3874702),
    (22, 23, 'RNAA', 'RNAA', 4.1716404),
))
def test_wcbase(u, i, bp, seg1, seg2, expected_value):
    sel = [(f'segid {seg1} and resid {i} and N1', f'segid {seg1} and resid {bp} and N3')]
    WC = nucleicacids.WCBASE(u, sel)
    WC.run()
    assert_almost_equal(WC.results[''], expected_value, decimal=3)