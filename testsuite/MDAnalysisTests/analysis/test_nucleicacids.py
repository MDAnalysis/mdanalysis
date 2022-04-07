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
from MDAnalysis.analysis.nucleicacids import WCDist, MinorDist, MajorDist, BaseSelect
from MDAnalysisTests.datafiles import RNA_PSF, RNA_PDB
from numpy.testing import (
    assert_almost_equal,
    assert_allclose,
)


@pytest.fixture(scope='module')
def u():
    return mda.Universe(RNA_PSF, RNA_PDB)


def test_wc_dist(u):
    sel = [(BaseSelect('RNAA', 1), BaseSelect('RNAA', 2)),
           (BaseSelect('RNAA', 22), BaseSelect('RNAA', 23))]
    WC = WCDist(u, sel)
    WC.run()

    assert_almost_equal(WC.results['distance'][0], 4.3874702, decimal=3)
    assert_almost_equal(WC.results['distance'][1], 4.1716404, decimal=3)


def test_minor_pair(u):
    sel = [(BaseSelect('RNAA', 3), BaseSelect('RNAA', 17)),
           (BaseSelect('RNAA', 20), BaseSelect('RNAA', 5))]
    MP = MinorDist(u, sel)
    MP.run()

    assert_almost_equal(MP.results['distance'][0], 15.06506, decimal=3)
    assert_almost_equal(MP.results['distance'][1], 3.219116, decimal=3)


def test_major_pair(u):
    sel = [(BaseSelect('RNAA', 2), BaseSelect('RNAA', 12)),
           (BaseSelect('RNAA', 5), BaseSelect('RNAA', 9))]
    MP = MajorDist(u, sel)
    MP.run()

    assert_almost_equal(MP.results['distance'][0], 26.884272, decimal=3)
    assert_almost_equal(MP.results['distance'][1], 13.578535, decimal=3)
