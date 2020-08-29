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
from six.moves import zip

import pytest

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF

from MDAnalysisTests.datafiles import two_water_gro


@pytest.fixture(scope='module')
def u():
    return mda.Universe(two_water_gro)


@pytest.fixture(scope='module')
def sels(u):
    for at, (x, y) in zip(u.atoms, zip([1] * 3 + [2] * 3, [2, 1, 3] * 2)):
        at.position = x, y, 0.0
    s1 = u.select_atoms('name OW')
    s2 = u.select_atoms('name HW1 HW2')
    return s1, s2


def test_nbins(u):
    s1 = u.atoms[:3]
    s2 = u.atoms[3:]
    rdf = InterRDF(s1, s2, nbins=412).run()

    assert len(rdf.bins) == 412


def test_range(u):
    s1 = u.atoms[:3]
    s2 = u.atoms[3:]
    rmin, rmax = 1.0, 13.0
    rdf = InterRDF(s1, s2, range=(rmin, rmax)).run()

    assert rdf.edges[0] == rmin
    assert rdf.edges[-1] == rmax


def test_count_sum(sels):
    # OW vs HW
    # should see 8 comparisons in count
    s1, s2 = sels
    rdf = InterRDF(s1, s2).run()
    assert rdf.count.sum() == 8


def test_count(sels):
    # should see two distances with 4 counts each
    s1, s2 = sels
    rdf = InterRDF(s1, s2).run()
    assert len(rdf.count[rdf.count == 4]) == 2


def test_double_run(sels):
    # running rdf twice should give the same result
    s1, s2 = sels
    rdf = InterRDF(s1, s2).run()
    rdf.run()
    assert len(rdf.count[rdf.count == 4]) == 2


def test_exclusion(sels):
    # should see two distances with 4 counts each
    s1, s2 = sels
    rdf = InterRDF(s1, s2, exclusion_block=(1, 2)).run()
    assert rdf.count.sum() == 4
