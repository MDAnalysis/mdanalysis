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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
from six.moves import zip

import pytest

import MDAnalysis as mda
from MDAnalysis.analysis.rdf_s import InterRDF_s

from MDAnalysisTests.datafiles import YiiP_lipids.gro, YiiP_lipids.xtc


@pytest.fixture(scope='module')
def u():
    return mda.Universe('YiiP_lipids.gro', 'YiiP_lipids.xtc')


@pytest.fixture(scope='module')
def sels(u):
    s1 = u.select_atoms('name ZND and resid 289')
    s2 = u.select_atoms('(name OD1 or name OD2) and resid 51 and sphzone 5.0 (resid 289)')
    s3 = u.select_atoms('name ZND and (resid 291 or resid 292)')
    s4 = u.select_atoms('(name OD1 or name OD2) and sphzone 5.0 (resid 291)')
    ag = [[s1, s2], [s3, s4]]
    return ag


def test_nbins(u):
    ag = sels(u)
    rdf = InterRDF_s(u, ag, nbins=412).run()

    assert len(rdf.bins) == 412


def test_range(u):
    ag = sels(u)
    rmin, rmax = 1.0, 13.0
    rdf = InterRDF_s(u, ag, range=(rmin, rmax)).run()

    assert rdf.edges[0] == rmin
    assert rdf.edges[-1] == rmax


def test_count_size(u):
    # ZND vs OD1 & OD2
    # should see 2 elements in rdf.count
    # 
    ag = sels(u)
    rdf = InterRDF_s(u, ag).run()
    assert len(rdf.count) == 2
    assert len(rdf.count[0]) == 1
    assert len(rdf.count[0][0]) == 2
    assert len(rdf.count[1]) == 2
    assert len(rdf.count[1][0]) == 2
    assert len(rdf.count[1][1]) == 2


def test_count(u):
    # should see one distance with 5 counts in count[0][0][1]
    # should see one distance with 3 counts in count[1][1][0]
    ag = sels(u)
    rdf = InterRDF_s(u, ag).run()
    assert len(rdf.count[0][0][1][rdf.count[0][0][1] == 5]) == 1
    assert len(rdf.count[1][1][0][rdf.count[1][1][0] == 3]) == 1


def test_double_run(u):
    # running rdf twice should give the same result
    ag = sels(u)
    rdf = InterRDF_s(u, ag).run()
    assert len(rdf.count[0][0][1][rdf.count[0][0][1] == 5]) == 1
    assert len(rdf.count[1][1][0][rdf.count[1][1][0] == 3]) == 1


def test_cdf(u)
    ag = sels(u)
    rdf = InterRDF_s(u, ag).run()    
    rdf.cdf_s()
    assert rdf.cdf_s[0][0][0][-1] == rdf.count[0][0][0].sum()/rdf.n_frames
    

def test_density(u)
    ag = sels(u)
    rdf_density_on = InterRDF_s(u, ag, density=True).run()
    rdf_density_off = InterRDF_s(u, ag, density=False).run()
    assert max(rdf_density_on.rdf_s[0][0][0]) == 13275.775528444701
    assert max(rdf_density_off.rdf_s[0][0][0]) == 0.021915460340071267
