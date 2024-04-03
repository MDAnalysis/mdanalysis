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
import pytest

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF

from MDAnalysisTests.datafiles import two_water_gro

from numpy.testing import assert_allclose


@pytest.fixture()
def u():
    u = mda.Universe(two_water_gro, in_memory=True)
    u.add_TopologyAttr('chainIDs', u.atoms.resids)
    return u


@pytest.fixture()
def sels(u):
    # modify the coordinates to produce specific test results
    # (NOTE: requires in-memory coordinates to make them permanent)
    for at, (x, y) in zip(u.atoms, zip([1] * 3 + [2] * 3, [2, 1, 3] * 2)):
        at.position = x, y, 0.0
    s1 = u.select_atoms('name OW')
    s2 = u.select_atoms('name HW1 HW2')
    return s1, s2


def test_nbins(u):
    s1 = u.atoms[:3]
    s2 = u.atoms[3:]
    rdf = InterRDF(s1, s2, nbins=412).run()

    assert len(rdf.results.bins) == 412


def test_range(u):
    s1 = u.atoms[:3]
    s2 = u.atoms[3:]
    rmin, rmax = 1.0, 13.0
    rdf = InterRDF(s1, s2, range=(rmin, rmax)).run()

    assert rdf.results.edges[0] == rmin
    assert rdf.results.edges[-1] == rmax


def test_count_sum(sels):
    # OW vs HW
    # should see 8 comparisons in count
    s1, s2 = sels
    rdf = InterRDF(s1, s2).run()
    assert rdf.results.count.sum() == 8


def test_count(sels):
    # should see two distances with 4 counts each
    s1, s2 = sels
    rdf = InterRDF(s1, s2).run()
    assert len(rdf.results.count[rdf.results.count == 4]) == 2


def test_double_run(sels):
    # running rdf twice should give the same result
    s1, s2 = sels
    rdf = InterRDF(s1, s2).run()
    rdf.run()
    assert len(rdf.results.count[rdf.results.count == 4]) == 2


def test_exclusion(sels):
    # should see two distances with 4 counts each
    s1, s2 = sels
    rdf = InterRDF(s1, s2, exclusion_block=(1, 2)).run()
    assert rdf.results.count.sum() == 4


@pytest.mark.parametrize("attr, count", [
    ("residue", 8),
    ("segment", 0),
    ("chain", 8)])
def test_ignore_same_residues(sels, attr, count):
    # should see two distances with 4 counts each
    s1, s2 = sels
    rdf = InterRDF(s2, s2, exclude_same=attr).run()
    assert rdf.rdf[0] == 0
    assert rdf.results.count.sum() == count


def test_ignore_same_residues_fails(sels):
    s1, s2 = sels
    with pytest.raises(ValueError, match="The exclude_same argument to InterRDF must be"):
        InterRDF(s2, s2, exclude_same="unsupported").run()

    with pytest.raises(ValueError, match="The exclude_same argument to InterRDF cannot be used with"):
        InterRDF(s2, s2, exclude_same="residue", exclusion_block=tuple()).run()
        
        
@pytest.mark.parametrize("attr", ("rdf", "bins", "edges", "count"))
def test_rdf_attr_warning(sels, attr):
    s1, s2 = sels
    rdf = InterRDF(s1, s2).run()
    wmsg = f"The `{attr}` attribute was deprecated in MDAnalysis 2.0.0"
    with pytest.warns(DeprecationWarning, match=wmsg):
        getattr(rdf, attr) is rdf.results[attr]


@pytest.mark.parametrize("norm, value", [
    ("density", 1.956823),
    ("rdf", 244602.88385),
    ("none", 4)])
def test_norm(sels, norm, value):
    s1, s2 = sels
    rdf = InterRDF(s1, s2, norm=norm).run()
    assert_allclose(max(rdf.results.rdf), value)


@pytest.mark.parametrize("norm, norm_required", [
    ("Density", "density"), (None, "none")])
def test_norm_values(sels, norm, norm_required):
    s1, s2 = sels
    rdf = InterRDF(s1, s2, norm=norm).run()
    assert rdf.norm == norm_required


def test_unknown_norm(sels):
    s1, s2 = sels
    with pytest.raises(ValueError, match="invalid norm"):
        InterRDF(s1, s2, sels, norm="foo")
