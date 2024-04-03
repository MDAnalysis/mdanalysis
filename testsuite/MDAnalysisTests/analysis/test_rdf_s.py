# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2018 The MDAnalysis Development Team and contributors
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

from numpy.testing import assert_allclose, assert_almost_equal

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF_s, InterRDF

from MDAnalysisTests.datafiles import GRO_MEMPROT, XTC_MEMPROT


@pytest.fixture(scope='module')
def u():
    return mda.Universe(GRO_MEMPROT, XTC_MEMPROT)


@pytest.fixture(scope='module')
def sels(u):
    s1 = u.select_atoms('name ZND and resid 289')
    s2 = u.select_atoms(
         '(name OD1 or name OD2) and resid 51 and sphzone 5.0 (resid 289)')
    s3 = u.select_atoms('name ZND and (resid 291 or resid 292)')
    s4 = u.select_atoms('(name OD1 or name OD2) and sphzone 5.0 (resid 291)')
    ags = [[s1, s2], [s3, s4]]
    return ags


@pytest.fixture(scope='module')
def rdf(u, sels):
    return InterRDF_s(u, sels).run()


def test_nbins(u, sels):
    rdf = InterRDF_s(u, sels, nbins=412).run()

    assert len(rdf.results.bins) == 412


def test_range(u, sels):
    rmin, rmax = 1.0, 13.0
    rdf = InterRDF_s(u, sels, range=(rmin, rmax)).run()

    assert rdf.results.edges[0] == rmin
    assert rdf.results.edges[-1] == rmax


def test_count_size(rdf):
    # ZND vs OD1 & OD2
    # should see 2 elements in rdf.results.count
    # 1 element in rdf.results.count[0]
    # 2 elements in rdf.results.count[0][0]
    # 2 elements in rdf.results.count[1]
    # 2 elements in rdf.results.count[1][0]
    # 2 elements in rdf.results.count[1][1]
    assert len(rdf.results.count) == 2
    assert len(rdf.results.count[0]) == 1
    assert len(rdf.results.count[0][0]) == 2
    assert len(rdf.results.count[1]) == 2
    assert len(rdf.results.count[1][0]) == 2
    assert len(rdf.results.count[1][1]) == 2


def test_count(rdf):
    # should see one distance with 5 counts in count[0][0][1]
    # should see one distance with 3 counts in count[1][1][0]
    sel0 = rdf.results.count[0][0][1] == 5
    sel1 = rdf.results.count[1][1][0] == 3
    assert len(rdf.results.count[0][0][1][sel0]) == 1
    assert len(rdf.results.count[1][1][0][sel1]) == 1


def test_double_run(rdf):
    # running rdf twice should give the same result
    sel0 = rdf.results.count[0][0][1] == 5
    sel1 = rdf.results.count[1][1][0] == 3
    assert len(rdf.results.count[0][0][1][sel0]) == 1
    assert len(rdf.results.count[1][1][0][sel1]) == 1


def test_cdf(rdf):
    rdf.get_cdf()
    ref = rdf.results.count[0][0][0].sum()/rdf.n_frames
    assert rdf.results.cdf[0][0][0][-1] == ref


@pytest.mark.parametrize("density, value", [
    (None, 26551.55088100731),    # default, like False (no kwarg, see below)
    (False, 26551.55088100731),
    (True, 0.021915460340071267)])
def test_density(u, sels, density, value):
    kwargs = {'density': density} if density is not None else {}
    rdf = InterRDF_s(u, sels, **kwargs).run()
    assert_almost_equal(max(rdf.results.rdf[0][0][0]), value)
    if not density:
        s1 = u.select_atoms('name ZND and resid 289')
        s2 = u.select_atoms(
                'name OD1 and resid 51 and sphzone 5.0 (resid 289)')
        rdf_ref = InterRDF(s1, s2).run()
        assert_almost_equal(rdf_ref.results.rdf, rdf.results.rdf[0][0][0])


def test_overwrite_norm(u, sels):
    rdf = InterRDF_s(u, sels, norm="rdf", density=True)
    assert rdf.norm == "density"


@pytest.mark.parametrize("norm, value", [
    ("density", 0.021915460340071267),
    ("rdf", 26551.55088100731),
    ("none", 0.6)])
def test_norm(u, sels, norm, value):
    rdf = InterRDF_s(u, sels, norm=norm).run()
    assert_allclose(max(rdf.results.rdf[0][0][0]), value)
    if norm == "rdf":
        s1 = u.select_atoms('name ZND and resid 289')
        s2 = u.select_atoms(
                'name OD1 and resid 51 and sphzone 5.0 (resid 289)')
        rdf_ref = InterRDF(s1, s2).run()
        assert_almost_equal(rdf_ref.results.rdf, rdf.results.rdf[0][0][0])


@pytest.mark.parametrize("norm, norm_required", [
    ("Density", "density"), (None, "none")])
def test_norm_values(u, sels, norm, norm_required):
    rdf = InterRDF_s(u, sels, norm=norm).run()
    assert rdf.norm == norm_required


def test_unknown_norm(u, sels):
    with pytest.raises(ValueError, match="invalid norm"):
        InterRDF_s(u, sels, norm="foo")


@pytest.mark.parametrize("attr", ("rdf", "bins", "edges", "count", "cdf"))
def test_rdf_attr_warning(rdf, attr):
    if attr == "cdf":
        rdf.get_cdf()
    wmsg = f"The `{attr}` attribute was deprecated in MDAnalysis 2.0.0"
    with pytest.warns(DeprecationWarning, match=wmsg):
        getattr(rdf, attr) is rdf.results[attr]
