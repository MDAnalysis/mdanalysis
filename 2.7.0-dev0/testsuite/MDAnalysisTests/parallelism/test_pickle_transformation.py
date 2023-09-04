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
import pickle
from numpy.testing import assert_almost_equal

import MDAnalysis as mda

from MDAnalysis.transformations.fit import fit_translation, fit_rot_trans
from MDAnalysis.transformations.positionaveraging import PositionAverager
from MDAnalysis.transformations.rotate import rotateby
from MDAnalysis.transformations.translate import translate, center_in_box
from MDAnalysis.transformations.wrap import wrap, unwrap

from MDAnalysisTests.datafiles import PSF_TRICLINIC, DCD_TRICLINIC


@pytest.fixture(params=[
    (PSF_TRICLINIC, DCD_TRICLINIC),
])
def u(request):
    top, traj = request.param
    return mda.Universe(top, traj)


@pytest.fixture()
def fit_translation_transformation(u):
    ag = u.atoms[0:10]
    return fit_translation(ag, ag)


@pytest.fixture()
def fit_rot_trans_transformation(u):
    ag = u.atoms[0:10]
    return fit_rot_trans(ag, ag)


@pytest.fixture()
def PositionAverager_transformation(u):
    return PositionAverager(3)


@pytest.fixture()
def rotateby_transformation(u):
    return rotateby(90, [0, 0, 1], [1, 2, 3])


@pytest.fixture()
def translate_transformation(u):
    return translate([1, 2, 3])


@pytest.fixture()
def center_in_box_transformation(u):
    ag = u.atoms[0:10]
    return center_in_box(ag)


@pytest.fixture()
def wrap_transformation(u):
    ag = u.atoms
    return wrap(ag)


@pytest.fixture()
def unwrap_transformation(u):
    ag = u.atoms
    return unwrap(ag)


def test_add_fit_translation_pickle(fit_translation_transformation, u):
    u.trajectory.add_transformations(fit_translation_transformation)
    u_p = pickle.loads(pickle.dumps(u))
    u.trajectory[0]
    for u_ts, u_p_ts in zip(u.trajectory[:5], u_p.trajectory[:5]):
        assert_almost_equal(u_ts.positions, u_p_ts.positions)


def test_add_fit_rot_trans_pickle(fit_rot_trans_transformation,
                                  u):
    u.trajectory.add_transformations(fit_rot_trans_transformation)
    u_p = pickle.loads(pickle.dumps(u))
    u.trajectory[0]
    for u_ts, u_p_ts in zip(u.trajectory[:5], u_p.trajectory[:5]):
        assert_almost_equal(u_ts.positions, u_p_ts.positions)


def test_add_PositionAverager_pickle(PositionAverager_transformation, u):
    u.trajectory.add_transformations(PositionAverager_transformation)
    u_p = pickle.loads(pickle.dumps(u))
    u.trajectory[0]
    for u_ts, u_p_ts in zip(u.trajectory[:5], u_p.trajectory[:5]):
        assert_almost_equal(u_ts.positions, u_p_ts.positions)


def test_add_rotateby_pickle(rotateby_transformation, u):
    u.trajectory.add_transformations(rotateby_transformation)
    u_p = pickle.loads(pickle.dumps(u))
    u.trajectory[0]
    for u_ts, u_p_ts in zip(u.trajectory[:5], u_p.trajectory[:5]):
        assert_almost_equal(u_ts.positions, u_p_ts.positions)


def test_add_translate_pickle(translate_transformation, u):
    u.trajectory.add_transformations(translate_transformation)
    u_p = pickle.loads(pickle.dumps(u))
    u.trajectory[0]
    for u_ts, u_p_ts in zip(u.trajectory[:5], u_p.trajectory[:5]):
        assert_almost_equal(u_ts.positions, u_p_ts.positions)


def test_add_center_in_box_pickle(center_in_box_transformation, u):
    u.trajectory.add_transformations(center_in_box_transformation)
    u_p = pickle.loads(pickle.dumps(u))
    u.trajectory[0]
    for u_ts, u_p_ts in zip(u.trajectory[:5], u_p.trajectory[:5]):
        assert_almost_equal(u_ts.positions, u_p_ts.positions)


def test_add_wrap_pickle(wrap_transformation, u):
    u.trajectory.add_transformations(wrap_transformation)
    u_p = pickle.loads(pickle.dumps(u))
    u.trajectory[0]
    for u_ts, u_p_ts in zip(u.trajectory[:5], u_p.trajectory[:5]):
        assert_almost_equal(u_ts.positions, u_p_ts.positions)


def test_add_unwrap_pickle(unwrap_transformation, u):
    u.trajectory.add_transformations(unwrap_transformation)
    u_p = pickle.loads(pickle.dumps(u))
    u.trajectory[0]
    for u_ts, u_p_ts in zip(u.trajectory[:5], u_p.trajectory[:5]):
        assert_almost_equal(u_ts.positions, u_p_ts.positions)
