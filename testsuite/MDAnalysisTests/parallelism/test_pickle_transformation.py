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
from numpy.testing import assert_equal

import MDAnalysis as mda

from MDAnalysis.transformations.fit import fit_translation, fit_rot_trans
from MDAnalysis.transformations.positionaveraging import PositionAverager
from MDAnalysis.transformations.rotate import rotateby
from MDAnalysis.transformations.translate import translate, center_in_box
from MDAnalysis.transformations.wrap import wrap, unwrap

from MDAnalysisTests.datafiles import TPR, XTC


@pytest.fixture()
def test_u():
    u = mda.Universe(TPR, XTC)
    return u


u = mda.Universe(TPR, XTC)
ag = u.atoms[0:10]


@pytest.fixture(params=[
    (fit_translation, [ag, ag]),
    (fit_rot_trans, [ag, ag]),
    (PositionAverager, [3]),
    (rotateby, [90, [0, 0, 1], [1, 2, 3]]),
    (translate, [[1, 2, 3]]),
    (center_in_box, [ag]),
    (wrap, [u.atoms]),
    (unwrap, [u.atoms])
])
def transformation(request):
    transform = request.param[0](*request.param[1])
    return transform


def test_transformation_pickle(transformation, test_u):
    ref_result = transformation(test_u.trajectory[5]).positions
    transformation_p = pickle.loads(pickle.dumps(transformation))
    result = transformation_p(test_u.trajectory[5]).positions
    assert_equal(ref_result, result)


@pytest.mark.skip(reason="`Universe` cannot be picked now")
def test_add_transformation_pickle(transformation, test_u):
    test_u.trajectory.add_transformations(transformation)
    test_u_p = pickle.loads(pickle.dumps(test_u))
    test_u.trajectory[0]
    for u_ts, u_p_ts in zip(test_u.trajectory[:5], test_u_p.trajectory[:5]):
        assert_equal(u_ts.positions, u_p_ts.positions)
