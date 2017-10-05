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

import pytest
from numpy.testing import assert_almost_equal

from MDAnalysisTests.datafiles import TXYZ, ARC

import MDAnalysis as mda


@pytest.fixture
def TXYZ_U():
    return mda.Universe(TXYZ)


@pytest.fixture
def ARC_U():
    return mda.Universe(ARC)


def test_txyz_positions(TXYZ_U):
    assert_almost_equal(TXYZ_U.atoms.positions[0],
                        [-6.553398, -1.854369, 0.000000])


def test_arc_positions(ARC_U):
    assert_almost_equal(ARC_U.atoms.positions[0],
                        [-6.553398, -1.854369, 0.000000])


def test_arc_positions_frame_2(ARC_U):
    ARC_U.trajectory[1]

    assert_almost_equal(ARC_U.atoms.positions[0],
                        [-0.231579, -0.350841, -0.037475])


def test_arc_traj_legnth(ARC_U):
    assert len(ARC_U.trajectory) == 2
