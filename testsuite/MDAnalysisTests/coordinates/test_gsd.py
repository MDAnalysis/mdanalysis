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
import os

import pytest
from numpy.testing import assert_almost_equal

import MDAnalysis as mda
from MDAnalysis.coordinates.GSD import GSDReader
from MDAnalysisTests.datafiles import GSD


@pytest.fixture
def GSD_U():
    return mda.Universe(GSD)


def test_gsd_positions(GSD_U):
    # first frame first particle
    ts = GSD_U.trajectory[0]
    assert_almost_equal(GSD_U.atoms.positions[0],
                        [ -5.4000001 , -10.19999981, -10.19999981])
    # second frame first particle
    ts = GSD_U.trajectory[1]
    assert_almost_equal(GSD_U.atoms.positions[0],
                        [ -5.58348083,  -9.98546982, -10.17657185])


def test_gsd_n_frames(GSD_U):
    assert len(GSD_U.trajectory) == 2


def test_gsd_dimensions(GSD_U):
    ts = GSD_U.trajectory[0]
    assert_almost_equal(ts.dimensions,
                        [ 21.60000038,21.60000038,21.60000038,90.,90.,90.])


def test_gsd_data_step(GSD_U):
    assert GSD_U.trajectory[0].data['step'] == 0
    assert GSD_U.trajectory[1].data['step'] == 500


def test_copy():
    reader = GSDReader(GSD, convert_units=False, dt=2,
                       time_offset=10, foo="bar")

    # test that variables have been allocated properly
    assert reader.convert_units is False
    assert reader._ts_kwargs['dt'] == reader.ts.data['dt'] == 2
    assert reader._ts_kwargs['time_offset'] == 10
    assert reader.ts.data['time_offset'] == 10

    # copy the reader and check that variables are the same
    new_reader = reader.copy()

    assert new_reader.convert_units is False
    assert new_reader._ts_kwargs['dt'] == new_reader.ts.data['dt'] == 2
    assert new_reader._ts_kwargs['time_offset'] == 10
    assert new_reader.ts.data['time_offset'] == 10
    assert new_reader._kwargs['foo'] == 'bar'
    assert new_reader.filename == reader.filename
    # n_atoms does not get passed, but check that it's the same anyways
    assert new_reader.n_atoms == reader.n_atoms
