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

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

import MDAnalysis as mdanalysis
from MDAnalysis.transformations import set_dimensions
from MDAnalysisTests import make_Universe


@pytest.fixture()
def boxdimensions_universe():
    # create Universe objects for tests
    new_u = make_Universe(trajectory=True)
    return new_u


def test_boxdimensions_dims(boxdimensions_universe):
    new_dims = np.float32([2, 2, 2, 90, 90, 90])
    set_dimensions(new_dims)(boxdimensions_universe.trajectory.ts)
    assert_array_almost_equal(boxdimensions_universe.dimensions,
                              new_dims, decimal=6)


@pytest.mark.parametrize('dim_vector_shapes', (
    [1, 1, 1, 90, 90],
    [1, 1, 1, 1, 90, 90, 90],
    np.array([[1], [1], [90], [90], [90]]),
    np.array([1, 1, 1, 90, 90]),
    np.array([1, 1, 1, 1, 90, 90, 90]),
    [1, 1, 1, 90, 90],
    111909090)
    )
def test_dimensions_vector(boxdimensions_universe, dim_vector_shapes):
    # wrong box dimension vector shape
    ts = boxdimensions_universe.trajectory.ts
    with pytest.raises(ValueError, match='valid box dimension shape'):
        set_dimensions(dim_vector_shapes)(ts)


@pytest.mark.parametrize('dim_vector_forms_dtypes', (
    ['a', 'b', 'c', 'd', 'e', 'f'],
    np.array(['a', 'b', 'c', 'd', 'e', 'f']),
    'abcd')
    )
def test_dimensions_vector_asarray(boxdimensions_universe,
                                   dim_vector_forms_dtypes):
    # box dimension input type not convertible into array
    ts = boxdimensions_universe.trajectory.ts
    with pytest.raises(ValueError, match='cannot be converted'):
        set_dimensions(dim_vector_forms_dtypes)(ts)


def test_dimensions_transformations_api(boxdimensions_universe):
    # test if transformation works with on-the-fly transformations API
    new_dims = np.float32([2, 2, 2, 90, 90, 90])
    transform = set_dimensions(new_dims)
    boxdimensions_universe.trajectory.add_transformations(transform)
    for ts in boxdimensions_universe.trajectory:
        assert_array_almost_equal(boxdimensions_universe.dimensions,
                                  new_dims, decimal=6)
