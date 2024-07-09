
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
import pytest
from numpy.testing import assert_equal
from threadpoolctl import threadpool_info

import MDAnalysis as mda
from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysis.transformations.base import TransformationBase


class DefaultTransformation(TransformationBase):
    """Default values for max_threads and parallelizable"""
    def __init__(self):
        super().__init__()

    def _transform(self, ts):
        self.runtime_info = threadpool_info()
        ts.positions = ts.positions + 1
        return ts


class NoTransform_Transformation(TransformationBase):
    """Default values for max_threads and parallelizable"""
    def __init__(self):
        super().__init__()


class CustomTransformation(TransformationBase):
    """Custom value for max_threads and parallelizable"""
    def __init__(self, max_threads=1, parallelizable=False):
        super().__init__(max_threads=max_threads,
                         parallelizable=parallelizable)

    def _transform(self, ts):
        self.runtime_info = threadpool_info()
        ts.positions = ts.positions + 1
        return ts


@pytest.fixture(scope='module')
def u():
    return mda.Universe(PSF, DCD)


def test_default_value():
    new_trans = DefaultTransformation()
    assert new_trans.max_threads is None
    assert new_trans.parallelizable is True


def test_no_transform_function(u):
    new_trans = NoTransform_Transformation()
    with pytest.raises(NotImplementedError, match=r"Only implemented"):
        _ = new_trans._transform(u.trajectory.ts)


def test_custom_value():
    new_trans = CustomTransformation()
    assert new_trans.max_threads == 1
    assert new_trans.parallelizable is False


def test_setting_thread_limit_value():
    new_trans = CustomTransformation(max_threads=4)
    assert new_trans.max_threads == 4


def test_thread_limit_apply(u):
    default_thread_info = threadpool_info()
    default_num_thread_limit_list = [thread_info['num_threads']
                                     for thread_info in default_thread_info]

    new_trans = CustomTransformation(max_threads=2)
    _ = new_trans(u.trajectory.ts)
    for thread_info in new_trans.runtime_info:
        assert thread_info['num_threads'] == 2

    #  test the thread limit is only applied locally.
    new_thread_info = threadpool_info()
    new_num_thread_limit_list = [thread_info['num_threads']
                                 for thread_info in new_thread_info]
    assert_equal(default_num_thread_limit_list, new_num_thread_limit_list)
