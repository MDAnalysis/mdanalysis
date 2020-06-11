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

import pytest
from numpy.testing import assert_equal

import MDAnalysis.lib.util as util
from MDAnalysis.tests.datafiles import PDB

import pickle
import sys


@pytest.fixture(params=[
    # filename mode
    (PDB, 'r'),
    (PDB, 'rt'),
    (PDB, 'rb'),
])
def f(request):
    filename, mode = request.param
    return util.pickle_open(filename, mode)


@pytest.mark.xfail(sys.version_info < (3, 0), reason="pickle function not \
                     working in python 2")
def test_iopickle(f):
    f_pickled = pickle.loads(pickle.dumps(f))
    assert_equal(f.readline(), f_pickled.readline())


@pytest.mark.xfail(sys.version_info < (3, 0), reason="pickle function not \
                     working in python 2")
def test_offset(f):
    f.readline()
    f_pickled = pickle.loads(pickle.dumps(f))
    assert_equal(f.tell(), f_pickled.tell())
