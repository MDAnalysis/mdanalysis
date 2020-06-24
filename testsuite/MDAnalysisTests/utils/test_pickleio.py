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
import pickle

import pytest
from numpy.testing import assert_equal

from MDAnalysis.lib.picklable_file_io import pickle_open, FileIOPicklable
from MDAnalysis.tests.datafiles import PDB


@pytest.fixture(params=[
    # filename mode
    (PDB, 'r'),
    (PDB, 'rt'),
    (PDB, 'rb'),
])
def f(request):
    filename, mode = request.param
    return pickle_open(filename, mode)


def test_iopickle(f):
    f_pickled = pickle.loads(pickle.dumps(f))
    assert_equal(f.readline(), f_pickled.readline())


def test_offset(f):
    f.readline()
    f_pickled = pickle.loads(pickle.dumps(f))
    assert_equal(f.tell(), f_pickled.tell())


def test_context_manager_pickle():
    with pickle_open(PDB) as file:
        file_pickled = pickle.loads(pickle.dumps(file))
        assert_equal(file.readline(), file_pickled.readline())


def test_fileio_pickle():
    raw_io = FileIOPicklable(PDB)
    raw_io_pickled = pickle.loads(pickle.dumps(raw_io))
    assert_equal(raw_io.readline(), raw_io_pickled.readline())


@pytest.fixture(params=[
    # filename mode
    (PDB, 'w'),
    (PDB, 'x'),
    (PDB, 'a'),
])
def unpicklable_f(request):
    filename, mode = request.param
    return filename, mode


def test_unpicklable_open_mode(unpicklable_f):
    filename, mode = unpicklable_f
    with pytest.raises(ValueError, match=r"Only read mode *"):
        pickle_open(filename, mode)


def test_pickle_seek_fallback():
    with pickle_open(PDB) as file:
        file.__next__()
        file_pickled = pickle.loads(pickle.dumps(file))
        assert_equal(file_pickled.tell(), 0)
