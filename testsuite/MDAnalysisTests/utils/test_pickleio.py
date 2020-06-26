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

from MDAnalysis.lib.util import anyopen
from MDAnalysis.lib.picklable_file_io import (
    pickle_open,
    BufferIOPicklable,
    FileIOPicklable,
    TextIOPicklable,
    BZ2Picklable,
    GzipPicklable
)
from MDAnalysis.tests.datafiles import (
    PDB,
    XYZ_bz2,
    MMTF_gz,
    GMS_ASYMOPT
)


@pytest.fixture(params=[
    # filename mode
    (PDB, 'r'),
    (PDB, 'rt'),
    (XYZ_bz2, 'rt'),
    (GMS_ASYMOPT, 'rt')
])
def f_text(request):
    filename, mode = request.param
    return anyopen(filename, mode)


def test_get_right_open_handler_text(f_text):
    assert_equal(f_text.__class__, TextIOPicklable)


def test_iopickle_text(f_text):
    f_text_pickled = pickle.loads(pickle.dumps(f_text))
    assert_equal(f_text.readline(), f_text_pickled.readline())


def test_offset_text_to_0(f_text):
    f_text.readline()
    f_text_pickled = pickle.loads(pickle.dumps(f_text))
    assert_equal(f_text_pickled.tell(), 0)


@pytest.fixture(params=[
    # filename mode ref_class
    (PDB, 'rb', BufferIOPicklable),
    (XYZ_bz2, 'rb', BZ2Picklable),
    (MMTF_gz, 'rb', GzipPicklable)
])
def f_byte(request):
    filename, mode, ref_reader_class = request.param
    return (anyopen(filename, mode), ref_reader_class)


def test_get_right_open_handler_byte(f_byte):
    assert_equal(f_byte[0].__class__, f_byte[1])


def test_iopickle_byte(f_byte):
    f_byte_pickled = pickle.loads(pickle.dumps(f_byte[0]))
    assert_equal(f_byte[0].readline(), f_byte_pickled.readline())


def test_offset_byte_to_tell(f_byte):
    f_byte[0].readline()
    f_byte_pickled = pickle.loads(pickle.dumps(f_byte[0]))
    assert_equal(f_byte_pickled.tell(), f_byte[0].tell())


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
