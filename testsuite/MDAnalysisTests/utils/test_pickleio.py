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
    BufferIOPicklable,
    FileIOPicklable,
    TextIOPicklable,
    BZ2Picklable,
    GzipPicklable,
    pickle_open,
    bz2_pickle_open,
    gzip_pickle_open,
)
from MDAnalysis.coordinates.GSD import (
    gsd_pickle_open,
    HAS_GSD
)
from MDAnalysis.coordinates.TRJ import (
    NCDFPicklable,
)
from MDAnalysis.coordinates.chemfiles import (
    check_chemfiles_version
)
if check_chemfiles_version():
    from MDAnalysis.coordinates.chemfiles import (
        ChemfilesPicklable
    )
from MDAnalysis.coordinates.H5MD import HAS_H5PY
from MDAnalysis.coordinates.H5MD import (
    H5PYPicklable
)

from MDAnalysis.tests.datafiles import (
    PDB,
    XYZ,
    XYZ_bz2,
    MMTF_gz,
    GMS_ASYMOPT,
    GSD,
    NCDF,
    H5MD_xvf
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
    assert_equal(f_text.readlines(), f_text_pickled.readlines())


def test_offset_text_same(f_text):
    f_text.readline()
    f_text_pickled = pickle.loads(pickle.dumps(f_text))
    assert_equal(f_text_pickled.tell(), f_text.tell())


@pytest.fixture(params=[
    # filename mode ref_class
    (PDB, 'rb', BufferIOPicklable),
    (XYZ_bz2, 'rb', BZ2Picklable),
    (MMTF_gz, 'rb', GzipPicklable)
])
def f_byte(request):
    filename, mode, ref_reader_class = request.param
    return anyopen(filename, mode), ref_reader_class


def test_get_right_open_handler_byte(f_byte):
    assert_equal(f_byte[0].__class__, f_byte[1])


def test_iopickle_byte(f_byte):
    file = f_byte[0]
    f_byte_pickled = pickle.loads(pickle.dumps(file))
    assert_equal(file.readlines(), f_byte_pickled.readlines())


def test_offset_byte_to_tell(f_byte):
    file = f_byte[0]
    file.readline()
    f_byte_pickled = pickle.loads(pickle.dumps(file))
    assert_equal(f_byte_pickled.tell(), file.tell())


def test_context_manager_pickle():
    with pickle_open(PDB) as file:
        file_pickled = pickle.loads(pickle.dumps(file))
        assert_equal(file.readlines(), file_pickled.readlines())


def test_fileio_pickle():
    raw_io = FileIOPicklable(PDB)
    raw_io_pickled = pickle.loads(pickle.dumps(raw_io))
    assert_equal(raw_io.readlines(), raw_io_pickled.readlines())


@pytest.fixture(params=[
    # filename mode open_func open_class
    ('test.pdb', 'w', pickle_open, FileIOPicklable),
    ('test.pdb', 'x', pickle_open, FileIOPicklable),
    ('test.pdb', 'a', pickle_open, FileIOPicklable),
    ('test.bz2', 'w', bz2_pickle_open, BZ2Picklable),
    ('test.gz', 'w', gzip_pickle_open, GzipPicklable),
])
def unpicklable_f(request):
    filename, mode, open_func, open_class = request.param
    return filename, mode, open_func, open_class


def test_unpicklable_open_mode(unpicklable_f, tmpdir):
    filename, mode, open_func, open_class = unpicklable_f
    with pytest.raises(ValueError, match=r"Only read mode"):
        open_func(tmpdir.mkdir("pickle").join(filename), mode)


def test_pickle_with_write_mode(unpicklable_f, tmpdir):
    filename, mode, open_func, open_class = unpicklable_f
    f_open_by_class = open_class(tmpdir.mkdir("pickle").join(filename), mode)
    with pytest.raises(RuntimeError, match=r"Can only pickle"):
        f_pickled = pickle.loads(pickle.dumps(f_open_by_class))


@pytest.mark.skipif(not HAS_GSD, reason='gsd not installed')
def test_GSD_pickle():
    gsd_io = gsd_pickle_open(GSD, mode='r')
    gsd_io_pickled = pickle.loads(pickle.dumps(gsd_io))
    assert_equal(gsd_io[0].particles.position,
                 gsd_io_pickled[0].particles.position)


@pytest.mark.skipif(not HAS_GSD, reason='gsd not installed')
def test_GSD_with_write_mode(tmpdir):
    with pytest.raises(ValueError, match=r"Only read mode"):
        gsd_io = gsd_pickle_open(tmpdir.mkdir("gsd").join('t.gsd'),
                                 mode='w')


def test_NCDF_pickle():
    ncdf_io = NCDFPicklable(NCDF, mmap=None)
    ncdf_io_pickled = pickle.loads(pickle.dumps(ncdf_io))
    assert_equal(ncdf_io.variables['coordinates'][0],
                 ncdf_io_pickled.variables['coordinates'][0])


def test_NCDF_mmap_pickle():
    ncdf_io = NCDFPicklable(NCDF, mmap=False)
    ncdf_io_pickled = pickle.loads(pickle.dumps(ncdf_io))
    assert_equal(ncdf_io_pickled.use_mmap, False)


@pytest.mark.skipif(not check_chemfiles_version(),
                    reason="Wrong version of chemfiles")
def test_Chemfiles_pickle():
    chemfiles_io = ChemfilesPicklable(XYZ)
    chemfiles_io_pickled = pickle.loads(pickle.dumps(chemfiles_io))
    #  frame has to be first saved to get the right position value.
    #  As opposed to `chemfiles_io.read().positions)
    frame = chemfiles_io.read()
    frame_pickled = chemfiles_io_pickled.read()
    assert_equal(frame.positions[:],
                 frame_pickled.positions[:])


@pytest.mark.skipif(not check_chemfiles_version(),
                    reason="Wrong version of chemfiles")
def test_Chemfiles_with_write_mode(tmpdir):
    with pytest.raises(ValueError, match=r"Only read mode"):
        chemfiles_io = ChemfilesPicklable(tmpdir.mkdir("xyz").join('t.xyz'),
                                          mode='w')


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_H5MD_pickle():
    h5md_io = H5PYPicklable(H5MD_xvf, 'r')
    h5md_io_pickled = pickle.loads(pickle.dumps(h5md_io))
    assert_equal(h5md_io['particles/trajectory/position/value'][0],
                 h5md_io_pickled['particles/trajectory/position/value'][0])


@pytest.mark.skipif(not HAS_H5PY, reason="h5py not installed")
def test_H5MD_pickle_with_driver():
    h5md_io = H5PYPicklable(H5MD_xvf, 'r', driver='core')
    h5md_io_pickled = pickle.loads(pickle.dumps(h5md_io))
    assert_equal(h5md_io['particles/trajectory/position/value'][0],
                 h5md_io_pickled['particles/trajectory/position/value'][0])
