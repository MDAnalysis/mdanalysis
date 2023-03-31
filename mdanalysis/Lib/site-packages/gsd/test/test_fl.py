# Copyright (c) 2016-2023 The Regents of the University of Michigan
# Part of GSD, released under the BSD 2-Clause License.

"""Test gsd.fl."""

import gsd.fl
import gsd.pygsd
import numpy
import platform
import pytest
import random
import pathlib
import os
import shutil
import sys

test_path = pathlib.Path(os.path.realpath(__file__)).parent


def test_create(tmp_path):
    """Test creation of GSD files."""
    gsd.fl.open(mode='xb',
                name=tmp_path / "test_create.gsd",
                application="test_create",
                schema="none",
                schema_version=[1, 2])


@pytest.mark.parametrize('typ', [
    numpy.uint8,
    numpy.uint16,
    numpy.uint32,
    numpy.uint64,
    numpy.int8,
    numpy.int16,
    numpy.int32,
    numpy.int64,
    numpy.float32,
    numpy.float64,
])
def test_dtype(tmp_path, typ):
    """Test all supported data types."""
    data1d = numpy.array([1, 2, 3, 4, 5, 127], dtype=typ)
    data2d = numpy.array([[10, 20], [30, 40], [50, 80]], dtype=typ)
    data_zero = numpy.array([], dtype=typ)

    gsd.fl.open(mode='xb',
                name=tmp_path / "test_dtype.gsd",
                application="test_dtype",
                schema="none",
                schema_version=[1, 2])

    with gsd.fl.open(name=tmp_path / "test_dtype.gsd",
                     mode='wb',
                     application="test_dtype",
                     schema="none",
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='data1d', data=data1d)
        f.write_chunk(name='data2d', data=data2d)
        f.write_chunk(name='data_zero', data=data_zero)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / "test_dtype.gsd",
                     mode='rb',
                     application="test_dtype",
                     schema="none",
                     schema_version=[1, 2]) as f:
        read_data1d = f.read_chunk(frame=0, name='data1d')
        read_data2d = f.read_chunk(frame=0, name='data2d')
        read_data_zero = f.read_chunk(frame=0, name='data_zero')

        assert data1d.dtype == read_data1d.dtype
        numpy.testing.assert_array_equal(data1d, read_data1d)
        assert data2d.dtype == read_data2d.dtype
        numpy.testing.assert_array_equal(data2d, read_data2d)
        assert data_zero.dtype == read_data_zero.dtype
        assert data_zero.shape == (0,)

    # test again with pygsd
    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path / "test_dtype.gsd"), mode='rb')) as f:
        read_data1d = f.read_chunk(frame=0, name='data1d')
        read_data2d = f.read_chunk(frame=0, name='data2d')

        assert data1d.dtype == read_data1d.dtype
        numpy.testing.assert_array_equal(data1d, read_data1d)
        assert data2d.dtype == read_data2d.dtype
        numpy.testing.assert_array_equal(data2d, read_data2d)


def test_metadata(tmp_path, open_mode):
    """Test file metadata."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)

    with gsd.fl.open(name=tmp_path / 'test_metadata.gsd',
                     mode=open_mode.write,
                     application='test_metadata',
                     schema='none',
                     schema_version=[1, 2]) as f:
        assert f.mode == open_mode.write
        for i in range(150):
            f.write_chunk(name='data', data=data)
            f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test_metadata.gsd',
                     mode=open_mode.read,
                     application='test_metadata',
                     schema='none',
                     schema_version=[1, 2]) as f:
        assert f.name == str(tmp_path / 'test_metadata.gsd')
        assert f.mode == open_mode.read
        assert f.application == 'test_metadata'
        assert f.schema == 'none'
        assert f.schema_version == (1, 2)
        assert f.nframes == 150
        assert f.gsd_version == (2, 0)

    # test again with pygsd
    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path / 'test_metadata.gsd'), mode='rb')) as f:
        assert f.name == str(tmp_path / 'test_metadata.gsd')
        assert f.mode == 'rb'
        assert f.application == 'test_metadata'
        assert f.schema == 'none'
        assert f.schema_version == (1, 2)
        assert f.nframes == 150
        assert f.gsd_version == (2, 0)


def test_append(tmp_path, open_mode):
    """Test that data chunks can be appended to existing files."""
    with gsd.fl.open(name=tmp_path / 'test_append.gsd',
                     mode=open_mode.write,
                     application='test_append',
                     schema='none',
                     schema_version=[1, 2]):
        pass

    data = numpy.array([10], dtype=numpy.int64)
    nframes = 1024

    with gsd.fl.open(name=tmp_path / 'test_append.gsd',
                     mode='ab',
                     application='test_append',
                     schema='none',
                     schema_version=[1, 2]) as f:
        assert f.mode == 'ab'
        for i in range(nframes):
            data[0] = i
            f.write_chunk(name='data1', data=data)
            data[0] = i * 10
            f.write_chunk(name='data10', data=data)
            f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test_append.gsd',
                     mode=open_mode.read,
                     application='test_append',
                     schema='none',
                     schema_version=[1, 2]) as f:
        assert f.nframes == nframes
        for i in range(nframes):
            data1 = f.read_chunk(frame=i, name='data1')
            data10 = f.read_chunk(frame=i, name='data10')
            assert data1[0] == i
            assert data10[0] == i * 10

    # test again with pygsd
    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path
                          / 'test_append.gsd'), mode=open_mode.read)) as f:
        assert f.nframes == nframes
        for i in range(nframes):
            data1 = f.read_chunk(frame=i, name='data1')
            data10 = f.read_chunk(frame=i, name='data10')
            assert data1[0] == i
            assert data10[0] == i * 10


def test_chunk_exists(tmp_path, open_mode):
    """Test the chunk_exists API."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)
    with gsd.fl.open(name=tmp_path / 'test_chunk_exists.gsd',
                     mode=open_mode.write,
                     application='test_chunk_exists',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()
        f.write_chunk(name='abcdefg', data=data)
        f.end_frame()
        f.write_chunk(name='test', data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test_chunk_exists.gsd',
                     mode=open_mode.read,
                     application='test_chunk_exists',
                     schema='none',
                     schema_version=[1, 2]) as f:
        assert f.chunk_exists(frame=0, name='chunk1')
        read_data = f.read_chunk(frame=0, name='chunk1')
        assert f.chunk_exists(frame=1, name='abcdefg')
        read_data = f.read_chunk(frame=1, name='abcdefg')
        assert f.chunk_exists(frame=2, name='test')
        read_data = f.read_chunk(frame=2, name='test')

        assert not f.chunk_exists(frame=1, name='chunk1')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=1, name='chunk1')
        assert not f.chunk_exists(frame=2, name='abcdefg')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=2, name='abcdefg')
        assert not f.chunk_exists(frame=0, name='test')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=0, name='test')

        assert not f.chunk_exists(frame=2, name='chunk1')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=2, name='chunk1')
        assert not f.chunk_exists(frame=0, name='abcdefg')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=0, name='abcdefg')
        assert not f.chunk_exists(frame=1, name='test')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=1, name='test')

    # test again with pygsd
    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path / 'test_chunk_exists.gsd'), mode='rb')) as f:
        assert f.chunk_exists(frame=0, name='chunk1')
        read_data = f.read_chunk(frame=0, name='chunk1')
        assert f.chunk_exists(frame=1, name='abcdefg')
        read_data = f.read_chunk(frame=1, name='abcdefg')
        assert f.chunk_exists(frame=2, name='test')
        read_data = f.read_chunk(frame=2, name='test')

        assert not f.chunk_exists(frame=1, name='chunk1')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=1, name='chunk1')
        assert not f.chunk_exists(frame=2, name='abcdefg')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=2, name='abcdefg')
        assert not f.chunk_exists(frame=0, name='test')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=0, name='test')

        assert not f.chunk_exists(frame=2, name='chunk1')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=2, name='chunk1')
        assert not f.chunk_exists(frame=0, name='abcdefg')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=0, name='abcdefg')
        assert not f.chunk_exists(frame=1, name='test')
        with pytest.raises(Exception):
            read_data = f.read_chunk(frame=1, name='test')  # noqa


def test_readonly_errors(tmp_path, open_mode):
    """Test that read only files provide the appropriate errors."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)
    with gsd.fl.open(name=tmp_path / 'test_readonly_errors.gsd',
                     mode=open_mode.write,
                     application='test_readonly_errors',
                     schema='none',
                     schema_version=[1, 2]) as f:
        for i in range(10):
            f.write_chunk(name='chunk1', data=data)
            f.end_frame()

    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)
    with gsd.fl.open(name=tmp_path / 'test_readonly_errors.gsd',
                     mode='rb',
                     application='test_readonly_errors',
                     schema='none',
                     schema_version=[1, 2]) as f:
        with pytest.raises(Exception):
            f.end_frame()

        with pytest.raises(Exception):
            f.write_chunk(name='chunk1', data=data)

    # test again with pygsd
    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path
                          / 'test_readonly_errors.gsd'), mode='rb')) as f:
        with pytest.raises(Exception):
            f.end_frame()

        with pytest.raises(Exception):
            f.write_chunk(name='chunk1', data=data)


def test_fileio_errors(tmp_path, open_mode):
    """Test that OS file I/O errors pass through."""
    # These test cause python to crash on windows....
    if platform.system() != "Windows":
        with pytest.raises(Exception):
            gsd.fl.open(name='/this/file/does/not/exist',
                        application='test_readonly_errors',
                        schema='none',
                        schema_version=[1, 2])

        with open(str(tmp_path / 'test_fileio_errors.gsd'),
                  open_mode.write) as f:
            f.write(b'test')

        with pytest.raises(RuntimeError):
            f = gsd.fl.open(name=tmp_path / 'test_fileio_errors.gsd',
                            mode=open_mode.read,
                            application='test_readonly_errors',
                            schema='none',
                            schema_version=[1, 2])


def test_dtype_errors(tmp_path, open_mode):
    """Test that unsupported data types result in errors."""
    with pytest.raises(Exception):
        data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.bool_)

        with gsd.fl.open(name=tmp_path / 'test_dtype_errors.gsd',
                         mode=open_mode.write,
                         application='test_dtype_errors',
                         schema='none',
                         schema_version=[1, 2]) as f:
            f.write_chunk(name='chunk1', data=data)
            f.end_frame()

    with pytest.raises(Exception):
        data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.float16)

        with gsd.fl.open(name=tmp_path / 'test_dtype_errors.gsd',
                         mode=open_mode.write,
                         application='test_dtype_errors',
                         schema='none',
                         schema_version=[1, 2]) as f:
            f.write_chunk(name='chunk1', data=data)
            f.end_frame()

    with pytest.raises(Exception):
        data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.complex64)

        with gsd.fl.open(name=tmp_path / 'test_dtype_errors.gsd',
                         mode=open_mode.write,
                         application='test_dtype_errors',
                         schema='none',
                         schema_version=[1, 2]) as f:
            f.write_chunk(name='chunk1', data=data)
            f.end_frame()

    with pytest.raises(Exception):
        data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.complex128)

        with gsd.fl.open(name=tmp_path / 'test_dtype_errors.gsd',
                         mode=open_mode.write,
                         application='test_dtype_errors',
                         schema='none',
                         schema_version=[1, 2]) as f:
            f.write_chunk(name='chunk1', data=data)
            f.end_frame()


def test_truncate(tmp_path):
    """Test that the truncate method functions."""
    data = numpy.ascontiguousarray(numpy.random.random(size=(1000, 3)),
                                   dtype=numpy.float32)
    with gsd.fl.open(name=tmp_path / 'test_truncate.gsd',
                     mode='wb',
                     application='test_truncate',
                     schema='none',
                     schema_version=[1, 2]) as f:
        assert f.mode == 'wb'
        for i in range(10):
            f.write_chunk(name='data', data=data)
            f.end_frame()

        assert f.nframes == 10

        f.truncate()
        assert f.nframes == 0
        assert f.application == 'test_truncate'
        assert f.schema == 'none'
        assert f.schema_version == (1, 2)

        f.write_chunk(name='data', data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test_truncate.gsd',
                     mode='rb',
                     application='test_truncate',
                     schema='none',
                     schema_version=[1, 2]) as f:
        assert f.name == str(tmp_path / 'test_truncate.gsd')
        assert f.mode == 'rb'
        assert f.application == 'test_truncate'
        assert f.schema == 'none'
        assert f.schema_version == (1, 2)
        assert f.nframes == 1


def test_namelen(tmp_path, open_mode):
    """Test that long names are truncated as documented."""
    app_long = 'abcdefga' * 100
    schema_long = 'ijklmnop' * 100
    chunk_long = '12345678' * 100

    with gsd.fl.open(name=tmp_path / 'test_namelen.gsd',
                     mode=open_mode.write,
                     application=app_long,
                     schema=schema_long,
                     schema_version=[1, 2]) as f:
        assert f.application == app_long[0:63]
        assert f.schema == schema_long[0:63]

        data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)
        f.write_chunk(name=chunk_long, data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test_namelen.gsd',
                     mode=open_mode.read,
                     application=app_long,
                     schema=schema_long,
                     schema_version=[1, 2]) as f:
        data_read = f.read_chunk(0, name=chunk_long)
        numpy.testing.assert_array_equal(data, data_read)

    # test again with pygsd
    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path / 'test_namelen.gsd'), mode='rb')) as f:
        data_read = f.read_chunk(0, name=chunk_long)
        numpy.testing.assert_array_equal(data, data_read)


def test_open(tmp_path):
    """Test the open() API."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode='xb',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test_2.gsd',
                     mode='xb+',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()
        f.read_chunk(0, name='chunk1')

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode='wb',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode='wb+',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()
        f.read_chunk(0, name='chunk1')

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode='ab',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode='rb',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.read_chunk(0, name='chunk1')
        f.read_chunk(1, name='chunk1')

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode='rb+',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()
        f.read_chunk(0, name='chunk1')
        f.read_chunk(1, name='chunk1')
        f.read_chunk(2, name='chunk1')


def test_find_matching_chunk_names(tmp_path, open_mode):
    """Test the find_matching_chunk_names API."""
    data = numpy.array([1, 2, 3, 4, 5], dtype=numpy.float32)

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode=open_mode.write,
                     application='test_find_matching_chunk_names',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='log/A', data=data)
        f.write_chunk(name='log/chunk2', data=data)
        f.end_frame()
        f.write_chunk(name='data/B', data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode=open_mode.read,
                     application='test_find_matching_chunk_names',
                     schema='none',
                     schema_version=[1, 2]) as f:
        all_chunks = f.find_matching_chunk_names('')
        assert len(all_chunks) == 3
        assert 'log/A' in all_chunks
        assert 'log/chunk2' in all_chunks
        assert 'data/B' in all_chunks

        log_chunks = f.find_matching_chunk_names('log/')
        assert len(log_chunks) == 2
        assert 'log/A' in log_chunks
        assert 'log/chunk2' in log_chunks

        data_chunks = f.find_matching_chunk_names('data/')
        assert len(data_chunks) == 1
        assert 'data/B' in data_chunks

        other_chunks = f.find_matching_chunk_names('other/')
        assert len(other_chunks) == 0

    # test again with pygsd
    with gsd.pygsd.GSDFile(file=open(str(tmp_path
                                         / "test.gsd"), mode='rb')) as f:
        all_chunks = f.find_matching_chunk_names('')
        assert len(all_chunks) == 3
        assert 'log/A' in all_chunks
        assert 'log/chunk2' in all_chunks
        assert 'data/B' in all_chunks

        log_chunks = f.find_matching_chunk_names('log/')
        assert len(log_chunks) == 2
        assert 'log/A' in log_chunks
        assert 'log/chunk2' in log_chunks

        data_chunks = f.find_matching_chunk_names('data/')
        assert len(data_chunks) == 1
        assert 'data/B' in data_chunks

        other_chunks = f.find_matching_chunk_names('other/')
        assert len(other_chunks) == 0


def test_chunk_name_limit(tmp_path, open_mode):
    """Test that providing more than the maximum allowed chunk names errors."""
    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode=open_mode.write,
                     application='test_chunk_name_limit',
                     schema='none',
                     schema_version=[1, 2]) as f:
        for i in range(65535):
            f.write_chunk(name=str(i), data=numpy.array([i], dtype=numpy.int32))

        # The GSD specification limits to 65535 names:
        with pytest.raises(RuntimeError):
            f.write_chunk(name='65536',
                          data=numpy.array([i], dtype=numpy.int32))


def test_many_names(tmp_path, open_mode):
    """Test that many chunk names can be written to a file."""
    values = list(range(1000))

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode=open_mode.write,
                     application='test_many_names',
                     schema='none',
                     schema_version=[1, 2]) as f:
        for frame in range(5):
            random.shuffle(values)
            for value in values:
                f.write_chunk(name=str(value),
                              data=numpy.array([value * 13], dtype=numpy.int32))
            f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test.gsd',
                     mode=open_mode.read,
                     application='test_many_names',
                     schema='none',
                     schema_version=[1, 2]) as f:

        for frame in range(5):
            random.shuffle(values)
            for value in values:
                data = numpy.array([value * 13], dtype=numpy.int32)
                data_read = f.read_chunk(frame=frame, name=str(value))
                numpy.testing.assert_array_equal(data, data_read)

    with gsd.pygsd.GSDFile(file=open(str(tmp_path
                                         / 'test.gsd'), mode='rb')) as f:
        for frame in range(5):
            random.shuffle(values)
            for value in values:
                data = numpy.array([value * 13], dtype=numpy.int32)
                data_read = f.read_chunk(frame=frame, name=str(value))
                numpy.testing.assert_array_equal(data, data_read)


def test_gsd_v1_read():
    """Test that the GSD v2 API can read v1 files."""
    values = list(range(127))
    values_str = [str(v) for v in values]
    values_str.sort()

    # test that we can:
    # 1) Read chunk values correctly
    # 2) Iterate through chunk names correctly
    def check_v1_file_read(f):
        assert f.gsd_version == (1, 0)

        for frame in range(5):
            random.shuffle(values)
            for value in values:
                data = numpy.array([value * 13], dtype=numpy.int32)
                data_read = f.read_chunk(frame=frame, name=str(value))
                numpy.testing.assert_array_equal(data, data_read)

        chunk_names = f.find_matching_chunk_names('')
        chunk_names.sort()
        assert chunk_names == values_str

    # test with the C implemantation
    with gsd.fl.open(name=test_path / 'test_gsd_v1.gsd',
                     mode='rb',
                     application='test_gsd_v1',
                     schema='none',
                     schema_version=[1, 2]) as f:

        check_v1_file_read(f)

    # and the pure Python implementation
    with gsd.pygsd.GSDFile(
            file=open(str(test_path / 'test_gsd_v1.gsd'), mode='rb')) as f:

        assert f.gsd_version == (1, 0)

        check_v1_file_read(f)


def test_gsd_v1_upgrade_read(tmp_path, open_mode):
    """Test that v1 files can be upgraded to v2."""
    values = list(range(127))
    values_str = [str(v) for v in values]
    values_str.sort()

    # test that we can:
    # 1) Read chunk values correctly
    # 2) Iterate through chunk names correctly
    def check_v1_file_read(f):
        for frame in range(5):
            random.shuffle(values)
            for value in values:
                data = numpy.array([value * 13], dtype=numpy.int32)
                data_read = f.read_chunk(frame=frame, name=str(value))
                numpy.testing.assert_array_equal(data, data_read)

        chunk_names = f.find_matching_chunk_names('')
        chunk_names.sort()
        assert chunk_names == values_str

    shutil.copy(test_path / 'test_gsd_v1.gsd', tmp_path / 'test_gsd_v1.gsd')

    with gsd.fl.open(name=tmp_path / 'test_gsd_v1.gsd',
                     mode='rb+',
                     application='test_gsd_v1',
                     schema='none',
                     schema_version=[1, 2]) as f:

        assert f.gsd_version == (1, 0)

        f.upgrade()

        # check that we can read the file contents after the upgrade in memory
        check_v1_file_read(f)

    # and the same tests again after closing and opening the file
    with gsd.fl.open(name=tmp_path / 'test_gsd_v1.gsd',
                     mode=open_mode.read,
                     application='test_gsd_v1',
                     schema='none',
                     schema_version=[1, 2]) as f:

        assert f.gsd_version == (2, 0)

        check_v1_file_read(f)

    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path / 'test_gsd_v1.gsd'), mode='rb')) as f:

        assert f.gsd_version == (2, 0)

        check_v1_file_read(f)


def test_gsd_v1_write(tmp_path, open_mode):
    """Test that v2 can write to v1 files."""
    values = list(range(256))
    # include a very long chunk name to check that the name is truncated
    # properly for the v1 format limitations
    long_name = 'abcdefg' * 1000
    values.append(long_name)

    values_str = []
    for v in values:
        if type(v) == str and len(v) > 63:
            # v1 files truncate names to 63 chars
            v = v[0:63]
        values_str.append(str(v))
    values_str.sort()

    shutil.copy(test_path / 'test_gsd_v1.gsd', tmp_path / 'test_gsd_v1.gsd')

    # test that we can:
    # 1) Read chunk values correctly
    # 2) Iterate through chunk names correctly
    def check_v1_file_read(f):
        assert f.gsd_version == (1, 0)

        chunk_names = f.find_matching_chunk_names('')
        chunk_names.sort()
        assert chunk_names == values_str

        frame = 5
        random.shuffle(values)
        for value in values:
            if type(value) == int:
                data = numpy.array([value * 13], dtype=numpy.int32)
            else:
                data = numpy.array([hash(value)], dtype=numpy.int64)
                # v1 files truncate names to 63 chars
                if len(value) > 63:
                    value = value[0:63]

            data_read = f.read_chunk(frame=frame, name=str(value))
            numpy.testing.assert_array_equal(data, data_read)

    # test that we can write new entries to the file
    with gsd.fl.open(name=tmp_path / 'test_gsd_v1.gsd',
                     mode='rb+',
                     application='test_gsd_v1',
                     schema='none',
                     schema_version=[1, 2]) as f:

        assert f.gsd_version == (1, 0)

        for value in values:
            if type(value) == int:
                data = numpy.array([value * 13], dtype=numpy.int32)
            else:
                data = numpy.array([hash(value)], dtype=numpy.int64)
            f.write_chunk(name=str(value), data=data)
        f.end_frame()

        check_v1_file_read(f)

    # test opening again with the C implemantation
    with gsd.fl.open(name=tmp_path / 'test_gsd_v1.gsd',
                     mode=open_mode.read,
                     application='test_gsd_v1',
                     schema='none',
                     schema_version=[1, 2]) as f:

        check_v1_file_read(f)

    # and the pure Python implementation
    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path / 'test_gsd_v1.gsd'), mode='rb')) as f:

        assert f.gsd_version == (1, 0)

        check_v1_file_read(f)


def test_gsd_v1_upgrade_write(tmp_path, open_mode):
    """Test that upgraded files can be written to after upgraded."""
    values = list(range(256))
    # include a very long chunk name to check that the name can be written
    # after the upgrade
    long_name = 'abcdefg' * 1000
    values.append(long_name)

    values_str = [str(v) for v in values]
    values_str.sort()

    shutil.copy(test_path / 'test_gsd_v1.gsd', tmp_path / 'test_gsd_v1.gsd')

    # test that we can:
    # 1) Read chunk values correctly
    # 2) Iterate through chunk names correctly
    def check_v1_file_read(f):
        chunk_names = f.find_matching_chunk_names('')
        chunk_names.sort()
        assert chunk_names == values_str

        frame = 5
        random.shuffle(values)
        for value in values:
            if type(value) == int:
                data = numpy.array([value * 13], dtype=numpy.int32)
            else:
                data = numpy.array([hash(value)], dtype=numpy.int64)

            data_read = f.read_chunk(frame=frame, name=str(value))
            numpy.testing.assert_array_equal(data, data_read)

    # test that we can write new entries to the file
    with gsd.fl.open(name=tmp_path / 'test_gsd_v1.gsd',
                     mode='rb+',
                     application='test_gsd_v1',
                     schema='none',
                     schema_version=[1, 2]) as f:

        assert f.gsd_version == (1, 0)

        f.upgrade()

        assert f.gsd_version == (2, 0)

        for value in values:
            if type(value) == int:
                data = numpy.array([value * 13], dtype=numpy.int32)
            else:
                data = numpy.array([hash(value)], dtype=numpy.int64)
            f.write_chunk(name=str(value), data=data)
        f.end_frame()

        check_v1_file_read(f)

    # test opening again with the C implemantation
    with gsd.fl.open(name=tmp_path / 'test_gsd_v1.gsd',
                     mode=open_mode.read,
                     application='test_gsd_v1',
                     schema='none',
                     schema_version=[1, 2]) as f:

        assert f.gsd_version == (2, 0)

        check_v1_file_read(f)

    # and the pure Python implementation
    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path / 'test_gsd_v1.gsd'), mode='rb')) as f:

        assert f.gsd_version == (2, 0)

        check_v1_file_read(f)


def test_zero_size(tmp_path, open_mode):
    """Test that zero-size data chunks are allowed."""
    data = numpy.array([], dtype=numpy.float32)

    with gsd.fl.open(name=tmp_path / 'test_zero.gsd',
                     mode=open_mode.write,
                     application='test_zero',
                     schema='none',
                     schema_version=[1, 2]) as f:

        f.write_chunk(name='data', data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test_zero.gsd',
                     mode=open_mode.read,
                     application='test_zero',
                     schema='none',
                     schema_version=[1, 2]) as f:
        assert f.nframes == 1
        data_read = f.read_chunk(frame=0, name='data')
        assert data_read.shape == (0,)
        assert data_read.dtype == numpy.float32

    # test again with pygsd
    with gsd.pygsd.GSDFile(
            file=open(str(tmp_path
                          / 'test_zero.gsd'), mode=open_mode.read)) as f:
        assert f.nframes == 1
        data_read = f.read_chunk(frame=0, name='data')
        assert data_read.shape == (0,)
        assert data_read.dtype == numpy.float32


@pytest.mark.skipif(sys.version_info < (3, 7),
                    reason="Python 3.6 fails to handle non-ascii characters.")
def test_utf8(tmp_path):
    """Test that the API handles UTF-8 encoding for the filename."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)

    fname = '中文.gsd'

    with gsd.fl.open(name=tmp_path / fname,
                     mode='xb',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    dir_list = os.listdir(tmp_path)
    print(dir_list)
    assert fname in dir_list

    with gsd.fl.open(name=tmp_path / fname,
                     mode='wb',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / fname,
                     mode='rb',
                     application='test_open',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.read_chunk(0, name='chunk1')


@pytest.mark.parametrize('mode', ['wb', 'wb+', 'ab', 'rb+', 'xb', 'xb+'])
def test_read_write(tmp_path, mode):
    """Test that data chunks can read from files opened in all write modes."""
    if mode[0] == 'r' or mode[0] == 'a':
        with gsd.fl.open(name=tmp_path / 'test_read_write.gsd',
                         mode='wb',
                         application='test_read_write',
                         schema='none',
                         schema_version=[1, 2]):
            pass

    data = numpy.array([10], dtype=numpy.int64)
    nframes = 1024

    with gsd.fl.open(name=tmp_path / 'test_read_write.gsd',
                     mode=mode,
                     application='test_read_write',
                     schema='none',
                     schema_version=[1, 2]) as f:
        assert f.mode == mode
        for i in range(nframes):
            data[0] = i
            f.write_chunk(name='data1', data=data)
            data[0] = i * 10
            f.write_chunk(name='data10', data=data)
            f.end_frame()

        for i in range(nframes):
            data1 = f.read_chunk(frame=i, name='data1')
            data10 = f.read_chunk(frame=i, name='data10')
            assert data1[0] == i
            assert data10[0] == i * 10
