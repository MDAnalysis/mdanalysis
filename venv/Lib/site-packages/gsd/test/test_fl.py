# Copyright (c) 2016-2025 The Regents of the University of Michigan
# Part of GSD, released under the BSD 2-Clause License.

"""Test gsd.fl."""

import os
import pathlib
import platform
import random
import shutil
import sys

import gsd.fl
import numpy
import pytest

import gsd.pygsd

test_path = pathlib.Path(os.path.realpath(__file__)).parent
current_gsd_version = (2, 1)


def test_create(tmp_path, open_mode):
    """Test creation of GSD files."""
    gsd.fl.open(
        mode=open_mode.write,
        name=tmp_path / 'test_create.gsd',
        application='test_create',
        schema='none',
        schema_version=[1, 2],
    )


@pytest.mark.parametrize(
    'typ',
    [
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
    ],
)
def test_nonstring_dtypes(tmp_path, typ):
    """Test all supported data types except for strings."""
    data1d = numpy.array([1, 2, 3, 4, 5, 127], dtype=typ)
    data2d = numpy.array([[10, 20], [30, 40], [50, 80]], dtype=typ)
    data_zero = numpy.array([], dtype=typ)

    gsd.fl.open(
        mode='x',
        name=tmp_path / 'test_dtype.gsd',
        application='test_dtype',
        schema='none',
        schema_version=[1, 2],
    )

    with gsd.fl.open(
        name=tmp_path / 'test_dtype.gsd',
        mode='w',
        application='test_dtype',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='data1d', data=data1d)
        f.write_chunk(name='data2d', data=data2d)
        f.write_chunk(name='data_zero', data=data_zero)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test_dtype.gsd',
        mode='r',
        application='test_dtype',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        read_data1d = f.read_chunk(frame=0, name='data1d')
        read_data2d = f.read_chunk(frame=0, name='data2d')
        read_data_zero = f.read_chunk(frame=0, name='data_zero')

        assert data1d.dtype.type == read_data1d.dtype.type
        numpy.testing.assert_array_equal(data1d, read_data1d)
        assert data2d.dtype.type == read_data2d.dtype.type
        numpy.testing.assert_array_equal(data2d, read_data2d)
        assert data_zero.dtype.type == read_data_zero.dtype.type
        assert data_zero.shape == (0,)

    # test again with pygsd
    with gsd.pygsd.GSDFile(file=open(str(tmp_path / 'test_dtype.gsd'), mode='rb')) as f:
        read_data1d = f.read_chunk(frame=0, name='data1d')
        read_data2d = f.read_chunk(frame=0, name='data2d')

        assert data1d.dtype.type == read_data1d.dtype.type
        numpy.testing.assert_array_equal(data1d, read_data1d)
        assert data2d.dtype.type == read_data2d.dtype.type
        numpy.testing.assert_array_equal(data2d, read_data2d)


def test_string_dtype(tmp_path):
    """Test string datatype.

    Note that the string datatype does not support 0-D or 2-D data.
    """
    data1d = 'test'

    gsd.fl.open(
        mode='x',
        name=tmp_path / 'test_dtype.gsd',
        application='test_dtype',
        schema='none',
        schema_version=[1, 2],
    )

    with gsd.fl.open(
        name=tmp_path / 'test_dtype.gsd',
        mode='w',
        application='test_dtype',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='data1d', data=data1d)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test_dtype.gsd',
        mode='r',
        application='test_dtype',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        read_data1d = f.read_chunk(frame=0, name='data1d')

        numpy.testing.assert_string_equal(data1d, read_data1d)

    # test again with pygsd
    with gsd.pygsd.GSDFile(file=open(str(tmp_path / 'test_dtype.gsd'), mode='rb')) as f:
        read_data1d = f.read_chunk(frame=0, name='data1d')

        numpy.testing.assert_string_equal(data1d, read_data1d)


def test_metadata(tmp_path, open_mode):
    """Test file metadata."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)

    with gsd.fl.open(
        name=tmp_path / 'test_metadata.gsd',
        mode=open_mode.write,
        application='test_metadata',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.mode == open_mode.write
        for _i in range(150):
            f.write_chunk(name='data', data=data)
            f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test_metadata.gsd',
        mode=open_mode.read,
        application='test_metadata',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.name == str(tmp_path / 'test_metadata.gsd')
        assert f.mode == open_mode.read
        assert f.application == 'test_metadata'
        assert f.schema == 'none'
        assert f.schema_version == (1, 2)
        assert f.nframes == 150
        assert f.gsd_version == current_gsd_version

    # test again with pygsd
    with gsd.pygsd.GSDFile(
        file=open(str(tmp_path / 'test_metadata.gsd'), mode='rb')
    ) as f:
        assert f.name == str(tmp_path / 'test_metadata.gsd')
        assert f.mode == 'r'
        assert f.application == 'test_metadata'
        assert f.schema == 'none'
        assert f.schema_version == (1, 2)
        assert f.nframes == 150
        assert f.gsd_version == current_gsd_version


def test_append(tmp_path, open_mode):
    """Test that data chunks can be appended to existing files."""
    with gsd.fl.open(
        name=tmp_path / 'test_append.gsd',
        mode=open_mode.write,
        application='test_append',
        schema='none',
        schema_version=[1, 2],
    ):
        pass

    data = numpy.array([10], dtype=numpy.int64)
    nframes = 1024

    with gsd.fl.open(
        name=tmp_path / 'test_append.gsd',
        mode='a',
        application='test_append',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.mode == 'a'
        for i in range(nframes):
            data[0] = i
            f.write_chunk(name='data1', data=data)
            data[0] = i * 10
            f.write_chunk(name='data10', data=data)
            f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test_append.gsd',
        mode=open_mode.read,
        application='test_append',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.nframes == nframes
        for i in range(nframes):
            data1 = f.read_chunk(frame=i, name='data1')
            data10 = f.read_chunk(frame=i, name='data10')
            assert data1[0] == i
            assert data10[0] == i * 10

    # test again with pygsd
    with gsd.pygsd.GSDFile(
        file=open(str(tmp_path / 'test_append.gsd'), mode='rb')
    ) as f:
        assert f.nframes == nframes
        for i in range(nframes):
            data1 = f.read_chunk(frame=i, name='data1')
            data10 = f.read_chunk(frame=i, name='data10')
            assert data1[0] == i
            assert data10[0] == i * 10


def test_chunk_exists(tmp_path, open_mode):
    """Test the chunk_exists API."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)
    with gsd.fl.open(
        name=tmp_path / 'test_chunk_exists.gsd',
        mode=open_mode.write,
        application='test_chunk_exists',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()
        f.write_chunk(name='abcdefg', data=data)
        f.end_frame()
        f.write_chunk(name='test', data=data)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test_chunk_exists.gsd',
        mode=open_mode.read,
        application='test_chunk_exists',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.chunk_exists(frame=0, name='chunk1')
        read_data = f.read_chunk(frame=0, name='chunk1')
        assert f.chunk_exists(frame=1, name='abcdefg')
        read_data = f.read_chunk(frame=1, name='abcdefg')
        assert f.chunk_exists(frame=2, name='test')
        read_data = f.read_chunk(frame=2, name='test')

        assert not f.chunk_exists(frame=1, name='chunk1')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=1, name='chunk1')
        assert not f.chunk_exists(frame=2, name='abcdefg')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=2, name='abcdefg')
        assert not f.chunk_exists(frame=0, name='test')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=0, name='test')

        assert not f.chunk_exists(frame=2, name='chunk1')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=2, name='chunk1')
        assert not f.chunk_exists(frame=0, name='abcdefg')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=0, name='abcdefg')
        assert not f.chunk_exists(frame=1, name='test')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=1, name='test')

    # test again with pygsd
    with gsd.pygsd.GSDFile(
        file=open(str(tmp_path / 'test_chunk_exists.gsd'), mode='rb')
    ) as f:
        assert f.chunk_exists(frame=0, name='chunk1')
        read_data = f.read_chunk(frame=0, name='chunk1')
        assert f.chunk_exists(frame=1, name='abcdefg')
        read_data = f.read_chunk(frame=1, name='abcdefg')
        assert f.chunk_exists(frame=2, name='test')
        read_data = f.read_chunk(frame=2, name='test')

        assert not f.chunk_exists(frame=1, name='chunk1')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=1, name='chunk1')
        assert not f.chunk_exists(frame=2, name='abcdefg')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=2, name='abcdefg')
        assert not f.chunk_exists(frame=0, name='test')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=0, name='test')

        assert not f.chunk_exists(frame=2, name='chunk1')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=2, name='chunk1')
        assert not f.chunk_exists(frame=0, name='abcdefg')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=0, name='abcdefg')
        assert not f.chunk_exists(frame=1, name='test')
        with pytest.raises(KeyError):
            read_data = f.read_chunk(frame=1, name='test')  # noqa


def test_readonly_errors(tmp_path, open_mode):
    """Test that read only files provide the appropriate errors."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)
    with gsd.fl.open(
        name=tmp_path / 'test_readonly_errors.gsd',
        mode=open_mode.write,
        application='test_readonly_errors',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        for _i in range(10):
            f.write_chunk(name='chunk1', data=data)
            f.end_frame()

    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)
    with gsd.fl.open(
        name=tmp_path / 'test_readonly_errors.gsd',
        mode='r',
        application='test_readonly_errors',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        with pytest.raises(RuntimeError):
            f.end_frame()

        with pytest.raises(RuntimeError):
            f.write_chunk(name='chunk1', data=data)

    # test again with pygsd
    with gsd.pygsd.GSDFile(
        file=open(str(tmp_path / 'test_readonly_errors.gsd'), mode='rb')
    ) as f:
        with pytest.raises(NotImplementedError):
            f.end_frame()

        with pytest.raises(NotImplementedError):
            f.write_chunk(name='chunk1', data=data)


def test_fileio_errors(tmp_path, open_mode):
    """Test that OS file I/O errors pass through."""
    # These test cause python to crash on windows....
    if platform.system() != 'Windows':
        with pytest.raises(FileNotFoundError):
            gsd.fl.open(
                name='/this/file/does/not/exist',
                mode='r',
                application='test_readonly_errors',
                schema='none',
                schema_version=[1, 2],
            )

        with open(str(tmp_path / 'test_fileio_errors.gsd'), 'wb') as f:
            f.write(b'test')

        with pytest.raises(RuntimeError):
            f = gsd.fl.open(
                name=tmp_path / 'test_fileio_errors.gsd',
                mode=open_mode.read,
                application='test_readonly_errors',
                schema='none',
                schema_version=[1, 2],
            )


def test_dtype_errors(tmp_path, open_mode):
    """Test that unsupported data types result in errors."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.bool_)

    with gsd.fl.open(
        name=tmp_path / 'test_dtype_errors1.gsd',
        mode=open_mode.write,
        application='test_dtype_errors',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        with pytest.raises(ValueError):
            f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.float16)

    with gsd.fl.open(
        name=tmp_path / 'test_dtype_errors2.gsd',
        mode=open_mode.write,
        application='test_dtype_errors',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        with pytest.raises(ValueError):
            f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.complex64)

    with gsd.fl.open(
        name=tmp_path / 'test_dtype_errors3.gsd',
        mode=open_mode.write,
        application='test_dtype_errors',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        with pytest.raises(ValueError):
            f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.complex128)

    with gsd.fl.open(
        name=tmp_path / 'test_dtype_errors4.gsd',
        mode=open_mode.write,
        application='test_dtype_errors',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        with pytest.raises(ValueError):
            f.write_chunk(name='chunk1', data=data)

        f.end_frame()


def test_truncate(tmp_path):
    """Test that the truncate method functions."""
    rng = numpy.random.default_rng()
    data = numpy.ascontiguousarray(rng.random(size=(1000, 3)), dtype=numpy.float32)
    with gsd.fl.open(
        name=tmp_path / 'test_truncate.gsd',
        mode='w',
        application='test_truncate',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.mode == 'w'
        for _i in range(10):
            f.write_chunk(name='data', data=data)
            f.end_frame()

        f.flush()

        assert f.nframes == 10

        f.truncate()
        assert f.nframes == 0
        assert f.application == 'test_truncate'
        assert f.schema == 'none'
        assert f.schema_version == (1, 2)

        f.write_chunk(name='data', data=data)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test_truncate.gsd',
        mode='r',
        application='test_truncate',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.name == str(tmp_path / 'test_truncate.gsd')
        assert f.mode == 'r'
        assert f.application == 'test_truncate'
        assert f.schema == 'none'
        assert f.schema_version == (1, 2)
        assert f.nframes == 1


def test_namelen(tmp_path, open_mode):
    """Test that long names are truncated as documented."""
    app_long = 'abcdefga' * 100
    schema_long = 'ijklmnop' * 100
    chunk_long = '12345678' * 100

    with gsd.fl.open(
        name=tmp_path / 'test_namelen.gsd',
        mode=open_mode.write,
        application=app_long,
        schema=schema_long,
        schema_version=[1, 2],
    ) as f:
        assert f.application == app_long[0:63]
        assert f.schema == schema_long[0:63]

        data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)
        f.write_chunk(name=chunk_long, data=data)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test_namelen.gsd',
        mode=open_mode.read,
        application=app_long,
        schema=schema_long,
        schema_version=[1, 2],
    ) as f:
        data_read = f.read_chunk(0, name=chunk_long)
        numpy.testing.assert_array_equal(data, data_read)

    # test again with pygsd
    with gsd.pygsd.GSDFile(
        file=open(str(tmp_path / 'test_namelen.gsd'), mode='rb')
    ) as f:
        data_read = f.read_chunk(0, name=chunk_long)
        numpy.testing.assert_array_equal(data, data_read)


def test_open(tmp_path):
    """Test the open() API."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode='x',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test_2.gsd',
        mode='x',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()
        f.flush()
        f.read_chunk(0, name='chunk1')

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode='w',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode='w',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()
        f.flush()
        f.read_chunk(0, name='chunk1')

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode='a',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode='r',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.read_chunk(0, name='chunk1')
        f.read_chunk(1, name='chunk1')

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode='r+',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()
        f.read_chunk(0, name='chunk1')
        f.read_chunk(1, name='chunk1')
        f.flush()
        f.read_chunk(2, name='chunk1')


def test_find_matching_chunk_names(tmp_path, open_mode):
    """Test the find_matching_chunk_names API."""
    data = numpy.array([1, 2, 3, 4, 5], dtype=numpy.float32)

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode=open_mode.write,
        application='test_find_matching_chunk_names',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='log/A', data=data)
        f.write_chunk(name='log/chunk2', data=data)
        f.end_frame()
        f.write_chunk(name='data/B', data=data)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode=open_mode.read,
        application='test_find_matching_chunk_names',
        schema='none',
        schema_version=[1, 2],
    ) as f:
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
    with gsd.pygsd.GSDFile(file=open(str(tmp_path / 'test.gsd'), mode='rb')) as f:
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
    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode=open_mode.write,
        application='test_chunk_name_limit',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        for i in range(65535):
            f.write_chunk(name=str(i), data=numpy.array([i], dtype=numpy.int32))

        # The GSD specification limits to 65535 names:
        with pytest.raises(RuntimeError):
            f.write_chunk(name='65536', data=numpy.array([i], dtype=numpy.int32))


def test_many_names(tmp_path, open_mode):
    """Test that many chunk names can be written to a file."""
    values = list(range(1000))

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode=open_mode.write,
        application='test_many_names',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        for _ in range(5):
            random.shuffle(values)
            for value in values:
                f.write_chunk(
                    name=str(value), data=numpy.array([value * 13], dtype=numpy.int32)
                )
            f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test.gsd',
        mode=open_mode.read,
        application='test_many_names',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        for frame in range(5):
            random.shuffle(values)
            for value in values:
                data = numpy.array([value * 13], dtype=numpy.int32)
                data_read = f.read_chunk(frame=frame, name=str(value))
                numpy.testing.assert_array_equal(data, data_read)

    with gsd.pygsd.GSDFile(file=open(str(tmp_path / 'test.gsd'), mode='rb')) as f:
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

    # test with the C implementation
    with gsd.fl.open(
        name=test_path / 'test_gsd_v1.gsd',
        mode='r',
        application='test_gsd_v1',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        check_v1_file_read(f)

    # and the pure Python implementation
    with gsd.pygsd.GSDFile(
        file=open(str(test_path / 'test_gsd_v1.gsd'), mode='rb')
    ) as f:
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

    with gsd.fl.open(
        name=tmp_path / 'test_gsd_v1.gsd',
        mode='r+',
        application='test_gsd_v1',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.gsd_version == (1, 0)

        f.upgrade()

        # check that we can read the file contents after the upgrade in memory
        check_v1_file_read(f)

    # and the same tests again after closing and opening the file
    with gsd.fl.open(
        name=tmp_path / 'test_gsd_v1.gsd',
        mode=open_mode.read,
        application='test_gsd_v1',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.gsd_version == current_gsd_version

        check_v1_file_read(f)

    with gsd.pygsd.GSDFile(
        file=open(str(tmp_path / 'test_gsd_v1.gsd'), mode='rb')
    ) as f:
        assert f.gsd_version == current_gsd_version

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
        check_v = v
        if isinstance(v, str) and len(v) > 63:
            # v1 files truncate names to 63 chars
            check_v = v[0:63]
        values_str.append(str(check_v))
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
            check_value = value

            if isinstance(value, int):
                data = numpy.array([value * 13], dtype=numpy.int32)
            else:
                data = numpy.array([hash(value)], dtype=numpy.int64)
                # v1 files truncate names to 63 chars
                if len(value) > 63:
                    check_value = value[0:63]

            data_read = f.read_chunk(frame=frame, name=str(check_value))
            numpy.testing.assert_array_equal(data, data_read)

    # test that we can write new entries to the file
    with gsd.fl.open(
        name=tmp_path / 'test_gsd_v1.gsd',
        mode='r+',
        application='test_gsd_v1',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.gsd_version == (1, 0)

        for value in values:
            if isinstance(value, int):
                data = numpy.array([value * 13], dtype=numpy.int32)
            else:
                data = numpy.array([hash(value)], dtype=numpy.int64)
            f.write_chunk(name=str(value), data=data)
        f.end_frame()
        f.flush()

        check_v1_file_read(f)

    # test opening again with the C implementation
    with gsd.fl.open(
        name=tmp_path / 'test_gsd_v1.gsd',
        mode=open_mode.read,
        application='test_gsd_v1',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        check_v1_file_read(f)

    # and the pure Python implementation
    with gsd.pygsd.GSDFile(
        file=open(str(tmp_path / 'test_gsd_v1.gsd'), mode='rb')
    ) as f:
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
            if isinstance(value, int):
                data = numpy.array([value * 13], dtype=numpy.int32)
            else:
                data = numpy.array([hash(value)], dtype=numpy.int64)

            data_read = f.read_chunk(frame=frame, name=str(value))
            numpy.testing.assert_array_equal(data, data_read)

    # test that we can write new entries to the file
    with gsd.fl.open(
        name=tmp_path / 'test_gsd_v1.gsd',
        mode='r+',
        application='test_gsd_v1',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.gsd_version == (1, 0)

        f.upgrade()

        assert f.gsd_version == current_gsd_version

        for value in values:
            if isinstance(value, int):
                data = numpy.array([value * 13], dtype=numpy.int32)
            else:
                data = numpy.array([hash(value)], dtype=numpy.int64)
            f.write_chunk(name=str(value), data=data)
        f.end_frame()
        f.flush()

        check_v1_file_read(f)

    # test opening again with the C implementation
    with gsd.fl.open(
        name=tmp_path / 'test_gsd_v1.gsd',
        mode=open_mode.read,
        application='test_gsd_v1',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.gsd_version == current_gsd_version

        check_v1_file_read(f)

    # and the pure Python implementation
    with gsd.pygsd.GSDFile(
        file=open(str(tmp_path / 'test_gsd_v1.gsd'), mode='rb')
    ) as f:
        assert f.gsd_version == current_gsd_version

        check_v1_file_read(f)


def test_zero_size(tmp_path, open_mode):
    """Test that zero-size data chunks are allowed."""
    data = numpy.array([], dtype=numpy.float32)

    with gsd.fl.open(
        name=tmp_path / 'test_zero.gsd',
        mode=open_mode.write,
        application='test_zero',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='data', data=data)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / 'test_zero.gsd',
        mode=open_mode.read,
        application='test_zero',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.nframes == 1
        data_read = f.read_chunk(frame=0, name='data')
        assert data_read.shape == (0,)
        assert data_read.dtype == numpy.float32

    # test again with pygsd
    with gsd.pygsd.GSDFile(file=open(str(tmp_path / 'test_zero.gsd'), mode='rb')) as f:
        assert f.nframes == 1
        data_read = f.read_chunk(frame=0, name='data')
        assert data_read.shape == (0,)
        assert data_read.dtype == numpy.float32


@pytest.mark.skipif(
    sys.version_info < (3, 7), reason='Python 3.6 fails to handle non-ascii characters.'
)
def test_utf8(tmp_path):
    """Test that the API handles UTF-8 encoding for the filename."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)

    fname = '中文.gsd'

    with gsd.fl.open(
        name=tmp_path / fname,
        mode='x',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    dir_list = os.listdir(tmp_path)
    print(dir_list)
    assert fname in dir_list

    with gsd.fl.open(
        name=tmp_path / fname,
        mode='w',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()

    with gsd.fl.open(
        name=tmp_path / fname,
        mode='r',
        application='test_open',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.read_chunk(0, name='chunk1')


@pytest.mark.parametrize('mode', ['w', 'x', 'a', 'r+'])
def test_read_write(tmp_path, mode):
    """Test that data chunks can read from files opened in all write modes."""
    if mode[0] == 'r' or mode[0] == 'a':
        with gsd.fl.open(
            name=tmp_path / 'test_read_write.gsd',
            mode='w',
            application='test_read_write',
            schema='none',
            schema_version=[1, 2],
        ):
            pass

    data = numpy.array([10], dtype=numpy.int64)
    nframes = 1024

    with gsd.fl.open(
        name=tmp_path / 'test_read_write.gsd',
        mode=mode,
        application='test_read_write',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.mode == mode
        for i in range(nframes):
            data[0] = i
            f.write_chunk(name='data1', data=data)
            data[0] = i * 10
            f.write_chunk(name='data10', data=data)
            f.end_frame()

        assert f.nframes == 0
        f.flush()
        assert f.nframes == nframes

        for i in range(nframes):
            data1 = f.read_chunk(frame=i, name='data1')
            data10 = f.read_chunk(frame=i, name='data10')
            assert data1[0] == i
            assert data10[0] == i * 10


@pytest.mark.parametrize('mode', ['w', 'x', 'a', 'r+'])
def test_read_after_write(tmp_path, mode):
    """Test that chunk names and data are not present until after flush."""
    if mode[0] == 'r' or mode[0] == 'a':
        with gsd.fl.open(
            name=tmp_path / 'test_read_write.gsd',
            mode='w',
            application='test_read_write',
            schema='none',
            schema_version=[1, 2],
        ):
            pass

    data = numpy.array([10], dtype=numpy.int64)

    with gsd.fl.open(
        name=tmp_path / 'test_read_write.gsd',
        mode=mode,
        application='test_read_write',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        data[0] = 1
        f.write_chunk(name='data1', data=data)
        data[0] = 2
        f.write_chunk(name='data2', data=data)
        f.end_frame()

        assert f.nframes == 0
        assert not f.chunk_exists(frame=0, name='data1')
        assert not f.chunk_exists(frame=0, name='data2')
        assert len(f.find_matching_chunk_names('')) == 0
        with pytest.raises(KeyError):
            f.read_chunk(frame=0, name='data1')
        with pytest.raises(KeyError):
            f.read_chunk(frame=0, name='data2')

        f.flush()

        data[0] = 1
        f.write_chunk(name='data1', data=data)
        data[0] = 2
        f.write_chunk(name='data2', data=data)
        data[0] = 3
        f.write_chunk(name='data3', data=data)
        f.end_frame()

        # frame 0 should now be available, but not frame 1
        assert f.nframes == 1
        assert f.chunk_exists(frame=0, name='data1')
        assert f.chunk_exists(frame=0, name='data2')
        assert len(f.find_matching_chunk_names('')) == 2
        data1 = f.read_chunk(frame=0, name='data1')
        assert data1[0] == 1
        data2 = f.read_chunk(frame=0, name='data2')
        assert data2[0] == 2

        assert not f.chunk_exists(frame=1, name='data1')
        assert not f.chunk_exists(frame=1, name='data2')
        assert not f.chunk_exists(frame=1, name='data3')
        with pytest.raises(KeyError):
            f.read_chunk(frame=1, name='data1')
        with pytest.raises(KeyError):
            f.read_chunk(frame=1, name='data2')
        with pytest.raises(KeyError):
            f.read_chunk(frame=1, name='data3')

        f.flush()
        assert f.nframes == 2
        assert f.chunk_exists(frame=1, name='data1')
        assert f.chunk_exists(frame=1, name='data2')
        assert f.chunk_exists(frame=1, name='data3')
        assert len(f.find_matching_chunk_names('')) == 3
        data1 = f.read_chunk(frame=1, name='data1')
        assert data1[0] == 1
        data2 = f.read_chunk(frame=1, name='data2')
        assert data2[0] == 2
        data3 = f.read_chunk(frame=1, name='data3')
        assert data3[0] == 3


@pytest.mark.parametrize('n_flush', [0, 1, 2])
def test_flush(tmp_path, open_mode, n_flush):
    """Test flush."""
    data = numpy.array([1, 2, 3, 4, 5, 10012], dtype=numpy.int64)
    with gsd.fl.open(
        name=tmp_path / 'test_flush.gsd',
        mode=open_mode.write,
        application='test_flush',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        f.write_chunk(name='chunk1', data=data)
        f.end_frame()
        f.write_chunk(name='chunk2', data=data)
        f.end_frame()
        f.write_chunk(name='chunk3', data=data)

        # Ensure that the data is buffered by opening the file with a 2nd
        # handle read-only and checking it.
        with gsd.fl.open(name=tmp_path / 'test_flush.gsd', mode='r') as f_readonly:
            assert not f_readonly.chunk_exists(frame=0, name='chunk1')
            assert not f_readonly.chunk_exists(frame=1, name='chunk2')
            assert f_readonly.nframes == 0

        # 0 calls to flush tests the implicit flush on close, 2 calls to flush
        # tests that repeated calls are handled properly.
        for _i in range(n_flush):
            f.flush()

    with gsd.fl.open(name=tmp_path / 'test_flush.gsd', mode=open_mode.read) as f:
        assert f.chunk_exists(frame=0, name='chunk1')
        assert f.chunk_exists(frame=1, name='chunk2')

        # The third chunk is not written because end_frame is not called a
        # third time.

        assert not f.chunk_exists(frame=2, name='chunk3')
        assert f.nframes == 2


def test_maximum_write_buffer_size(tmp_path, open_mode):
    """Test maximum_write_buffer_size."""
    file_name = tmp_path / 'test_maximum_write_buffer_size.gsd'

    data_256 = numpy.full(shape=(256,), fill_value=1, dtype=numpy.uint8)
    data_1024 = numpy.full(shape=(1024,), fill_value=2, dtype=numpy.uint8)

    with gsd.fl.open(
        name=file_name,
        mode=open_mode.write,
        application='test_maximum_write_buffer_size',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        assert f.maximum_write_buffer_size > 0

        with pytest.raises(RuntimeError):
            f.maximum_write_buffer_size = 0

        f.maximum_write_buffer_size = 1024
        assert f.maximum_write_buffer_size == 1024

        initial_size = os.path.getsize(file_name)
        f.write_chunk(name='data', data=data_256)
        assert os.path.getsize(file_name) == initial_size
        f.end_frame()

        f.write_chunk(name='data', data=data_256)
        f.end_frame()
        f.write_chunk(name='data', data=data_256)
        f.end_frame()
        f.write_chunk(name='data', data=data_256)
        f.end_frame()
        f.write_chunk(name='data', data=data_256)
        assert os.path.getsize(file_name) == initial_size + 1024
        f.end_frame()
        f.write_chunk(name='data', data=data_1024)
        assert os.path.getsize(file_name) == initial_size + 1024 + 256 + 1024
        f.end_frame()

        f.flush()
        for i in range(5):
            numpy.testing.assert_array_equal(f.read_chunk(i, 'data'), data_256)
        numpy.testing.assert_array_equal(f.read_chunk(5, 'data'), data_1024)


def test_file_exists_error():
    """Test that IO errors throw the correct Python Exception."""
    with pytest.raises(FileExistsError):
        with gsd.fl.open(
            name=test_path / 'test_gsd_v1.gsd',
            mode='x',
            application='test_gsd_v1',
            schema='none',
            schema_version=[1, 2],
        ):
            pass


def test_pending_index_entries(tmp_path):
    """Ensure that gsd preserves pending index entries."""
    with gsd.fl.open(
        tmp_path / 'test_pending_index_entries.gsd',
        'w',
        application='My application',
        schema='My Schema',
        schema_version=[1, 0],
    ) as f:
        # Frame 0 must be complete to trigger the bug.
        f.write_chunk(name='0', data=numpy.array([0]))
        f.end_frame()

        for i in range(16):
            f.write_chunk(name=str(i), data=numpy.array([i]))

        # Flush with pending chunks in the frame index.
        f.flush()

        f.end_frame()
        f.flush()

        # All test chunks should be present in the file.
        for i in range(16):
            assert f.chunk_exists(name=str(i), frame=1)


@pytest.mark.parametrize('flush', [False, True])
def test_expand(tmp_path, open_mode, flush):
    """Test index expansion."""
    with gsd.fl.open(
        name=tmp_path / 'test_expand.gsd',
        mode=open_mode.write,
        application='test_expand',
        schema='none',
        schema_version=[1, 2],
    ) as f:
        N_ENTRIES = 1024

        for i in range(N_ENTRIES):
            data = numpy.array([i] * i, dtype=numpy.int64)
            f.write_chunk(name='data', data=data)
            f.end_frame()
            if flush and (i & 0xF == 0):
                f.flush()

        f.flush()
        for i in range(N_ENTRIES):
            expected_data = numpy.array([i] * i, dtype=numpy.int64)
            read_data = f.read_chunk(frame=i, name='data')
            numpy.testing.assert_array_equal(read_data, expected_data)
