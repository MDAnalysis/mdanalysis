# Copyright (c) 2016-2023 The Regents of the University of Michigan
# Part of GSD, released under the BSD 2-Clause License.

"""Test the gsd.fl API with large files."""

import gsd.hoomd
import numpy
import pytest
import gc


@pytest.mark.validate
@pytest.mark.parametrize("N", [2**27, 2**28, 2**29 + 1])
def test_large_n(tmp_path, N):
    """Test data chunks and files larger than 2 GB."""
    gc.collect()

    data = numpy.linspace(0, N, num=N, endpoint=False, dtype=numpy.uint32)
    with gsd.fl.open(name=tmp_path / 'test_large_N.gsd',
                     mode='xb',
                     application='test_large_N',
                     schema='none',
                     schema_version=[1, 2]) as f:
        f.write_chunk(name='data', data=data)
        f.end_frame()

    with gsd.fl.open(name=tmp_path / 'test_large_N.gsd',
                     mode='rb',
                     application='test_large_N',
                     schema='none',
                     schema_version=[1, 2]) as f:
        read_data = f.read_chunk(frame=0, name='data')

        # compare the array with memory usage so this test can pass on CI
        # platforms
        diff = (data - read_data)
        data = None
        read_data = None
        gc.collect()
        diff = diff**2

        assert numpy.sum(diff) == 0
