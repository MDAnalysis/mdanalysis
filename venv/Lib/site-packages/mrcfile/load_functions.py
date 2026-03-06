# Copyright (c) 2018, Science and Technology Facilities Council
# This software is distributed under a BSD licence. See LICENSE.txt.

"""
load_functions
--------------

Module for top-level functions that open MRC files and form the main API of
the package.

"""

# Import Python 3 features for future-proofing
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import io
import os

import numpy as np

from . import utils
from .bzip2mrcfile import Bzip2MrcFile
from .constants import MAP_ID, MAP_ID_OFFSET_BYTES
from .future_mrcfile import FutureMrcFile
from .gzipmrcfile import GzipMrcFile
from .mrcfile import MrcFile
from .mrcmemmap import MrcMemmap


def new(name, data=None, compression=None, overwrite=False):
    """Create a new MRC file.
    
    Args:
        name: The file name to use, as a string or :class:`~pathlib.Path`.
        data: Data to put in the file, as a :class:`numpy array
            <numpy.ndarray>`. The default is :data:`None`, to create an empty
            file.
        compression: The compression format to use. Acceptable values are:
            :data:`None` (the default; for no compression), ``'gzip'`` or
            ``'bzip2'``.
            It's good practice to name compressed files with an appropriate
            extension (for example, ``.mrc.gz`` for gzip) but this is not
            enforced.
        overwrite: Flag to force overwriting of an existing file. If
            :data:`False` and a file of the same name already exists, the file
            is not overwritten and an exception is raised.
    
    Returns:
        An :class:`~mrcfile.mrcfile.MrcFile` object (or a
        subclass of it if ``compression`` is specified).
    
    Raises:
        :exc:`ValueError`: If the file already exists and overwrite is
            :data:`False`.
        :exc:`ValueError`: If the compression format is not recognised.

    Warns:
        RuntimeWarning: If the data array contains Inf or NaN values.
    """
    if compression == 'gzip':
        NewMrc = GzipMrcFile
    elif compression == 'bzip2':
        NewMrc = Bzip2MrcFile
    elif compression is not None:
        raise ValueError("Unknown compression format '{0}'"
                         .format(compression))
    else:
        NewMrc = MrcFile
    mrc = NewMrc(name, mode='w+', overwrite=overwrite)
    if data is not None:
        mrc.set_data(data)
    return mrc


def open(name, mode='r', permissive=False, header_only=False):  # @ReservedAssignment
    """Open an MRC file.
    
    This function opens both normal and compressed MRC files. Supported
    compression formats are: gzip, bzip2.
    
    It is possible to use this function to create new MRC files (using mode
    ``w+``) but the :func:`new` function is more flexible.
    
    This function offers a permissive read mode for attempting to open corrupt
    or invalid files. In permissive mode, :mod:`warnings` are issued instead of
    exceptions if problems with the file are encountered. See
    :class:`~mrcfile.mrcinterpreter.MrcInterpreter` or the
    :doc:`usage guide <../usage_guide>` for more information.
    
    Args:
        name: The file name to open, as a string or :class:`~pathlib.Path`.
        mode: The file mode to use. This should be one of the following: ``r``
            for read-only, ``r+`` for read and write, or ``w+`` for a new empty
            file. The default is ``r``.
        permissive: Read the file in permissive mode. The default is
            :data:`False`.
        header_only: Only read the header (and extended header) from the file.
            The default is :data:`False`.
    
    Returns:
        An :class:`~mrcfile.mrcfile.MrcFile` object (or a
        :class:`~mrcfile.gzipmrcfile.GzipMrcFile` object if the file is
        gzipped).
    
    Raises:
        :exc:`ValueError`: If the mode is not one of ``r``, ``r+`` or ``w+``.
        :exc:`ValueError`: If the file is not a valid MRC file and
            ``permissive`` is :data:`False`.
        :exc:`ValueError`: If the mode is ``w+`` and the file already exists.
            (Call :func:`new` with ``overwrite=True`` to deliberately overwrite
            an existing file.)
        :exc:`OSError`: If the mode is ``r`` or ``r+`` and the file does not
            exist.
    
    Warns:
        RuntimeWarning: If the file appears to be a valid MRC file but the data
            block is longer than expected from the dimensions in the header.
        RuntimeWarning: If the file is not a valid MRC file and ``permissive``
            is :data:`True`.
        RuntimeWarning: If the header's ``exttyp`` field is set to a known
            value but the extended header's size is not a multiple of the
            number of bytes in the corresponding dtype.
    """
    NewMrc = MrcFile
    name = str(name)  # in case name is a pathlib Path
    if os.path.exists(name):
        if "w" in mode:
            raise ValueError("File '{0}' already exists; use a different name, or"
                             " delete it first or call 'mrcfile.new()' with"
                             " 'overwrite=True' to overwrite it".format(name))
        with io.open(name, 'rb') as f:
            start = f.read(MAP_ID_OFFSET_BYTES + len(MAP_ID))
        # Check for map ID string to avoid trying to decompress normal files
        # where the nx value happens to include the magic number for a
        # compressed format. (This still risks failing to correctly decompress
        # compressed files which happen to have 'MAP ' at position 208, but
        # that is less likely and if it does occur, the CompressedMrcFile
        # class can always be used directly instead.)
        if start[-len(MAP_ID):] != MAP_ID:
            if start[:2] == b'\x1f\x8b':
                NewMrc = GzipMrcFile
            elif start[:2] == b'BZ':
                NewMrc = Bzip2MrcFile
    return NewMrc(name, mode=mode, permissive=permissive,
                  header_only=header_only)


def read(name):
    """Read an MRC file's data into a numpy array.

    This is a convenience function to read the data from an MRC file when there is no
    need for the file's header information. To read the headers as well, or if you need
    access to an :class:`~mrcfile.mrcfile.MrcFile` object representing the file, use
    :func:`mrcfile.open` instead.

    Args:
        name: The file name to read, as a string or :class:`~pathlib.Path`.

    Returns:
        A :class:`numpy array<numpy.ndarray>` containing the data from the file.
    """
    with open(name, mode='r', permissive=True) as mrc:
        data = mrc.data.copy()
    return data


def write(name, data=None, overwrite=False, voxel_size=None):
    """Write a new MRC file.

    This is a convenience function to allow data to be quickly written to a file (with
    optional compression) using just a single function call. However, there is no
    control over the file's metadata except for optionally setting the voxel size. For
    more control, or if you need access to an :class:`~mrcfile.mrcfile.MrcFile` object
    representing the new file, use :func:`mrcfile.new` instead.

    Args:
        name: The file name to use, as a string or :class:`~pathlib.Path`. If the name
            ends with ``.gz`` or ``.bz2``, the file will be compressed using gzip or
            bzip2 respectively.
        data: Data to put in the file, as a :class:`numpy array
            <numpy.ndarray>`. The default is :data:`None`, to create an empty
            file.
        overwrite: Flag to force overwriting of an existing file. If
            :data:`False` and a file of the same name already exists, the file
            is not overwritten and an exception is raised.
        voxel_size: float | 3-tuple
            The voxel size to be written in the file header.

    Raises:
        :exc:`ValueError`: If the file already exists and overwrite is
            :data:`False`.

    Warns:
        RuntimeWarning: If the data array contains Inf or NaN values.
    """
    name = str(name)  # in case name is a pathlib Path
    compression = None
    if name.endswith('.gz'):
        compression = 'gzip'
    elif name.endswith('.bz2'):
        compression = 'bzip2'
    with new(name, data, compression, overwrite) as mrc:
        if voxel_size is not None:
            mrc.voxel_size = voxel_size


def open_async(name, mode='r', permissive=False):
    """Open an MRC file asynchronously in a separate thread.

    This allows a file to be opened in the background while the main thread
    continues with other work. This can be a good way to improve performance if
    the main thread is busy with intensive computation, but will be less
    effective if the main thread is itself busy with disk I/O.

    Multiple files can be opened in the background simultaneously. However,
    this implementation is relatively crude; each call to this function will
    start a new thread and immediately use it to start opening a file. If you
    try to open many large files at the same time, performance will decrease as
    all of the threads attempt to access the disk at once. You'll also risk
    running out of memory to store the data from all the files.

    This function returns a :class:`~mrcfile.future_mrcfile.FutureMrcFile`
    object, which deliberately mimics the API of the
    :class:`~concurrent.futures.Future` object from Python 3's
    :mod:`concurrent.futures` module. (Future versions of this library might
    return genuine :class:`~concurrent.futures.Future` objects instead.)

    To get the real :class:`~mrcfile.mrcfile.MrcFile` object from a
    :class:`~mrcfile.future_mrcfile.FutureMrcFile`, call
    :meth:`~mrcfile.future_mrcfile.FutureMrcFile.result`. This will block until
    the file has been read and the :class:`~mrcfile.mrcfile.MrcFile` object is
    ready. To check if the :class:`~mrcfile.mrcfile.MrcFile` is ready without
    blocking, call :meth:`~mrcfile.future_mrcfile.FutureMrcFile.running` or
    :meth:`~mrcfile.future_mrcfile.FutureMrcFile.done`.

    Args:
        name: The file name to open, as a string or :class:`~pathlib.Path`.
        mode: The file mode (one of ``r``, ``r+`` or ``w+``).
        permissive: Read the file in permissive mode. The default is
            :data:`False`.

    Returns:
        A :class:`~mrcfile.future_mrcfile.FutureMrcFile` object.

    """
    return FutureMrcFile(open, (name,), dict(mode=mode, permissive=permissive))


def mmap(name, mode='r', permissive=False):
    """Open a memory-mapped MRC file.
    
    This allows much faster opening of large files, because the data is only
    accessed on disk when a slice is read or written from the data array. See
    the :class:`~mrcfile.mrcmemmap.MrcMemmap` class documentation for more
    information.
    
    Because the memory-mapped data array accesses the disk directly, compressed
    files cannot be opened with this function. In all other ways, :func:`mmap`
    behaves in exactly the same way as :func:`open`. The
    :class:`~mrcfile.mrcmemmap.MrcMemmap` object returned by this function can
    be used in exactly the same way as a normal
    :class:`~mrcfile.mrcfile.MrcFile` object.
    
    Args:
        name: The file name to open, as a string or :class:`~pathlib.Path`.
        mode: The file mode (one of ``r``, ``r+`` or ``w+``).
        permissive: Read the file in permissive mode. The default is
            :data:`False`.
    
    Returns:
        An :class:`~mrcfile.mrcmemmap.MrcMemmap` object.
    """
    return MrcMemmap(name, mode=mode, permissive=permissive)


def new_mmap(name, shape, mrc_mode=0, fill=None, overwrite=False, extended_header=None, exttyp=None):
    """Create a new, empty memory-mapped MRC file.

    This function is useful for creating very large files. The initial contents
    of the data array can be set with the ``fill`` parameter if needed, but be
    aware that filling a large array can take a long time.

    If ``fill`` is not set, the new data array's contents are unspecified and
    system-dependent. (Some systems fill a new empty mmap with zeros, others
    fill it with the bytes from the disk at the newly-mapped location.) If you
    are definitely going to fill the entire array with new data anyway you can
    safely leave ``fill`` as :data:`None`, otherwise it is advised to use a
    sensible fill value (or ensure you are on a system that fills new mmaps
    with a reasonable default value).

    Args:
        name: The file name to use, as a string or :class:`~pathlib.Path`.
        shape: The shape of the data array to open, as a 2-, 3- or 4-tuple of
            ints. For example, ``(nz, ny, nx)`` for a new 3D volume, or
            ``(ny, nx)`` for a new 2D image.
        mrc_mode: The MRC mode to use for the new file. One of 0, 1, 2, 4 or 6,
            which correspond to numpy dtypes as follows:

            * mode 0 -> int8
            * mode 1 -> int16
            * mode 2 -> float32
            * mode 4 -> complex64
            * mode 6 -> uint16

            The default is 0.
        fill: An optional value to use to fill the new data array. If
            :data:`None`, the data array will not be filled and its contents
            are unspecified. Numpy's usual rules for rounding or rejecting
            values apply, according to the dtype of the array.
        overwrite: Flag to force overwriting of an existing file. If
            :data:`False` and a file of the same name already exists, the file
            is not overwritten and an exception is raised.
        extended_header: The extended header object
        exttyp: The extended header type

    Returns:
        A new :class:`~mrcfile.mrcmemmap.MrcMemmap` object.

    Raises:
        :exc:`ValueError`: If the MRC mode is invalid.
        :exc:`ValueError`: If the file already exists and overwrite is
            :data:`False`.
        :exc:`ValueError`: If the extended header or any of the dimensions of
                the shape have a size greater than 2,147,483,647 (and
                therefore the size cannot be stored in the header).
    """
    for dim in shape:
        if dim > np.iinfo(np.int32).max:
            raise ValueError("New shape is too large! Found a dimension of"
                             " size {}. The maximum allowed is {}."
                             .format(dim, np.iinfo(np.int32).max))

    mrc = MrcMemmap(name, mode='w+', overwrite=overwrite)

    # Add the extended header and type. We need to do this before creating the
    # memory mapped file to avoid having to copy lots of data
    if extended_header is not None:
        mrc.set_extended_header(extended_header)
    if exttyp is not None:
        mrc.header.exttyp = exttyp

    dtype = utils.dtype_from_mode(mrc_mode)
    mrc._open_memmap(dtype, shape)
    mrc.update_header_from_data()
    if fill is not None:
        mrc.data[...] = fill
    mrc.flush()
    return mrc
