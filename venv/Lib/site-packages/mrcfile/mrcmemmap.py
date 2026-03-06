# Copyright (c) 2016, Science and Technology Facilities Council
# This software is distributed under a BSD licence. See LICENSE.txt.
"""
mrcmemmap
---------

Module which exports the :class:`MrcMemmap` class.

Classes:
    :class:`MrcMemmap`: An MrcFile subclass that uses a memory-mapped data
    array.

"""

# Import Python 3 features for future-proofing
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import mmap
import os
import warnings

import numpy as np

import mrcfile.utils as utils
from .mrcfile import MrcFile


class MrcMemmap(MrcFile):
    
    """MrcFile subclass that uses a :class:`numpy memmap array <numpy.memmap>`
    for the data.
    
    Using a memmap means that the disk access is done lazily: the data array
    will only be read or written in small chunks when required. To access the
    contents of the array, use the array slice operator.
    
    Usage is the same as for :class:`~mrcfile.mrcfile.MrcFile`.
    
    Note that memmap arrays use a fairly small chunk size and so performance
    could be poor on file systems that are optimised for infrequent large I/O
    operations.
    
    If required, it is possible to create a very large empty file by creating a
    new MrcMemmap and then calling :meth:`_open_memmap` to create the memmap
    array, which can then be filled slice-by-slice. Be aware that the contents
    of a new, empty memmap array depend on your platform: the data values
    could be garbage or zeros.
    
    """
    
    def __repr__(self):
        return "MrcMemmap('{0}', mode='{1}')".format(self._iostream.name,
                                                     self._mode)
    
    def set_extended_header(self, extended_header):
        """Replace the file's extended header.
        
        Note that the file's entire data block must be moved if the extended
        header size changes. Setting a new extended header can therefore be
        very time consuming with large files, if the new extended header
        occupies a different number of bytes than the previous one.
        """
        old_ext_header_size = self._extended_header.nbytes
        super(MrcMemmap, self).set_extended_header(extended_header)
        if extended_header.nbytes != old_ext_header_size:
            data_copy = self._data.copy()
            self._close_data()
            self._extended_header = extended_header
            self.header.nsymbt = extended_header.nbytes
            header_nbytes = self.header.nbytes + extended_header.nbytes
            total_nbytes = header_nbytes + data_copy.nbytes

            # Workaround for https://github.com/ccpem/mrcfile/issues/65
            if data_copy.nbytes == 0 and total_nbytes % mmap.ALLOCATIONGRANULARITY == 0:
                # Add one extra byte here to avoid triggering mmap error
                total_nbytes += 1

            self._iostream.truncate(total_nbytes)
            self._open_memmap(data_copy.dtype, data_copy.shape)
            np.copyto(self._data, data_copy)
    
    def flush(self):
        """Flush the header and data arrays to the file buffer."""
        if not self._read_only:
            self._iostream.seek(0)
            self._iostream.write(self.header)
            self._iostream.write(self.extended_header)
            
            # Flushing the file before the mmap makes the mmap flush faster
            self._iostream.flush()
            self._data.flush()
            self._iostream.flush()
            
            # Seek to end of data block so stream is left in the same position
            # as normal
            self._iostream.seek(self._data.nbytes, os.SEEK_CUR)
    
    def _read_data(self):
        """Read the data block from the file.
        
        This method first calculates the parameters needed to read the data
        (block start position, endian-ness, file mode, array shape) and then
        opens the data as a numpy memmap array.
        """
        try:
            dtype = utils.data_dtype_from_header(self.header)
        except ValueError as err:
            if self._permissive:
                warnings.warn("{0} - data block not read".format(err),
                              RuntimeWarning)
                self._data = None
                return
            else:
                raise
        
        shape = utils.data_shape_from_header(self.header)
        
        self._open_memmap(dtype, shape)
    
    def _open_memmap(self, dtype, shape):
        """Open a new memmap array pointing at the file's data block."""
        acc_mode = 'r' if self._read_only else 'r+'
        # Need to use self.header.nsymbt rather than self.extended_header.nbytes because
        # self.extended_header might be None in permissive read mode. Need to convert to
        # Python int (rather than numpy int32) to avoid possible overflow.
        header_nbytes = self.header.nbytes + int(self.header.nsymbt)
        
        self._iostream.flush()
        try:
            self._data = np.memmap(self._iostream,
                                   dtype=dtype,
                                   mode=acc_mode,
                                   offset=header_nbytes,
                                   shape=shape)
        except Exception as ex:
            if self._permissive:
                warnings.warn("Error opening memmap", RuntimeWarning)
                self._data = None
            else:
                raise ex

        # Check if the file is the expected size.
        if self.data is not None:
            file_size = self._get_file_size()
            remaining_file_size = file_size - header_nbytes
            data_size = self.data.nbytes

            # Workaround for https://github.com/ccpem/mrcfile/issues/65
            if data_size == 0 and header_nbytes % mmap.ALLOCATIONGRANULARITY == 0:
                # Expect a file one byte larger here to avoid triggering mmap error
                data_size = 1

            if data_size < remaining_file_size:
                msg = ("MRC file is {0} bytes larger than expected"
                       .format(remaining_file_size - data_size))
                warnings.warn(msg, RuntimeWarning)
    
    def _close_data(self):
        """Delete the existing memmap array, if it exists.
        
        The array is flagged as read-only before deletion, so if a reference to
        it has been kept elsewhere, changes to it should no longer be able to
        change the file contents.
        """
        if self._data is not None:
            self._data.flush()
            self._data.flags.writeable = False
            self._data = None
    
    def _set_new_data(self, data):
        """Override of :meth:`_set_new_data` to handle opening a new memmap and
        copying data into it."""
        # Need to use self.header.nsymbt rather than self.extended_header.nbytes because
        # self.extended_header might be None in permissive read mode. Need to convert to
        # Python int (rather than numpy int32) to avoid possible overflow.
        file_size = self.header.nbytes + int(self.header.nsymbt) + data.nbytes
        self._iostream.truncate(file_size)
        self._open_memmap(data.dtype, data.shape)
        np.copyto(self._data, data, casting='no')
