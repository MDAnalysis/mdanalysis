# Copyright (c) 2016, Science and Technology Facilities Council
# This software is distributed under a BSD licence. See LICENSE.txt.
"""
bzip2mrcfile
------------

Module which exports the :class:`Bzip2MrcFile` class.

Classes:
    :class:`Bzip2MrcFile`: An object which represents a bzip2-compressed MRC
    file.

"""

# Import Python 3 features for future-proofing
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import bz2
import os

from .mrcfile import MrcFile


class Bzip2MrcFile(MrcFile):
    
    """:class:`~mrcfile.mrcfile.MrcFile` subclass for handling bzip2-compressed
    files.
    
    Usage is the same as for :class:`~mrcfile.mrcfile.MrcFile`.
    
    """
    
    def __repr__(self):
        return "Bzip2MrcFile('{0}', mode='{1}')".format(self._fname,
                                                        self._mode)
    
    def _open_file(self, name):
        """Override _open_file() to open a bzip2 file."""
        self._fname = name
        if 'w' in self._mode and not os.path.exists(name):
            open(name, mode='w').close()
        self._iostream = bz2.BZ2File(name, mode='r')
    
    def _read(self, header_only=False):
        """Override _read() to ensure bzip2 file is in read mode."""
        self._ensure_readable_bzip2_stream()
        super(Bzip2MrcFile, self)._read(header_only)
    
    def _ensure_readable_bzip2_stream(self):
        """Make sure _iostream is a bzip2 stream that can be read."""
        if hasattr(self._iostream, "readable"):
            # Python 3
            readable = self._iostream.readable()
        else:
            # Python 2
            readable = (self._iostream.mode[0] == 'r')
        if not readable:
            self._iostream.close()
            self._iostream = bz2.BZ2File(self._fname, mode='r')
    
    def _get_file_size(self):
        """Override _get_file_size() to ensure stream is readable first."""
        self._ensure_readable_bzip2_stream()
        return super(Bzip2MrcFile, self)._get_file_size()

    def _read_bytearray_from_stream(self, number_of_bytes):
        """Override because BZ2File in Python 2 does not support
        :meth:`~io.BufferedIOBase.readinto`."""
        if hasattr(self._iostream, "readinto"):
            # Python 3 - BZ2File supports ``readinto()`` so we just use the normal implementation
            return super(Bzip2MrcFile, self)._read_bytearray_from_stream(number_of_bytes)
        else:
            # Python 2 - need to read as bytes then copy to a bytearray
            result_bytes = self._iostream.read(number_of_bytes)
            return bytearray(result_bytes), len(result_bytes)
    
    def flush(self):
        """Override :meth:`~mrcfile.mrcinterpreter.MrcInterpreter.flush` since
        BZ2File objects need special handling.
        """
        if not self._read_only:
            self._iostream.close()
            self._iostream = bz2.BZ2File(self._fname, mode='w')
            
            # Arrays converted to bytes so bz2 can calculate sizes correctly
            self._iostream.write(self.header.tobytes())
            self._iostream.write(self.extended_header.tobytes())
            self._iostream.write(self.data.tobytes())
            # no equivalent for flush() with BZ2File
