# Copyright (c) 2016, Science and Technology Facilities Council
# This software is distributed under a BSD licence. See LICENSE.txt.

"""
utils
-----

Utility functions used by the other modules in the mrcfile package.

Functions
---------

* :func:`data_dtype_from_header`: Work out the data :class:`dtype
  <numpy.dtype>` from an MRC header.
* :func:`data_shape_from_header`: Work out the data array shape from an MRC
  header
* :func:`mode_from_dtype`: Convert a :class:`numpy dtype <numpy.dtype>` to an
  MRC mode number.
* :func:`dtype_from_mode`: Convert an MRC mode number to a :class:`numpy dtype
  <numpy.dtype>`.
* :func:`pretty_machine_stamp`: Get a nicely-formatted string from a machine
  stamp.
* :func:`machine_stamp_from_byte_order`: Get a machine stamp from a byte order
  indicator.
* :func:`byte_orders_equal`: Compare two byte order indicators for equal
  endianness.
* :func:`normalise_byte_order`: Convert a byte order indicator to ``<`` or
  ``>``.
* :func:`spacegroup_is_volume_stack`: Identify if a space group number
  represents a volume stack.

"""

# Import Python 3 features for future-proofing
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import string
import sys

import numpy as np

from .constants import IMAGE_STACK_SPACEGROUP


def data_dtype_from_header(header):
    """Return the data dtype indicated by the given header.
    
    This function calls :func:`dtype_from_mode` to get the basic dtype, and
    then makes sure that the byte order of the new dtype matches the byte order
    of the header's ``mode`` field.
    
    Args:
        header: An MRC header as a :class:`numpy record array
            <numpy.recarray>`.
    
    Returns:
        The :class:`numpy dtype <numpy.dtype>` object for the data array
        corresponding to the given header.
    
    Raises:
        :exc:`ValueError`: If there is no corresponding dtype for the given
            mode.
    """
    mode = header.mode
    return dtype_from_mode(mode).newbyteorder(mode.dtype.byteorder)


def data_shape_from_header(header):
    """Return the data shape indicated by the given header.
    
    Args:
        header: An MRC header as a :class:`numpy record array
            <numpy.recarray>`.
    
    Returns:
        The shape tuple for the data array corresponding to the given header.
    """
    nx = int(header.nx)
    ny = int(header.ny)
    nz = int(header.nz)
    mz = int(header.mz)
    
    if spacegroup_is_volume_stack(header.ispg):
        shape = (nz // mz, mz, ny, nx)
    elif header.ispg == IMAGE_STACK_SPACEGROUP and nz == 1:
        # Use a 2D array for a single image
        shape = (ny, nx)
    else:
        shape = (nz, ny, nx)
    
    return shape


_dtype_to_mode = dict(f2=12, f4=2, i1=0, i2=1, u1=6, u2=6, c8=4)

def mode_from_dtype(dtype):
    """Return the MRC mode number corresponding to the given :class:`numpy
    dtype <numpy.dtype>`.
    
    The conversion is as follows:
    
    * float16   -> mode 12
    * float32   -> mode 2
    * int8      -> mode 0
    * int16     -> mode 1
    * uint8     -> mode 6 (data will be widened to 16 bits in the file)
    * uint16    -> mode 6
    * complex64 -> mode 4
    
    Note that there is no numpy dtype which corresponds to MRC mode 3.
    
    Args:
        dtype: A :class:`numpy dtype <numpy.dtype>` object.
    
    Returns:
        The MRC mode number.
    
    Raises:
        :exc:`ValueError`: If there is no corresponding MRC mode for the given
            dtype.
    """
    kind_and_size = dtype.kind + str(dtype.itemsize)
    if kind_and_size in _dtype_to_mode:
        return _dtype_to_mode[kind_and_size]
    raise ValueError("dtype '{0}' cannot be converted "
                     "to an MRC file mode".format(dtype))


_mode_to_dtype = { 0: np.int8,
                   1: np.int16,
                   2: np.float32,
                   4: np.complex64,
                   6: np.uint16,
                   12: np.float16 }

def dtype_from_mode(mode):
    """Return the :class:`numpy dtype <numpy.dtype>` corresponding to the given
    MRC mode number.
    
    The mode parameter may be given as a Python scalar, numpy scalar or
    single-item numpy array.
    
    The conversion is as follows:
    
    * mode 0 -> int8
    * mode 1 -> int16
    * mode 2 -> float32
    * mode 4 -> complex64
    * mode 6 -> uint16
    * mode 12 -> float16
    
    Note that modes 3 and 101 are not supported as there is no matching numpy dtype.
    
    Args:
        mode: The MRC mode number. This may be given as any type which can be
            converted to an int, for example a Python scalar (``int`` or
            ``float``), a numpy scalar or a single-item numpy array.
    
    Returns:
        The :class:`numpy dtype <numpy.dtype>` object corresponding to the
        given mode.
    
    Raises:
        :exc:`ValueError`: If there is no corresponding dtype for the given
            mode, or if ``mode`` is an array and does not contain exactly one
            item.
    """
    if isinstance(mode, np.ndarray):
        if mode.size != 1:
            raise ValueError("Mode array should contain exactly one item")
        mode = mode.item()
    if mode in _mode_to_dtype:
        return np.dtype(_mode_to_dtype[mode])
    else:
        raise ValueError("Unrecognised mode '{0}'".format(mode))


def pretty_machine_stamp(machst):
    """Return a human-readable hex string for a machine stamp."""
    return " ".join("0x{:02x}".format(byte) for byte in machst)


def byte_order_from_machine_stamp(machst):
    """Return the byte order corresponding to the given machine stamp.
    
    Args:
        machst: The machine stamp, as a :class:`bytearray` or a :class:`numpy
            array <numpy.ndarray>` of bytes.
    
    Returns:
        ``<`` if the machine stamp represents little-endian data, or ``>`` if
        it represents big-endian.
    
    Raises:
        :exc:`ValueError`: If the machine stamp is invalid.
    """
    if machst[0] == 0x44 and machst[1] in (0x44, 0x41):
        return '<'
    elif (machst[0] == 0x11 and machst[1] == 0x11):
        return '>'
    else:
        pretty_bytes = pretty_machine_stamp(machst)
        raise ValueError("Unrecognised machine stamp: " + pretty_bytes)


_byte_order_to_machine_stamp = {'<': bytearray((0x44, 0x44, 0, 0)),
                                '>': bytearray((0x11, 0x11, 0, 0))}

def machine_stamp_from_byte_order(byte_order='='):
    """Return the machine stamp corresponding to the given byte order
    indicator.
    
    Args:
        byte_order: The byte order indicator: one of ``=``, ``<`` or ``>``, as
            defined and used by numpy dtype objects.
    
    Returns:
        The machine stamp which corresponds to the given byte order, as a
        :class:`bytearray`. This will be either ``(0x44, 0x44, 0, 0)`` for
        little-endian or ``(0x11, 0x11, 0, 0)`` for big-endian. If the given
        byte order indicator is ``=``, the native byte order is used.
    
    Raises:
        :exc:`ValueError`: If the byte order indicator is unrecognised.
    """
    # If byte order is '=', replace it with the system-native order
    byte_order = normalise_byte_order(byte_order)
    return _byte_order_to_machine_stamp[byte_order]

def byte_orders_equal(a, b):
    """Work out if the byte order indicators represent the same endianness.
    
    Args:
        a: The first byte order indicator: one of ``=``, ``<`` or ``>``, as
            defined and used by :class:`numpy dtype <numpy.dtype>` objects.
        b: The second byte order indicator.
    
    Returns:
        :data:`True` if the byte order indicators represent the same
        endianness.
    
    Raises:
        :exc:`ValueError`: If the byte order indicator is not recognised.
    """
    return normalise_byte_order(a) == normalise_byte_order(b)

def normalise_byte_order(byte_order):
    """Convert a numpy byte order indicator to one of ``<`` or ``>``.
    
    Args:
        byte_order: One of ``=``, ``<`` or ``>``.
    
    Returns:
        ``<`` if the byte order indicator represents little-endian data, or
        ``>`` if it represents big-endian. Therefore on a little-endian
        machine, ``=`` will be converted to ``<``, but on a big-endian machine
        it will be converted to ``>``.
    
    Raises:
        :exc:`ValueError`: If ``byte_order`` is not one of ``=``, ``<`` or
            ``>``.
    """
    if byte_order not in ('<', '>', '='):
        raise ValueError("Unrecognised byte order indicator '{0}'"
                         .format(byte_order))
    if byte_order == '=':
        return '<' if sys.byteorder == 'little' else '>'
    return byte_order

def spacegroup_is_volume_stack(ispg):
    """Identify if the given space group number represents a volume stack.
    
    Args:
        ispg: The space group number, as an integer, numpy scalar or single-
            element numpy array.
    
    Returns:
        :data:`True` if the space group number is in the range 401--630.
    """
    return 401 <= ispg <= 630


printable_chars = ' ' + string.ascii_letters + string.digits + string.punctuation


def is_printable_ascii(string_):
    """Check if a string is entirely composed of printable ASCII characters."""
    try:
        # Python 3 version
        return str.isprintable(string_) and str.isascii(string_)
    except AttributeError:
        # Probably Python 2, fall back to checking characters individually
        try:
            return all(char in printable_chars for char in string_)
        except UnicodeDecodeError:
            return False


def printable_string_from_bytes(bytes_):
    """Convert bytes into a printable ASCII string by removing non-printable characters.
    """
    string_ = bytes.decode(bytes_, encoding='ascii', errors='ignore')
    if not is_printable_ascii(string_):
        string_ = "".join(list(s for s in string_ if is_printable_ascii(s)))
    return string_


def bytes_from_string(string_):
    """Convert a string to bytes.

    Even though this is a one-liner, the details are tricky to get right so things work
    properly in both Python 2 and 3. It's broken out as a separate function so it can be
    thoroughly tested.

    Raises:
        UnicodeError: If the input contains non-ASCII characters.
    """
    return str.encode(str(string_), encoding='ascii', errors='strict')
