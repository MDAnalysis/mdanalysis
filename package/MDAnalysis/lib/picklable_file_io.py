# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
"""
Picklable read-only I/O classes --- :mod:`MDAnalysis.lib.picklable_file_io`
===========================================================================

Provide with an interface for pickling read-only IO file object.
These classes are used for further pickling :class:`MDAnalysis.core.universe`
in a object composition approach.

.. autoclass:: FileIOPicklable
   :members:

.. autoclass:: BufferIOPicklable
   :members:

.. autoclass:: TextIOPicklable
   :members:

.. autoclass:: BZ2Picklable
   :members:

.. autoclass:: GzipPicklable
   :members:

.. autofunction:: pickle_open

.. autofunction:: bz2_pickle_open

.. autofunction:: gzip_pickle_open


.. versionadded:: 2.0.0
"""
import io
import os

import bz2
import gzip


class FileIOPicklable(io.FileIO):
    """File object (read-only) that can be pickled.

    This class provides a file-like object (as returned by :func:`open`,
    namely :class:`io.FileIO`) that, unlike standard Python file objects,
    can be pickled. Only read mode is supported.

    When the file is pickled, filename and position of the open file handle in
    the file are saved. On unpickling, the file is opened by filename,
    and the file is seeked to the saved position.
    This means that for a successful unpickle, the original file still has to
    be accessible with its filename.

    Parameters
    ----------
    name : str
        a filename given a text or byte string.

    Example
    -------
    ::

        >>> file = FileIOPicklable(PDB)
        >>> file.readline()
        >>> file_pickled = pickle.loads(pickle.dumps(file))
        >>> print(file.tell(), file_pickled.tell())
            55 55

    See Also
    ---------
    TextIOPicklable
    BufferIOPicklable


    .. versionadded:: 2.0.0
    """
    def __init__(self, name):
        super().__init__(name, mode='r')

    def __getstate__(self):
        return self.name, self.tell()

    def __setstate__(self, args):
        name = args[0]
        super().__init__(name, mode='r')
        self.seek(args[1])


class BufferIOPicklable(io.BufferedReader):
    """A picklable buffer object for read-only FileIO object.

    This class provides a buffered :class:`io.BufferedReader`
    that can be pickled.
    Note that this only works in read mode.

    Parameters
    ----------
    raw : FileIO object

    Example
    -------
    ::

        file = FileIOPicklable('filename')
        buffer_wrapped = BufferIOPicklable(file)

    See Also
    ---------
    FileIOPicklable
    TextIOPicklable


    .. versionadded:: 2.0.0
    """
    def __init__(self, raw):
        super().__init__(raw)
        self.raw_class = raw.__class__

    def __getstate__(self):
        return self.raw_class, self.name, self.tell()

    def __setstate__(self, args):
        raw_class = args[0]
        name = args[1]
        raw = raw_class(name)
        super().__init__(raw)
        self.seek(args[2])


class TextIOPicklable(io.TextIOWrapper):
    """Character and line based picklable file-like object.

    This class provides a file-like :class:`io.TextIOWrapper` object that can
    be pickled. Note that this only works in read mode.

    Parameters
    ----------
    raw : FileIO object

    Example
    -------
    ::

        file = FileIOPicklable('filename')
        text_wrapped = TextIOPicklable(file)

    See Also
    ---------
    FileIOPicklable
    BufferIOPicklable


    .. versionadded:: 2.0.0
    """
    def __init__(self, raw):
        super().__init__(raw)
        self.raw_class = raw.__class__

    def __getstate__(self):
        try:
            name = self.name
        except AttributeError:
            # This is kind of ugly--BZ2File does not save its name.
            name = self.buffer._fp.name
        return self.raw_class, name

    def __setstate__(self, args):
        raw_class = args[0]
        name = args[1]
        # raw_class is used for further expansion this functionality to
        # Gzip files, which also requires a text wrapper.
        raw = raw_class(name)
        super().__init__(raw)


class BZ2Picklable(bz2.BZ2File):
    """File object (read-only) for bzip2 (de)compression that can be pickled.

    This class provides a file-like object (as returned by :func:`bz2.open`,
    namely :class:`bz2.BZ2File`) that, unlike standard Python file objects,
    can be pickled. Only read mode is supported.

    When the file is pickled, filename and position of the open file handle in
    the file are saved. On unpickling, the file is opened by filename,
    and the file is seeked to the saved position.
    This means that for a successful unpickle, the original file still has to
    be accessible with its filename.

    Note
    ----
    This class only supprots opening files in binary mode. If you need to open
    to open a compressed file in text mode, use the :func:`bz2_pickle_open`.

    Parameters
    ----------
    name : str
        a filename given a text or byte string.
    mode : str
        can only be 'r', 'rb' to make pickle work.

    Example
    -------
    ::

        >>> file = BZ2Picklable(XYZ_bz2)
        >>> file.readline()
        >>> file_pickled = pickle.loads(pickle.dumps(file))
        >>> print(file.tell(), file_pickled.tell())
            5 5

    See Also
    ---------
    FileIOPicklable
    BufferIOPicklable
    TextIOPicklable
    GzipPicklable


    .. versionadded:: 2.0.0
    """
    def __init__(self, name, mode='rb'):
        super().__init__(name, mode)

    def __getstate__(self):
        return self._fp.name, self.tell()

    def __setstate__(self, args):
        super().__init__(args[0])
        self.seek(args[1])


class GzipPicklable(gzip.GzipFile):
    """File object (read-only) that can be pickled.

    This class provides a file-like object (as returned by :func:`gzip.open`,
    namely :class:`gzip.GzipFile`) that, unlike standard Python file objects,
    can be pickled. Only read mode is supported.

    When the file is pickled, filename and position of the open file handle in
    the file are saved. On unpickling, the file is opened by filename,
    and the file is seeked to the saved position.
    This means that for a successful unpickle, the original file still has to
    be accessible with its filename.

    Note
    ----
    This class only supprots opening files in binary mode. If you need to open
    to open a compressed file in text mode, use the :func:`gzip_pickle_open`.

    Parameters
    ----------
    name : str
        a filename given a text or byte string.
    mode : str
        can only be 'r', 'rb' to make pickle work.

    Example
    -------
    ::

        >>> file = GzipPicklable(MMTF_gz)
        >>> file.readline()
        >>> file_pickled = pickle.loads(pickle.dumps(file))
        >>> print(file.tell(), file_pickled.tell())
            1218 1218

    See Also
    ---------
    FileIOPicklable
    BufferIOPicklable
    TextIOPicklable
    BZ2Picklable


    .. versionadded:: 2.0.0
    """
    def __init__(self, name, mode='rb'):
        super().__init__(name, mode)

    def __getstate__(self):
        return self.name, self.tell()

    def __setstate__(self, args):
        super().__init__(args[0])
        self.seek(args[1])


def pickle_open(name, mode='rt'):
    """Open file and return a stream with pickle function implemented.

    This function returns either BufferIOPicklable or TextIOPicklable wrapped
    FileIOPicklable object given different reading mode. It can be used as a
    context manager, and replace the built-in :func:`open` function
    in read mode that only returns an unpicklable file object.
    In order to serialize a :class:`MDAnalysis.core.Universe`, this function
    can used to open trajectory/topology files--an object composition approach,
    as opposed to class inheritance, which is more flexible and easier for
    pickle implementation for new readers.

    Note
    ----
    Can be only used with read mode.

    Parameters
    ----------
    name : str
        a filename given a text or byte string.
    mode: {'r', 'rt', 'rb'} (optional)
        'r':  open for reading in text mode;
        'rt': read in text mode (default);
        'rb': read in binary mode;

    Returns
    -------
    stream-like object: BufferIOPicklable or TextIOPicklable
        when mode is 'r' or 'rt', returns TextIOPicklable;
        when mode is 'rb', returns BufferIOPicklable

    Raises
    ------
    ValueError
        if `mode` is not one of the allowed read modes

    Examples
    -------
    open as context manager::

        with pickle_open('filename') as f:
            line = f.readline()

    open as function::

        f = pickle_open('filename')
        line = f.readline()
        f.close()

    See Also
    --------
    :func:`MDAnalysis.lib.util.anyopen`
    :func:`io.open`


    .. versionadded:: 2.0.0
    """
    if mode not in {'r', 'rt', 'rb'}:
        raise ValueError("Only read mode ('r', 'rt', 'rb') \
                         files can be pickled.")
    name = os.fspath(name)
    raw = FileIOPicklable(name)
    if mode == 'rb':
        return BufferIOPicklable(raw)
    elif mode in {'r', 'rt'}:
        return TextIOPicklable(raw)


def bz2_pickle_open(name, mode='rb'):
    """Open a bzip2-compressed file in binary or text mode
    with pickle function implemented.

    This function returns either BZ2Picklable or TextIOPicklable wrapped
    BZ2Picklable object given different reading mode. It can be used as a
    context manager, and replace the built-in :func:`bz2.open` function
    in read mode that only returns an unpicklable file object.

    Note
    ----
    Can be only used with read mode.

    Parameters
    ----------
    name : str
        a filename given a text or byte string.
    mode: {'r', 'rt', 'rb'} (optional)
        'r':  open for reading in binary mode;
        'rt': read in text mode;
        'rb': read in binary mode; (default)

    Returns
    -------
    stream-like object: BZ2Picklable or TextIOPicklable
        when mode is 'rt', returns TextIOPicklable;
        when mode is 'r' or 'rb', returns BZ2Picklable

    Raises
    ------
    ValueError
        if `mode` is not one of the allowed read modes

    Examples
    -------
    open as context manager::

        with bz2_pickle_open('filename') as f:
            line = f.readline()

    open as function::

        f = bz2_pickle_open('filename')
        line = f.readline()
        f.close()

    See Also
    --------
    :func:`io.open`
    :func:`bz2.open`
    :func:`MDAnalysis.lib.util.anyopen`
    :func:`MDAnalysis.lib.picklable_file_io.pickle_open`
    :func:`MDAnalysis.lib.picklable_file_io.gzip_pickle_open`


    .. versionadded:: 2.0.0
    """
    if mode not in {'r', 'rt', 'rb'}:
        raise ValueError("Only read mode ('r', 'rt', 'rb') \
                         files can be pickled.")
    bz_mode = mode.replace("t", "")
    binary_file = BZ2Picklable(name, bz_mode)
    if "t" in mode:
        return TextIOPicklable(binary_file)
    else:
        return binary_file


def gzip_pickle_open(name, mode='rb'):
    """Open a gzip-compressed file in binary or text mode
    with pickle function implemented.

    This function returns either GzipPicklable or TextIOPicklable wrapped
    GzipPicklable object given different reading mode. It can be used as a
    context manager, and replace the built-in :func:`gzip.open` function
    in read mode that only returns an unpicklable file object.

    Note
    ----
    Can be only used with read mode.

    Parameters
    ----------
    name : str
        a filename given a text or byte string.
    mode: {'r', 'rt', 'rb'} (optional)
        'r':  open for reading in binary mode;
        'rt': read in text mode;
        'rb': read in binary mode; (default)

    Returns
    -------
    stream-like object: GzipPicklable or TextIOPicklable
        when mode is 'rt', returns TextIOPicklable;
        when mode is 'r' or 'rb', returns GzipPicklable

    Raises
    ------
    ValueError
        if `mode` is not one of the allowed read modes

    Examples
    -------
    open as context manager::

        with gzip_pickle_open('filename') as f:
            line = f.readline()

    open as function::

        f = gzip_pickle_open('filename')
        line = f.readline()
        f.close()

    See Also
    --------
    :func:`io.open`
    :func:`gzip.open`
    :func:`MDAnalysis.lib.util.anyopen`
    :func:`MDAnalysis.lib.picklable_file_io.pickle_open`
    :func:`MDAnalysis.lib.picklable_file_io.bz2_pickle_open`


    .. versionadded:: 2.0.0
    """
    if mode not in {'r', 'rt', 'rb'}:
        raise ValueError("Only read mode ('r', 'rt', 'rb') \
                         files can be pickled.")
    gz_mode = mode.replace("t", "")
    binary_file = GzipPicklable(name, gz_mode)
    if "t" in mode:
        return TextIOPicklable(binary_file)
    else:
        return binary_file
