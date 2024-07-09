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
r"""
Helper functions --- :mod:`MDAnalysis.lib.util`
====================================================

Small helper functions that don't fit anywhere else.

.. versionchanged:: 0.11.0
   Moved mathematical functions into lib.mdamath

.. versionchanged::2.0.0
   The following aliases, that existed for compatibility with python versions
   older than 3.6, were removed: `callable` for the built-in of the same name,
   `PathLike` for :class:`os.PathLike`, and `bz_open` for :func:`bz2.open`.


Files and directories
---------------------

.. autofunction:: filename
.. autofunction:: openany
.. autofunction:: anyopen
.. autofunction:: greedy_splitext
.. autofunction:: which
.. autofunction:: realpath
.. autofunction:: get_ext
.. autofunction:: check_compressed_format
.. autofunction:: format_from_filename_extension
.. autofunction:: guess_format

Streams
-------

Many of the readers are not restricted to just reading files. They can
also use gzip-compressed or bzip2-compressed files (through the
internal use of :func:`openany`). It is also possible to provide more
general streams as inputs, such as a :class:`io.StringIO`
instances (essentially, a memory buffer) by wrapping these instances
into a :class:`NamedStream`. This :class:`NamedStream` can then be
used in place of an ordinary file name (typically, with a
class:`~MDAnalysis.core.universe.Universe` but it is also possible to
*write* to such a stream using :func:`MDAnalysis.Writer`).

.. rubric: Examples

In the following example, we use a PDB stored as a string ``pdb_s``::

   import MDAnalysis
   from MDAnalysis.lib.util import NamedStream
   from io import StringIO

   pdb_s = "TITLE     Lonely Ion\\nATOM      1  NA  NA+     1      81.260  64.982  10.926  1.00  0.00\\n"
   u = MDAnalysis.Universe(NamedStream(StringIO(pdb_s), "ion.pdb"))
   print(u)
   #  <Universe with 1 atoms>
   print(u.atoms.positions)
   # [[ 81.26000214  64.98200226  10.92599964]]

It is important to provide a proper pseudo file name with the correct extension
(".pdb") to :class:`NamedStream` because the file type recognition uses the
extension of the file name to determine the file format or alternatively
provide the ``format="pdb"`` keyword argument to the
:class:`~MDAnalysis.core.universe.Universe`.

The use of streams becomes more interesting when MDAnalysis is used as glue
between different analysis packages and when one can arrange things so that
intermediate frames (typically in the PDB format) are not written to disk but
remain in memory via e.g. :class:`io.StringIO` buffers.


.. The following does *not* work because most readers need to
.. reopen files, which is not possible with http streams. Might
.. need to implement a buffer.
..
.. Read a test LAMMPS data file from the MDAnalysis repository::
..
..   import MDAnalysis
..   from MDAnalysis.lib.util import NamedStream
..   import urllib2
..   URI = "https://mdanalysis.googlecode.com/git-history/develop/testsuite/MDAnalysisTests/data/mini.data"
..   urldata = NamedStream(urllib2.urlopen(URI), "mini.data")
..   u = MDAnalysis.Universe(urldata)

.. Note::  A remote connection created by :func:`urllib2.urlopen` is not seekable
           and therefore will often not work as an input. But try it...

.. autoclass:: NamedStream
   :members:

.. autofunction:: isstream

Containers and lists
--------------------

.. autofunction:: iterable
.. autofunction:: asiterable
.. autofunction:: hasmethod
.. autoclass:: Namespace

Arrays
------

.. autofunction:: unique_int_1d(values)
.. autofunction:: unique_rows
.. autofunction:: blocks_of
.. autofunction:: group_same_or_consecutive_integers

File parsing
------------

.. autoclass:: FORTRANReader
   :members:
.. autodata:: FORTRAN_format_regex

Data manipulation and handling
------------------------------

.. autofunction:: fixedwidth_bins
.. autofunction:: get_weights
.. autofunction:: ltruncate_int
.. autofunction:: flatten_dict

Strings
-------

.. autofunction:: convert_aa_code
.. autofunction:: parse_residue
.. autofunction:: conv_float

Class decorators
----------------

.. autofunction:: cached
.. autofunction:: store_init_arguments

Function decorators
-------------------

.. autofunction:: static_variables
.. autofunction:: warn_if_not_unique
.. autofunction:: check_coords
.. autofunction:: check_atomgroup_not_empty

Code management
---------------

.. autofunction:: deprecate
.. autoclass:: _Deprecate
.. autofunction:: dedent_docstring

Data format checks
------------------

.. autofunction:: check_box

.. Rubric:: Footnotes

.. [#NamedStreamClose] The reason why :meth:`NamedStream.close` does
   not close a stream by default (but just rewinds it to the
   beginning) is so that one can use the class :class:`NamedStream` as
   a drop-in replacement for file names, which are often re-opened
   (e.g. when the same file is used as a topology and coordinate file
   or when repeatedly iterating through a trajectory in some
   implementations). The ``close=True`` keyword can be supplied in
   order to make :meth:`NamedStream.close` actually close the
   underlying stream and ``NamedStream.close(force=True)`` will also
   close it.
"""
import sys

__docformat__ = "restructuredtext en"


import os
import os.path
import errno
from contextlib import contextmanager
import bz2
import gzip
import re
import io
import warnings
import functools
from functools import wraps
import textwrap
import weakref
import itertools

import numpy as np

from numpy.testing import assert_equal
import inspect

from .picklable_file_io import pickle_open, bz2_pickle_open, gzip_pickle_open

from ..exceptions import StreamWarning, DuplicateWarning

try:
    from ._cutil import unique_int_1d # pylint: disable=unused-import
except ImportError:
    raise ImportError("MDAnalysis not installed properly. "
                      "This can happen if your C extensions "
                      "have not been built.")


def int_array_is_sorted(array):
    mask = array[:-1] <= array[1:]
    try:
        return mask[0] and mask.argmin() == 0
    except IndexError:
        # Empty arrays are sorted, I guess...
        return True


def unique_int_1d_unsorted(array):
    values, indices = np.unique(array, return_index=True)
    return array[np.sort(indices)]


def filename(name, ext=None, keep=False):
    """Return a new name that has suffix attached; replaces other extensions.

    Parameters
    ----------
    name : str or NamedStream
        filename; extension is replaced unless ``keep=True``;
        `name` can also be a :class:`NamedStream` (and its
        :attr:`NamedStream.name` will be changed accordingly)
    ext : None or str
        extension to use in the new filename
    keep : bool
        - ``False``: replace existing extension with `ext`;
        - ``True``: keep old extension if one existed


    .. versionchanged:: 0.9.0
       Also permits :class:`NamedStream` to pass through.
    """
    if ext is not None:
        ext = ext.lower()
        if not ext.startswith(os.path.extsep):
            ext = os.path.extsep + ext
        root, origext = os.path.splitext(name)
        if not keep or len(origext) == 0:
            newname = root + ext
            if isstream(name):
                name.name = newname
            else:
                name = newname
    return name if isstream(name) else str(name)


@contextmanager
def openany(datasource, mode='rt', reset=True):
    """Context manager for :func:`anyopen`.

    Open the `datasource` and close it when the context of the :keyword:`with`
    statement exits.

    `datasource` can be a filename or a stream (see :func:`isstream`). A stream
    is reset to its start if possible (via :meth:`~io.IOBase.seek` or
    :meth:`~cString.StringIO.reset`).

    The advantage of this function is that very different input sources
    ("streams") can be used for a "file", ranging from files on disk (including
    compressed files) to open file objects to sockets and strings---as long as
    they have a file-like interface.

    Parameters
    ----------
    datasource : a file or a stream
    mode : {'r', 'w'} (optional)
        open in r(ead) or w(rite) mode
    reset : bool (optional)
        try to read (`mode` 'r') the stream from the start [``True``]

    Examples
    --------
    Open a gzipped file and process it line by line::

        with openany("input.pdb.gz") as pdb:
            for line in pdb:
                if line.startswith('ATOM'):
                    print(line)

    Open a URL and read it::

       import urllib2
       with openany(urllib2.urlopen("https://www.mdanalysis.org/")) as html:
           print(html.read())


    See Also
    --------
    :func:`anyopen`
    """
    stream = anyopen(datasource, mode=mode, reset=reset)
    try:
        yield stream
    finally:
        stream.close()


def anyopen(datasource, mode='rt', reset=True):
    """Open datasource (gzipped, bzipped, uncompressed) and return a stream.

    `datasource` can be a filename or a stream (see :func:`isstream`). By
    default, a stream is reset to its start if possible (via
    :meth:`~io.IOBase.seek` or :meth:`~cString.StringIO.reset`).

    If possible, the attribute ``stream.name`` is set to the filename or
    "<stream>" if no filename could be associated with the *datasource*.

    Parameters
    ----------
    datasource
        a file (from :class:`file` or :func:`open`) or a stream (e.g. from
        :func:`urllib2.urlopen` or :class:`io.StringIO`)
    mode: {'r', 'w', 'a'} (optional)
        Open in r(ead), w(rite) or a(ppend) mode. This string is directly
        passed to the file opening handler (either Python's openfe, bz2, or
        gzip). More complex modes are supported if the file opening handler
        supports it.
    reset: bool (optional)
        try to read (`mode` 'r') the stream from the start

    Returns
    -------
    file-like object

    See Also
    --------
    :func:`openany`
      to be used with the :keyword:`with` statement.


    .. versionchanged:: 0.9.0
       Only returns the ``stream`` and tries to set ``stream.name = filename`` instead of the previous
       behavior to return a tuple ``(stream, filename)``.

    .. versionchanged:: 2.0.0
       New read handlers support pickle functionality
       if `datasource` is a filename.
       They return a custom picklable file stream in
       :class:`MDAnalysis.lib.picklable_file_io`.

    """
    read_handlers = {'bz2': bz2_pickle_open,
                     'gz': gzip_pickle_open,
                     '': pickle_open}
    write_handlers = {'bz2': bz2.open,
                      'gz': gzip.open,
                      '': open}

    if mode.startswith('r'):
        if isstream(datasource):
            stream = datasource
            try:
                filename = str(stream.name)  # maybe that does not always work?
            except AttributeError:
                filename = "<stream>"
            if reset:
                try:
                    stream.reset()
                except (AttributeError, IOError):
                    try:
                        stream.seek(0)
                    except (AttributeError, IOError):
                        warnings.warn("Stream {0}: not guaranteed to be at the beginning."
                                      "".format(filename),
                                      category=StreamWarning)
        else:
            stream = None
            filename = datasource
            for ext in ('bz2', 'gz', ''):  # file == '' should be last
                openfunc = read_handlers[ext]
                stream = _get_stream(datasource, openfunc, mode=mode)
                if stream is not None:
                    break
            if stream is None:
                raise IOError(errno.EIO, "Cannot open file or stream in mode={mode!r}.".format(**vars()), repr(filename))
    elif mode.startswith('w') or mode.startswith('a'):  # append 'a' not tested...
        if isstream(datasource):
            stream = datasource
            try:
                filename = str(stream.name)  # maybe that does not always work?
            except AttributeError:
                filename = "<stream>"
        else:
            stream = None
            filename = datasource
            name, ext = os.path.splitext(filename)
            if ext.startswith('.'):
                ext = ext[1:]
            if not ext in ('bz2', 'gz'):
                ext = ''  # anything else but bz2 or gz is just a normal file
            openfunc = write_handlers[ext]
            stream = openfunc(datasource, mode=mode)
            if stream is None:
                raise IOError(errno.EIO, "Cannot open file or stream in mode={mode!r}.".format(**vars()), repr(filename))
    else:
        raise NotImplementedError("Sorry, mode={mode!r} is not implemented for {datasource!r}".format(**vars()))
    try:
        stream.name = filename
    except (AttributeError, TypeError):
        pass  # can't set name (e.g. io.StringIO)
    return stream


def _get_stream(filename, openfunction=open, mode='r'):
    """Return open stream if *filename* can be opened with *openfunction* or else ``None``."""
    try:
        stream = openfunction(filename, mode=mode)
    except (IOError, OSError) as err:
        # An exception might be raised due to two reasons, first the openfunction is unable to open the file, in this
        # case we have to ignore the error and return None. Second is when openfunction can't open the file because
        # either the file isn't there or the permissions don't allow access.
        if errno.errorcode[err.errno] in ['ENOENT', 'EACCES']:
            raise sys.exc_info()[1] from err
        return None
    if mode.startswith('r'):
        # additional check for reading (eg can we uncompress) --- is this needed?
        try:
            stream.readline()
        except IOError:
            stream.close()
            stream = None
        except:
            stream.close()
            raise
        else:
            stream.close()
            stream = openfunction(filename, mode=mode)
    return stream


def greedy_splitext(p):
    """Split extension in path *p* at the left-most separator.

    Extensions are taken to be separated from the filename with the
    separator :data:`os.extsep` (as used by :func:`os.path.splitext`).

    Arguments
    ---------
    p : str
       path

    Returns
    -------
    (root, extension) : tuple
          where ``root`` is the full path and filename with all
          extensions removed whereas ``extension`` is the string of
          all extensions.

    Example
    -------

    >>> from MDAnalysis.lib.util import greedy_splitext
    >>> greedy_splitext("/home/joe/protein.pdb.bz2")
    ('/home/joe/protein', '.pdb.bz2')

    """
    path, root = os.path.split(p)
    extension = ''
    while True:
        root, ext = os.path.splitext(root)
        extension = ext + extension
        if not ext:
            break
    return os.path.join(path, root), extension


def hasmethod(obj, m):
    """Return ``True`` if object *obj* contains the method *m*."""
    return hasattr(obj, m) and callable(getattr(obj, m))


def isstream(obj):
    """Detect if `obj` is a stream.

    We consider anything a stream that has the methods

    - ``close()``

    and either set of the following

    - ``read()``, ``readline()``, ``readlines()``
    - ``write()``, ``writeline()``, ``writelines()``

    Parameters
    ----------
    obj : stream or str

    Returns
    -------
    bool
        ``True`` if `obj` is a stream, ``False`` otherwise

    See Also
    --------
    :mod:`io`


    .. versionadded:: 0.9.0
    """
    signature_methods = ("close",)
    alternative_methods = (
        ("read", "readline", "readlines"),
        ("write", "writeline", "writelines"))

    # Must have ALL the signature methods
    for m in signature_methods:
        if not hasmethod(obj, m):
            return False
    # Must have at least one complete set of alternative_methods
    alternative_results = [
        np.all([hasmethod(obj, m) for m in alternatives])
        for alternatives in alternative_methods]
    return np.any(alternative_results)


def which(program):
    """Determine full path of executable `program` on :envvar:`PATH`.

    (Jay at http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python)

    Parameters
    ----------
    programe : str
       name of the executable

    Returns
    -------
    path : str or None
       absolute path to the executable if it can be found, else ``None``


    .. deprecated:: 2.7.0
       This method is deprecated and will be removed in version 3.0.0.
       Please use shutil.which instead.
    """
    # Can't use decorator because it's declared after this method
    wmsg = ("This method is deprecated as of MDAnalysis version 2.7.0 "
            "and will be removed in version 3.0.0. Please use shutil.which "
            "instead.")
    warnings.warn(wmsg, DeprecationWarning)

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        real_program = realpath(program)
        if is_exe(real_program):
            return real_program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


@functools.total_ordering
class NamedStream(io.IOBase, os.PathLike):
    """Stream that also provides a (fake) name.

    By wrapping a stream `stream` in this class, it can be passed to
    code that uses inspection of the filename to make decisions. For
    instance. :func:`os.path.split` will work correctly on a
    :class:`NamedStream`.

    The class can be used as a context manager.

    :class:`NamedStream` is derived from :class:`io.IOBase` (to indicate that
    it is a stream). Many operations that normally expect a string will also
    work with a :class:`NamedStream`; for instance, most of the functions in
    :mod:`os.path` will work with the exception of :func:`os.path.expandvars`
    and :func:`os.path.expanduser`, which will return the :class:`NamedStream`
    itself instead of a string if no substitutions were made.

    Example
    -------
    Wrap a :class:`io.StringIO` instance to write to::

      from io import StringIO
      import os.path
      stream = StringIO()
      f = NamedStream(stream, "output.pdb")
      print(os.path.splitext(f))

    Wrap a :class:`file` instance to read from::

      stream = open("input.pdb")
      f = NamedStream(stream, stream.name)

    Use as a context manager (closes stream automatically when the
    :keyword:`with` block is left)::

      with NamedStream(open("input.pdb"), "input.pdb") as f:
         # use f
         print f.closed  # --> False
         # ...
      print f.closed     # --> True

    Note
    ----
    This class uses its own :meth:`__getitem__` method so if `stream`
    implements :meth:`stream.__getitem__` then that will be masked and
    this class should not be used.

    Warning
    -------
    By default, :meth:`NamedStream.close` will **not close the
    stream** but instead :meth:`~NamedStream.reset` it to the
    beginning. [#NamedStreamClose]_ Provide the ``force=True`` keyword
    to :meth:`NamedStream.close` to always close the stream.

    """

    def __init__(self, stream, filename, reset=True, close=False):
        """Initialize the :class:`NamedStream` from a `stream` and give it a `name`.

        The constructor attempts to rewind the stream to the beginning unless
        the keyword `reset` is set to ``False``. If rewinding fails, a
        :class:`MDAnalysis.StreamWarning` is issued.

        Parameters
        ----------
        stream : stream
            an open stream (e.g. :class:`file` or :class:`io.StringIO`)
        filename : str
            the filename that should be associated with the stream
        reset : bool (optional)
            start the stream from the beginning (either :meth:`reset` or :meth:`seek`)
            when the class instance is constructed
        close : bool (optional)
            close the stream when a :keyword:`with` block exits or when
            :meth:`close` is called; note that the default is **not to close
            the stream**

        Notes
        -----
        By default, this stream will *not* be closed by :keyword:`with` and
        :meth:`close` (see there) unless the `close` keyword is set to
        ``True``.


        .. versionadded:: 0.9.0
        """
        # constructing the class from an instance of itself has weird behavior
        # on __del__ and super on python 3. Let's warn the user and ensure the
        # class works normally.
        if isinstance(stream, NamedStream):
            warnings.warn("Constructed NamedStream from a NamedStream",
                          RuntimeWarning)
            stream = stream.stream
        self.stream = stream
        self.name = filename
        self.close_stream = close
        if reset:
            self.reset()

    def reset(self):
        """Move to the beginning of the stream"""
        # try to rewind
        try:
            self.stream.reset()  # e.g. StreamIO
        except (AttributeError, IOError):
            try:
                self.stream.seek(0)  # typical file objects
            except (AttributeError, IOError):
                warnings.warn("NamedStream {0}: not guaranteed to be at the beginning."
                              "".format(self.name),
                              category=StreamWarning)

    # access the stream
    def __getattr__(self, x):
        try:
            return getattr(self.stream, x)
        except AttributeError:
            return getattr(self.name, x)

    def __iter__(self):
        return iter(self.stream)

    def __next__(self):
        return self.stream.__next__()

    def __enter__(self):
        # do not call the stream's __enter__ because the stream is already open
        return self

    def __exit__(self, *args):
        # NOTE: By default (close=False) we only reset the stream and NOT close it; this makes
        #       it easier to use it as a drop-in replacement for a filename that might
        #       be opened repeatedly (at least in MDAnalysis)
        #try:
        #    return self.stream.__exit__(*args)
        #except AttributeError:
        #    super(NamedStream, self).__exit__(*args)
        self.close()

    def __fspath__(self):
        return self.name

    # override more IOBase methods, as these are provided by IOBase and are not
    # caught with __getattr__ (ugly...)
    def close(self, force=False):
        """Reset or close the stream.

        If :attr:`NamedStream.close_stream` is set to ``False`` (the default)
        then this method will *not close the stream* and only :meth:`reset` it.

        If the *force* = ``True`` keyword is provided, the stream will be
        closed.

        .. Note:: This ``close()`` method is non-standard. ``del NamedStream``
                  always closes the underlying stream.

        """
        if self.closed:
            return

        if self.close_stream or force:
            try:
                return self.stream.close()
            except AttributeError:
                return super(NamedStream, self).close()
        else:
            self.flush()
            self.reset()

    def __del__(self):
        """Always closes the stream."""
        self.close(force=True)

    @property
    def closed(self):
        """``True`` if stream is closed."""
        try:
            return self.stream.closed
        except AttributeError:
            return super(NamedStream, self).closed

    def seek(self, offset, whence=os.SEEK_SET):
        """Change the stream position to the given byte `offset` .

        Parameters
        ----------
        offset : int
             `offset` is interpreted relative to the position
             indicated by `whence`.
        whence : {0, 1, 2} (optional)
             Values for `whence` are:

               - :data:`io.SEEK_SET` or 0 – start of the stream (the default); `offset`
                 should be zero or positive
               - :data:`io.SEEK_CUR` or 1 – current stream position; `offset` may be
                 negative
               - :data:`io.SEEK_END` or 2 – end of the stream; `offset` is usually
                 negative

        Returns
        -------
        int
            the new absolute position in bytes.

        """
        try:
            return self.stream.seek(offset, whence)  # file.seek: no kw
        except AttributeError:
            return super(NamedStream, self).seek(offset, whence)

    def tell(self):
        """Return the current stream position."""
        try:
            return self.stream.tell()
        except AttributeError:
            return super(NamedStream, self).tell()

    def truncate(self, *size):
        """Truncate the stream's size to `size`.

        Parameters
        ----------
        size : int (optional)
             The `size` defaults to the current position (if no `size` argument
             is supplied). The current file position is not changed.

        """
        try:
            return self.stream.truncate(*size)
        except AttributeError:
            return super(NamedStream, self).truncate(*size)

    def seekable(self):
        """Return ``True`` if the stream supports random access.

        If ``False``, :meth:`seek`, :meth:`tell` and :meth:`truncate` will
        raise :exc:`IOError`.

        """
        try:
            return self.stream.seekable()
        except AttributeError:
            return super(NamedStream, self).seekable()

    def readable(self):
        """Return ``True`` if the stream can be read from.

        If ``False``, :meth:`read` will raise :exc:`IOError`.
        """
        try:
            return self.stream.readable()
        except AttributeError:
            return super(NamedStream, self).readable()

    def writable(self):
        """Return ``True`` if the stream can be written to.

        If ``False``, :meth:`write` will raise :exc:`IOError`.
        """
        try:
            return self.stream.writable()
        except AttributeError:
            return super(NamedStream, self).writable()

    def flush(self):
        """Flush the write buffers of the stream if applicable.

        This does nothing for read-only and non-blocking streams. For file
        objects one also needs to call :func:`os.fsync` to write contents to
        disk.
        """
        try:
            return self.stream.flush()
        except AttributeError:
            return super(NamedStream, self).flush()

    def fileno(self):
        """Return the underlying file descriptor (an integer) of the stream if it exists.

        An :exc:`IOError` is raised if the IO object does not use a file descriptor.
        """
        try:
            return self.stream.fileno()
        except AttributeError:
            # IOBase.fileno does not raise IOError as advertised so we do this here
            errmsg = "This NamedStream does not use a file descriptor."
            raise IOError(errmsg) from None

    def readline(self):
        try:
            return self.stream.readline()
        except AttributeError:
            return super(NamedStream, self).readline()

    # fake the important parts of the string API
    # (other methods such as rfind() are automatically dealt with via __getattr__)
    def __getitem__(self, x):
        return self.name[x]

    def __eq__(self, x):
        return self.name == x

    def __ne__(self, x):
        return not self == x

    def __lt__(self, x):
        return self.name < x

    def __len__(self):
        return len(self.name)

    def __add__(self, x):
        return self.name + x

    def __radd__(self, x):
        return x + self.name

    def __mul__(self, x):
        return self.name * x

    __rmul__ = __mul__

    def __format__(self, format_spec):
        return self.name.format(format_spec)

    def __str__(self):
        return self.name

    def __repr__(self):
        return "<NamedStream({0}, {1})>".format(self.stream, self.name)


def realpath(*args):
    """Join all args and return the real path, rooted at ``/``.

    Expands '~', '~user', and environment variables such as :envvar:`$HOME`.

    Returns ``None`` if any of the args is ``None``.
    """
    if None in args:
        return None
    return os.path.realpath(os.path.expanduser(os.path.expandvars(os.path.join(*args))))


def get_ext(filename):
    """Return the lower-cased extension of `filename` without a leading dot.

    Parameters
    ----------
    filename : str

    Returns
    -------
    root : str
    ext : str
    """
    root, ext = os.path.splitext(filename)

    if ext.startswith(os.extsep):
        ext = ext[1:]

    return root, ext.lower()


def check_compressed_format(root, ext):
    """Check if this is a supported gzipped/bzip2ed file format and return UPPERCASE format.

    Parameters
    ----------
    root : str
       path of a file, without extension `ext`
    ext : str
       extension (currently only "bz2" and "gz" are recognized as compressed formats)

    Returns
    -------
    format : str
       upper case format extension *if* the compression can be handled by
       :func:`openany`

    See Also
    --------
    openany : function that is used to open and decompress formats on the fly; only
              compression formats implemented in :func:`openany` are recognized

    """
    # XYZReader&others are setup to handle both plain and compressed (bzip2, gz) files
    # ..so if the first file extension is bzip2 or gz, look at the one to the left of it
    if ext.lower() in ("bz2", "gz"):
        try:
            root, ext = get_ext(root)
        except Exception:
            errmsg = f"Cannot determine coordinate format for '{root}.{ext}'"
            raise TypeError(errmsg) from None

    return ext.upper()


def format_from_filename_extension(filename):
    """Guess file format from the file extension.

    Parameters
    ----------
    filename : str

    Returns
    -------
    format : str

    Raises
    ------
    TypeError
        if the file format cannot be determined
    """
    try:
        root, ext = get_ext(filename)
    except Exception:
        errmsg = (f"Cannot determine file format for file '{filename}'.\n"
                  f"           You can set the format explicitly with "
                  f"'Universe(..., format=FORMAT)'.")
        raise TypeError(errmsg) from None
    format = check_compressed_format(root, ext)
    return format


def guess_format(filename):
    """Return the format of `filename`

    The current heuristic simply looks at the filename extension and can work
    around compressed format extensions.

    Parameters
    ----------
    filename : str or stream
        path to the file or a stream, in which case ``filename.name`` is looked
        at for a hint to the format

    Returns
    -------
    format : str
        format specifier (upper case)

    Raises
    ------
    ValueError
        if the heuristics are insufficient to guess a supported format


    .. versionadded:: 0.11.0
       Moved into lib.util

    """
    if isstream(filename):
        # perhaps StringIO or open stream
        try:
            format = format_from_filename_extension(filename.name)
        except AttributeError:
            # format is None so we need to complain:
            errmsg = (f"guess_format requires an explicit format specifier "
                      f"for stream {filename}")
            raise ValueError(errmsg) from None
    else:
        # iterator, list, filename: simple extension checking... something more
        # complicated is left for the ambitious.
        # Note: at the moment the upper-case extension *is* the format specifier
        # and list of filenames is handled by ChainReader
        format = (format_from_filename_extension(filename)
                  if not iterable(filename) else 'CHAIN')

    return format.upper()


def iterable(obj):
    """Returns ``True`` if `obj` can be iterated over and is *not* a  string
    nor a :class:`NamedStream`"""
    if isinstance(obj, (str, NamedStream)):
        return False  # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True  # any iterator will do
    try:
        len(obj)  # anything else that might work
    except (TypeError, AttributeError):
        return False
    return True


def asiterable(obj):
    """Returns `obj` so that it can be iterated over.

    A string is *not* detected as and iterable and is wrapped into a :class:`list`
    with a single element.

    See Also
    --------
    iterable

    """
    if not iterable(obj):
        obj = [obj]
    return obj


#: Regular expresssion (see :mod:`re`) to parse a simple `FORTRAN edit descriptor`_.
#: ``(?P<repeat>\d?)(?P<format>[IFELAX])(?P<numfmt>(?P<length>\d+)(\.(?P<decimals>\d+))?)?``
#:
#: .. _FORTRAN edit descriptor: http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
FORTRAN_format_regex = (r"(?P<repeat>\d+?)(?P<format>[IFEAX])"
                        r"(?P<numfmt>(?P<length>\d+)(\.(?P<decimals>\d+))?)?")
_FORTRAN_format_pattern = re.compile(FORTRAN_format_regex)


def strip(s):
    """Convert `s` to a string and return it white-space stripped."""
    return str(s).strip()


class FixedcolumnEntry(object):
    """Represent an entry at specific fixed columns.

    Reads from line[start:stop] and converts according to
    typespecifier.
    """
    convertors = {'I': int, 'F': float, 'E': float, 'A': strip}

    def __init__(self, start, stop, typespecifier):
        """
        Parameters
        ----------
        start : int
            first column
        stop : int
            last column + 1
        typespecifier : str
            'I': int, 'F': float, 'E': float, 'A': stripped string

        The start/stop arguments follow standard Python convention in that
        they are 0-based and that the *stop* argument is not included.
        """
        self.start = start
        self.stop = stop
        self.typespecifier = typespecifier
        self.convertor = self.convertors[typespecifier]

    def read(self, line):
        """Read the entry from `line` and convert to appropriate type."""
        try:
            return self.convertor(line[self.start:self.stop])
        except ValueError:
            errmsg = (f"{self}: Failed to read&convert "
                      f"{line[self.start:self.stop]}")
            raise ValueError(errmsg) from None

    def __len__(self):
        """Length of the field in columns (stop - start)"""
        return self.stop - self.start

    def __repr__(self):
        return "FixedcolumnEntry({0:d},{1:d},{2!r})".format(self.start, self.stop, self.typespecifier)


class FORTRANReader(object):
    """FORTRANReader provides a method to parse FORTRAN formatted lines in a file.

    The contents of lines in a file can be parsed according to FORTRAN format
    edit descriptors (see `Fortran Formats`_ for the syntax).

    Only simple one-character specifiers supported here: *I F E A X* (see
    :data:`FORTRAN_format_regex`).

    Strings are stripped of leading and trailing white space.

    .. _`Fortran Formats`: http://www.webcitation.org/5xbaWMV2x
    .. _`Fortran Formats (URL)`:
       http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html

    """

    def __init__(self, fmt):
        """Set up the reader with the FORTRAN format string.

        The string `fmt` should look like '2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10'.

        Parameters
        ----------
        fmt : str
           FORTRAN format edit descriptor for a line as described in `Fortran
           Formats`_

        Example
        -------
        Parsing of a standard CRD file::

           atomformat = FORTRANReader('2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10')
           for line in open('coordinates.crd'):
               serial,TotRes,resName,name,x,y,z,chainID,resSeq,tempFactor = atomformat.read(line)

        """
        self.fmt = fmt.split(',')
        descriptors = [self.parse_FORTRAN_format(descriptor) for descriptor in self.fmt]
        start = 0
        self.entries = []
        for d in descriptors:
            if d['format'] != 'X':
                for x in range(d['repeat']):
                    stop = start + d['length']
                    self.entries.append(FixedcolumnEntry(start, stop, d['format']))
                    start = stop
            else:
                start += d['totallength']

    def read(self, line):
        """Parse `line` according to the format string and return list of values.

        Values are converted to Python types according to the format specifier.

        Parameters
        ----------
        line : str

        Returns
        -------
        list
            list of entries with appropriate types

        Raises
        ------
        ValueError
            Any of the conversions cannot be made (e.g. space for an int)

        See Also
        --------
        :meth:`FORTRANReader.number_of_matches`
        """
        return [e.read(line) for e in self.entries]

    def number_of_matches(self, line):
        """Return how many format entries could be populated with legal values."""
        # not optimal, I suppose...
        matches = 0
        for e in self.entries:
            try:
                e.read(line)
                matches += 1
            except ValueError:
                pass
        return matches

    def parse_FORTRAN_format(self, edit_descriptor):
        """Parse the descriptor.


        Parameters
        ----------
        edit_descriptor : str
            FORTRAN format edit descriptor

        Returns
        -------
        dict
            dict with totallength (in chars), repeat, length, format, decimals

        Raises
        ------
        ValueError
            The `edit_descriptor` is not recognized and cannot be parsed

        Note
        ----
        Specifiers: *L ES EN T TL TR / r S SP SS BN BZ* are *not* supported,
        and neither are the scientific notation *Ew.dEe* forms.

        """

        m = _FORTRAN_format_pattern.match(edit_descriptor.upper())
        if m is None:
            try:
                m = _FORTRAN_format_pattern.match("1" + edit_descriptor.upper())
                if m is None:
                    raise ValueError  # really no idea what the descriptor is supposed to mean
            except:
                raise ValueError("unrecognized FORTRAN format {0!r}".format(edit_descriptor))
        d = m.groupdict()
        if d['repeat'] == '':
            d['repeat'] = 1
        if d['format'] == 'X':
            d['length'] = 1
        for k in ('repeat', 'length', 'decimals'):
            try:
                d[k] = int(d[k])
            except ValueError:  # catches ''
                d[k] = 0
            except TypeError:  # keep None
                pass
        d['totallength'] = d['repeat'] * d['length']
        return d

    def __len__(self):
        """Returns number of entries."""
        return len(self.entries)

    def __repr__(self):
        return self.__class__.__name__ + "(" + ",".join(self.fmt) + ")"


def fixedwidth_bins(delta, xmin, xmax):
    """Return bins of width `delta` that cover `xmin`, `xmax` (or a larger range).

    The bin parameters are computed such that the bin size `delta` is
    guaranteed. In order to achieve this, the range `[xmin, xmax]` can be
    increased.

    Bins can be calculated for 1D data (then all parameters are simple floats)
    or nD data (then parameters are supplied as arrays, with each entry
    correpsonding to one dimension).

    Parameters
    ----------
    delta : float or array_like
        desired spacing of the bins
    xmin : float or array_like
        lower bound (left boundary of first bin)
    xmax : float or array_like
        upper bound (right boundary of last bin)

    Returns
    -------
    dict
        The dict contains 'Nbins', 'delta', 'min', and 'max'; these are either
        floats or arrays, depending on the input.

    Example
    -------
    Use with :func:`numpy.histogram`::

       B = fixedwidth_bins(delta, xmin, xmax)
       h, e = np.histogram(data, bins=B['Nbins'], range=(B['min'], B['max']))

    """
    if not np.all(xmin < xmax):
        raise ValueError('Boundaries are not sane: should be xmin < xmax.')
    _delta = np.asarray(delta, dtype=np.float64)
    _xmin = np.asarray(xmin, dtype=np.float64)
    _xmax = np.asarray(xmax, dtype=np.float64)
    _length = _xmax - _xmin
    N = np.ceil(_length / _delta).astype(np.int_)  # number of bins
    dx = 0.5 * (N * _delta - _length)  # add half of the excess to each end
    return {'Nbins': N, 'delta': _delta, 'min': _xmin - dx, 'max': _xmax + dx}

def get_weights(atoms, weights):
    """Check that a `weights` argument is compatible with `atoms`.

    Parameters
    ----------
    atoms : AtomGroup or array_like
        The atoms that the `weights` should be applied to. Typically this
        is a :class:`AtomGroup` but because only the length is compared,
        any sequence for which ``len(atoms)`` is defined is acceptable.
    weights : {"mass", None} or array_like
        All MDAnalysis functions or classes understand "mass" and will then
        use ``atoms.masses``. ``None`` indicates equal weights for all atoms.
        Using an ``array_like`` assigns a custom weight to each element of
        `atoms`.

    Returns
    -------
    weights : array_like or None
         If "mass" was selected, ``atoms.masses`` is returned, otherwise the
         value of `weights` (which can be ``None``).

    Raises
    ------
    TypeError
        If `weights` is not one of the allowed values or if "mass" is
        selected but ``atoms.masses`` is not available.
    ValueError
        If `weights` is not a 1D array with the same length as
        `atoms`, then the exception is raised.  :exc:`TypeError` is
        also raised if ``atoms.masses`` is not defined.

    """
    if not iterable(weights) and weights == "mass":
        try:
            weights = atoms.masses
        except AttributeError:
            errmsg = "weights='mass' selected but atoms.masses is missing"
            raise TypeError(errmsg) from None

    if iterable(weights):
        if len(np.asarray(weights, dtype=object).shape) != 1:
            raise ValueError("weights must be a 1D array, not with shape "
                            "{0}".format(np.asarray(weights,
                             dtype=object).shape))
        elif len(weights) != len(atoms):
            raise ValueError("weights (length {0}) must be of same length as "
                             "the atoms ({1})".format(
                                 len(weights), len(atoms)))
    elif weights is not None:
        raise ValueError("weights must be {'mass', None} or an iterable of the "
                         "same size as the atomgroup.")

    return weights


# String functions
# ----------------

#: translation table for 3-letter codes --> 1-letter codes
#: .. SeeAlso:: :data:`alternative_inverse_aa_codes`
canonical_inverse_aa_codes = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
#: translation table for 1-letter codes --> *canonical* 3-letter codes.
#: The table is used for :func:`convert_aa_code`.
amino_acid_codes = {one: three for three,
                    one in canonical_inverse_aa_codes.items()}
#: non-default charge state amino acids or special charge state descriptions
#: (Not fully synchronized with :class:`MDAnalysis.core.selection.ProteinSelection`.)
alternative_inverse_aa_codes = {
    'HISA': 'H', 'HISB': 'H', 'HSE': 'H', 'HSD': 'H', 'HID': 'H', 'HIE': 'H', 'HIS1': 'H',
    'HIS2': 'H',
    'ASPH': 'D', 'ASH': 'D',
    'GLUH': 'E', 'GLH': 'E',
    'LYSH': 'K', 'LYN': 'K',
    'ARGN': 'R',
    'CYSH': 'C', 'CYS1': 'C', 'CYS2': 'C'}
#: lookup table from 3/4 letter resnames to 1-letter codes. Note that non-standard residue names
#: for tautomers or different protonation states such as HSE are converted to canonical 1-letter codes ("H").
#: The table is used for :func:`convert_aa_code`.
#: .. SeeAlso:: :data:`canonical_inverse_aa_codes` and :data:`alternative_inverse_aa_codes`
inverse_aa_codes = {}
inverse_aa_codes.update(canonical_inverse_aa_codes)
inverse_aa_codes.update(alternative_inverse_aa_codes)


def convert_aa_code(x):
    """Converts between 3-letter and 1-letter amino acid codes.

    Parameters
    ----------
    x : str
        1-letter or 3-letter amino acid code

    Returns
    -------
    str
        3-letter or 1-letter amino acid code

    Raises
    ------
    ValueError
        No conversion can be made; the amino acid code is not defined.

    Note
    ----
    Data are defined in :data:`amino_acid_codes` and :data:`inverse_aa_codes`.
    """
    if len(x) == 1:
        d = amino_acid_codes
    else:
        d = inverse_aa_codes

    try:
        return d[x.upper()]
    except KeyError:
        errmsg = (f"No conversion for {x} found (1 letter -> 3 letter or 3/4 "
                  f"letter -> 1 letter)")
        raise ValueError(errmsg) from None


#: Regular expression to match and parse a residue-atom selection; will match
#: "LYS300:HZ1" or "K300:HZ1" or "K300" or "4GB300:H6O" or "4GB300" or "YaA300".
RESIDUE = re.compile(r"""
                 (?P<aa>([ACDEFGHIKLMNPQRSTVWY])   # 1-letter amino acid
                        |                          #   or
                        ([0-9A-Z][a-zA-Z][A-Z][A-Z]?)    # 3-letter or 4-letter residue name
                 )
                 \s*                               # white space allowed
                 (?P<resid>\d+)                    # resid
                 \s*
                 (:                                # separator ':'
                   \s*
                   (?P<atom>\w+)                   # atom name
                 )?                                # possibly one
            """, re.VERBOSE | re.IGNORECASE)


# from GromacsWrapper cbook.IndexBuilder
def parse_residue(residue):
    """Process residue string.

    Parameters
    ----------
    residue: str
        The *residue* must contain a 1-letter or 3-letter or
        4-letter residue string, a number (the resid) and
        optionally an atom identifier, which must be separate
        from the residue with a colon (":"). White space is
        allowed in between.

    Returns
    -------
    tuple
        `(3-letter aa string, resid, atomname)`; known 1-letter
        aa codes are converted to 3-letter codes

    Examples
    --------
     - "LYS300:HZ1" --> ("LYS", 300, "HZ1")
     - "K300:HZ1" --> ("LYS", 300, "HZ1")
     - "K300" --> ("LYS", 300, None)
     - "4GB300:H6O" --> ("4GB", 300, "H6O")
     - "4GB300" --> ("4GB", 300, None)

    """

    # XXX: use _translate_residue() ....
    m = RESIDUE.match(residue)
    if not m:
        raise ValueError("Selection {residue!r} is not valid (only 1/3/4 letter resnames, resid required).".format(**vars()))
    resid = int(m.group('resid'))
    residue = m.group('aa')
    if len(residue) == 1:
        resname = convert_aa_code(residue)  # only works for AA
    else:
        resname = residue  # use 3-letter for any resname
    atomname = m.group('atom')
    return (resname, resid, atomname)


def conv_float(s):
    """Convert an object `s` to float if possible.

    Function to be passed into :func:`map` or a list comprehension. If
    the argument can be interpreted as a float it is converted,
    otherwise the original object is passed back.
    """
    try:
        return float(s)
    except ValueError:
        return s


# A dummy, empty, cheaply-hashable object class to use with weakref caching.
# (class object doesn't allow weakrefs to its instances, but user-defined
#  classes do)
class _CacheKey:
    pass


def cached(key, universe_validation=False):
    """Cache a property within a class.

    Requires the Class to have a cache dict :attr:`_cache` and, with
    `universe_validation`, a :attr:`universe` with a cache dict :attr:`_cache`.

    Example
    -------
    How to add a cache for a variable to a class by using the `@cached`
    decorator::

       class A(object):
           def__init__(self):
               self._cache = dict()

           @property
           @cached('keyname')
           def size(self):
               # This code gets run only if the lookup of keyname fails
               # After this code has been run once, the result is stored in
               # _cache with the key: 'keyname'
               return 10.0

           @property
           @cached('keyname', universe_validation=True)
           def othersize(self):
               # This code gets run only if the lookup
               # id(self) is not in the validation set under
               # self.universe._cache['_valid']['keyname']
               # After this code has been run once, id(self) is added to that
               # set. The validation set can be centrally invalidated at the
               # universe level (say, if a topology change invalidates specific
               # caches).
               return 20.0


    .. versionadded:: 0.9.0

    .. versionchanged::2.0.0
        Added the `universe_validation` keyword.
    """

    def cached_lookup(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            try:
                if universe_validation:  # Universe-level cache validation
                    u_cache = self.universe._cache.setdefault('_valid', dict())
                    # A WeakSet is used so that keys from out-of-scope/deleted
                    # objects don't clutter it.
                    valid_caches = u_cache.setdefault(key, weakref.WeakSet())
                    try:
                        if self._cache_key not in valid_caches:
                            raise KeyError
                    except AttributeError:  # No _cache_key yet
                        # Must create a reference key for the validation set.
                        # self could be used itself as a weakref but set()
                        # requires hashing it, which can be slow for AGs. Using
                        # id(self) fails because ints can't be weak-referenced.
                        self._cache_key = _CacheKey()
                        raise KeyError
                return self._cache[key]
            except KeyError:
                self._cache[key] = ret = func(self, *args, **kwargs)
                if universe_validation:
                    valid_caches.add(self._cache_key)
                return ret

        return wrapper

    return cached_lookup


def unique_rows(arr, return_index=False):
    """Return the unique rows of an array.

    Arguments
    ---------
    arr : numpy.ndarray
        Array of shape ``(n1, m)``.
    return_index : bool, optional
      If ``True``, returns indices of array that formed answer (see
      :func:`numpy.unique`)

    Returns
    -------
    unique_rows : numpy.ndarray
         Array of shape ``(n2, m)`` containing only the unique rows of `arr`.
    r_idx : numpy.ndarray (optional)
          Array containing the corresponding row indices (if `return_index`
          is ``True``).

    Examples
    --------
    Remove dupicate rows from an array:

    >>> import numpy as np
    >>> from MDAnalysis.lib.util import unique_rows
    >>> a = np.array([[0, 1], [1, 2], [1, 2], [0, 1], [2, 3]])
    >>> b = unique_rows(a)
    >>> b
    array([[0, 1],
           [1, 2],
           [2, 3]])

    See Also
    --------
    numpy.unique

    """
    # From here, but adapted to handle any size rows
    # https://mail.scipy.org/pipermail/scipy-user/2011-December/031200.html

    # This seems to fail if arr.flags['OWNDATA'] is False
    # this can occur when second dimension was created through broadcasting
    # eg: idx = np.array([1, 2])[None, :]
    if not arr.flags['OWNDATA']:
        arr = arr.copy()

    m = arr.shape[1]

    if return_index:
        u, r_idx = np.unique(arr.view(dtype=np.dtype([(str(i), arr.dtype)
                                                      for i in range(m)])),
                             return_index=True)
        return u.view(arr.dtype).reshape(-1, m), r_idx
    else:
        u = np.unique(arr.view(
            dtype=np.dtype([(str(i), arr.dtype) for i in range(m)])
        ))
        return u.view(arr.dtype).reshape(-1, m)


def blocks_of(a, n, m):
    """Extract a view of ``(n, m)`` blocks along the diagonal of the array `a`.

    Parameters
    ----------
    a : numpy.ndarray
        Input array, must be C contiguous and at least 2D.
    n : int
        Size of block in first dimension.
    m : int
        Size of block in second dimension.

    Returns
    -------
    view : numpy.ndarray
        A view of the original array with shape ``(nblocks, n, m)``, where
        ``nblocks`` is the number of times the miniblocks of shape ``(n, m)``
        fit in the original.

    Raises
    ------
    ValueError
        If the supplied `n` and `m` don't divide `a` into an integer number
        of blocks or if `a` is not C contiguous.

    Examples
    --------
    >>> import numpy as np
    >>> from MDAnalysis.lib.util import blocks_of
    >>> arr = np.arange(16).reshape(4, 4)
    >>> view = blocks_of(arr, 2, 2)
    >>> view[:] = 100
    >>> arr
    array([[100, 100,   2,   3],
           [100, 100,   6,   7],
           [  8,   9, 100, 100],
           [ 12,  13, 100, 100]])

    Notes
    -----
    `n`, `m` must divide `a` into an identical integer number of blocks. Please
    note that if the block size is larger than the input array, this number will
    be zero, resulting in an empty view!

    Uses strides and therefore requires that the array is C contiguous.

    Returns a view, so editing this modifies the original array.


    .. versionadded:: 0.12.0

    """
    # based on:
    # http://stackoverflow.com/a/10862636
    # but generalised to handle non square blocks.
    if not a.flags['C_CONTIGUOUS']:
        raise ValueError("Input array is not C contiguous.")

    nblocks = a.shape[0] // n
    nblocks2 = a.shape[1] // m

    if not nblocks == nblocks2:
        raise ValueError("Must divide into same number of blocks in both"
                         " directions.  Got {} by {}"
                         "".format(nblocks, nblocks2))

    new_shape = (nblocks, n, m)
    new_strides = (n * a.strides[0] + m * a.strides[1],
                   a.strides[0], a.strides[1])

    return np.lib.stride_tricks.as_strided(a, new_shape, new_strides)


def group_same_or_consecutive_integers(arr):
    """Split an array of integers into a list of same or consecutive
    sequences.

    Parameters
    ----------
    arr: :class:`numpy.ndarray`

    Returns
    -------
    list of :class:`numpy.ndarray`

    Examples
    >>> import numpy as np
    >>> arr = np.array([ 2,  3,  4,  7,  8,  9, 10, 11, 15, 16])
    >>> group_same_or_consecutive_integers(arr)
    [array([2, 3, 4]), array([ 7,  8,  9, 10, 11]), array([15, 16])]
    """
    return np.split(arr, np.where(np.ediff1d(arr)-1 > 0)[0] + 1)


class Namespace(dict):
    """Class to allow storing attributes in new namespace. """

    def __getattr__(self, key):
        # a.this causes a __getattr__ call for key = 'this'
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            errmsg = f'"{key}" is not known in the namespace.'
            raise AttributeError(errmsg) from None

    def __setattr__(self, key, value):
        dict.__setitem__(self, key, value)

    def __delattr__(self, key):
        try:
            dict.__delitem__(self, key)
        except KeyError:
            errmsg = f'"{key}" is not known in the namespace.'
            raise AttributeError(errmsg) from None

    def __eq__(self, other):
        try:
            # this'll allow us to compare if we're storing arrays
            assert_equal(self, other)
        except AssertionError:
            return False
        return True


def ltruncate_int(value, ndigits):
    """Truncate an integer, retaining least significant digits

    Parameters
    ----------
    value : int
      value to truncate
    ndigits : int
      number of digits to keep

    Returns
    -------
    truncated : int
      only the `ndigits` least significant digits from `value`

    Examples
    --------
    >>> from MDAnalysis.lib.util import ltruncate_int
    >>> ltruncate_int(123, 2)
    23
    >>> ltruncate_int(1234, 5)
    1234
    """
    return int(str(value)[-ndigits:])


def flatten_dict(d, parent_key=tuple()):
    """Flatten a nested dict `d` into a shallow dict with tuples as keys.

    Parameters
    ----------
    d : dict

    Returns
    -------
    dict

    Note
    -----
    Based on https://stackoverflow.com/a/6027615/
    by user https://stackoverflow.com/users/1897/imran

    .. versionadded:: 0.18.0
    """

    items = []
    for k, v in d.items():
        if type(k) != tuple:
            new_key = parent_key + (k, )
        else:
            new_key = parent_key + k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key).items())
        else:
            items.append((new_key, v))
    return dict(items)


def static_variables(**kwargs):
    """Decorator equipping functions or methods with static variables.

    Static variables are declared and initialized by supplying keyword arguments
    and initial values to the decorator.

    Example
    -------

    >>> from MDAnalysis.lib.util import static_variables
    >>> @static_variables(msg='foo calls', calls=0)
    ... def foo():
    ...     foo.calls += 1
    ...     print("{}: {}".format(foo.msg, foo.calls))
    ...
    >>> foo()
    foo calls: 1
    >>> foo()
    foo calls: 2


    .. note:: Based on https://stackoverflow.com/a/279586
        by `Claudiu <https://stackoverflow.com/users/15055/claudiu>`_

    .. versionadded:: 0.19.0
    """
    def static_decorator(func):
        for kwarg in kwargs:
            setattr(func, kwarg, kwargs[kwarg])
        return func
    return static_decorator


# In a lot of Atom/Residue/SegmentGroup methods such as center_of_geometry() and
# the like, results are biased if the calling group is not unique, i.e., if it
# contains duplicates.
# We therefore raise a `DuplicateWarning` whenever an affected method is called
# from a non-unique group. Since several of the affected methods involve calls
# to other affected methods, simply raising a warning in every affected method
# would potentially lead to a massive amount of warnings. This is exactly where
# the `warn_if_unique` decorator below comes into play. It ensures that a
# warning is only raised once for a method using this decorator, and suppresses
# all such warnings that would potentially be raised in methods called by that
# method. Of course, as it is generally the case with Python warnings, this is
# *not threadsafe*.

@static_variables(warned=False)
def warn_if_not_unique(groupmethod):
    """Decorator triggering a :class:`~MDAnalysis.exceptions.DuplicateWarning`
    if the underlying group is not unique.

    Assures that during execution of the decorated method only the first of
    potentially multiple warnings concerning the uniqueness of groups is shown.

    Raises
    ------
    :class:`~MDAnalysis.exceptions.DuplicateWarning`
        If the :class:`~MDAnalysis.core.groups.AtomGroup`,
        :class:`~MDAnalysis.core.groups.ResidueGroup`, or
        :class:`~MDAnalysis.core.groups.SegmentGroup` of which the decorated
        method is a member contains duplicates.


    .. versionadded:: 0.19.0
    """
    @wraps(groupmethod)
    def wrapper(group, *args, **kwargs):
        # Proceed as usual if the calling group is unique or a DuplicateWarning
        # has already been thrown:
        if group.isunique or warn_if_not_unique.warned:
            return groupmethod(group, *args, **kwargs)
        # Otherwise, throw a DuplicateWarning and execute the method.
        method_name = ".".join(
            (group.__class__.__name__, groupmethod.__name__))
        # Try to get the group's variable name(s):
        caller_locals = inspect.currentframe().f_back.f_locals.items()
        group_names = []
        for name, obj in caller_locals:
            try:
                if obj is group:
                    group_names.append("'{}'".format(name))
            except:
                pass
        if not group_names:
            group_name = "'unnamed {}'".format(group.__class__.__name__)
        elif len(group_names) == 1:
            group_name = group_names[0]
        else:
            group_name = " a.k.a. ".join(sorted(group_names))
        group_repr = repr(group)
        msg = ("{}(): {} {} contains duplicates. Results might be biased!"
               "".format(method_name, group_name, group_repr))
        warnings.warn(message=msg, category=DuplicateWarning, stacklevel=2)
        warn_if_not_unique.warned = True
        try:
            result = groupmethod(group, *args, **kwargs)
        finally:
            warn_if_not_unique.warned = False
        return result
    return wrapper


def check_coords(*coord_names, **options):
    """Decorator for automated coordinate array checking.

    This decorator is intended for use especially in
    :mod:`MDAnalysis.lib.distances`.
    It takes an arbitrary number of positional arguments which must correspond
    to names of positional arguments of the decorated function.
    It then checks if the corresponding values are valid coordinate arrays or
    an :class:`~MDAnalysis.core.groups.AtomGroup`.
    If the input is an array and all these arrays are single coordinates
    (i.e., their shape is ``(3,)``), the decorated function can optionally
    return a single coordinate (or angle) instead of an array of coordinates
    (or angles). This can be used to enable computations of single observables
    using functions originally designed to accept only 2-d coordinate arrays.

    If the input is an :class:`~MDAnalysis.core.groups.AtomGroup` it is
    converted into its corresponding position array via a call to
    `AtomGroup.positions`.

    The checks performed on each individual coordinate array are:

    * Check that coordinate arrays are of type :class:`numpy.ndarray`.
    * Check that coordinate arrays have a shape of ``(n, 3)`` (or ``(3,)`` if
      single coordinates are allowed; see keyword argument `allow_single`).
    * Automatic dtype conversion to ``numpy.float32``.
    * Optional replacement by a copy; see keyword argument `enforce_copy` .
    * If coordinate arrays aren't C-contiguous, they will be automatically
      replaced by a C-contiguous copy.
    * Optional check for equal length of all coordinate arrays; see optional
      keyword argument `check_lengths_match`.

    Parameters
    ----------
    *coord_names : tuple
        Arbitrary number of strings corresponding to names of positional
        arguments of the decorated function.
    **options : dict, optional
        * **enforce_copy** (:class:`bool`, optional) -- Enforce working on a
          copy of the coordinate arrays. This is useful to ensure that the input
          arrays are left unchanged. Default: ``True``
        * **enforce_dtype** (:class:`bool`, optional) -- Enforce a conversion
          to float32.  Default: ``True``
        * **allow_single** (:class:`bool`, optional) -- Allow the input
          coordinate array to be a single coordinate with shape ``(3,)``.
        * **convert_single** (:class:`bool`, optional) -- If ``True``, single
          coordinate arrays will be converted to have a shape of ``(1, 3)``.
          Only has an effect if `allow_single` is ``True``. Default: ``True``
        * **reduce_result_if_single** (:class:`bool`, optional) -- If ``True``
          and *all* input coordinates are single, a decorated function ``func``
          will return ``func()[0]`` instead of ``func()``. Only has an effect if
          `allow_single` is ``True``. Default: ``True``
        * **check_lengths_match** (:class:`bool`, optional) -- If ``True``, a
          :class:`ValueError` is raised if not all coordinate arrays contain the
          same number of coordinates. Default: ``True``
        * **allow_atomgroup** (:class:`bool`, optional) -- If ``False``, a
          :class:`TypeError` is raised if an :class:`AtomGroup` is supplied
          Default: ``False``

    Raises
    ------
    ValueError
        If the decorator is used without positional arguments (for development
        purposes only).

        If any of the positional arguments supplied to the decorator doesn't
        correspond to a name of any of the decorated function's positional
        arguments.

        If any of the coordinate arrays has a wrong shape.
    TypeError
        If any of the coordinate arrays is not a :class:`numpy.ndarray` or an
        :class:`~MDAnalysis.core.groups.AtomGroup`.

        If the dtype of any of the coordinate arrays is not convertible to
          ``numpy.float32``.

    Example
    -------

    >>> import numpy as np
    >>> import MDAnalysis as mda
    >>> from MDAnalysis.tests.datafiles import PSF, DCD
    >>> from MDAnalysis.lib.util import check_coords
    >>> @check_coords('coords1', 'coords2', allow_atomgroup=True)
    ... def coordsum(coords1, coords2):
    ...     assert coords1.dtype == np.float32
    ...     assert coords2.flags['C_CONTIGUOUS']
    ...     return coords1 + coords2
    ...
    >>> # automatic dtype conversion:
    >>> coordsum(np.zeros(3, dtype=np.int64), np.ones(3))
    array([1., 1., 1.], dtype=float32)
    >>>
    >>> # automatic handling of non-contiguous arrays:
    >>> coordsum(np.zeros(3), np.ones(6)[::2])
    array([1., 1., 1.], dtype=float32)
    >>>
    >>> # automatic handling of AtomGroups
    >>> u = mda.Universe(PSF, DCD)
    >>> try:
    ...     coordsum(u.atoms, u.select_atoms("index 1 to 10"))
    ... except ValueError as err:
    ...     err
    ValueError('coordsum(): coords1, coords2 must contain the same number of coordinates, got [3341, 10].')
    >>>
    >>> # automatic shape checking:
    >>> try:
    ...     coordsum(np.zeros(3), np.ones(6))
    ... except ValueError as err:
    ...     err
    ValueError('coordsum(): coords2.shape must be (3,) or (n, 3), got (6,)')


    .. versionadded:: 0.19.0
    .. versionchanged:: 2.3.0
       Can now accept an :class:`AtomGroup` as input, and added option
       allow_atomgroup with default False to retain old behaviour
    """
    enforce_copy = options.get('enforce_copy', True)
    enforce_dtype = options.get('enforce_dtype', True)
    allow_single = options.get('allow_single', True)
    convert_single = options.get('convert_single', True)
    reduce_result_if_single = options.get('reduce_result_if_single', True)
    check_lengths_match = options.get('check_lengths_match',
                                      len(coord_names) > 1)
    allow_atomgroup = options.get('allow_atomgroup', False)
    if not coord_names:
        raise ValueError("Decorator check_coords() cannot be used without "
                         "positional arguments.")

    def check_coords_decorator(func):
        fname = func.__name__
        code = func.__code__
        argnames = code.co_varnames
        nargs = len(code.co_varnames)
        ndefaults = len(func.__defaults__) if func.__defaults__ else 0
        # Create a tuple of positional argument names:
        nposargs = code.co_argcount - ndefaults
        posargnames = argnames[:nposargs]
        # The check_coords() decorator is designed to work only for positional
        # arguments:
        for name in coord_names:
            if name not in posargnames:
                raise ValueError("In decorator check_coords(): Name '{}' "
                                 "doesn't correspond to any positional "
                                 "argument of the decorated function {}()."
                                 "".format(name, func.__name__))

        def _check_coords(coords, argname):
            is_single = False
            if isinstance(coords, np.ndarray):
                if allow_single:
                    if (coords.ndim not in (1, 2)) or (coords.shape[-1] != 3):
                        errmsg = (f"{fname}(): {argname}.shape must be (3,) or "
                                  f"(n, 3), got {coords.shape}")
                        raise ValueError(errmsg)
                    if coords.ndim == 1:
                        is_single = True
                        if convert_single:
                            coords = coords[None, :]
                else:
                    if (coords.ndim != 2) or (coords.shape[1] != 3):
                        errmsg = (f"{fname}(): {argname}.shape must be (n, 3) "
                                  f"got {coords.shape}")
                        raise ValueError(errmsg)
                if enforce_dtype:
                    try:
                        coords = coords.astype(
                            np.float32, order='C', copy=enforce_copy)
                    except ValueError:
                        errmsg = (f"{fname}(): {argname}.dtype must be"
                                  f"convertible to float32, got"
                                  f" {coords.dtype}.")
                        raise TypeError(errmsg) from None
                # coordinates should now be the right shape
                ncoord = coords.shape[0]
            else:
                try:
                    coords = coords.positions  # homogenise to a numpy array
                    ncoord = coords.shape[0]
                    if not allow_atomgroup:
                        err = TypeError("AtomGroup or other class with a"
                                        "`.positions` method supplied as an"
                                        "argument, but allow_atomgroup is"
                                        " False")
                        raise err
                except AttributeError:
                    raise TypeError(f"{fname}(): Parameter '{argname}' must be"
                                    f" a numpy.ndarray or an AtomGroup,"
                                    f" got {type(coords)}.")

            return coords, is_single, ncoord

        @wraps(func)
        def wrapper(*args, **kwargs):
            # Check for invalid function call:
            if len(args) != nposargs:
                # set marker for testing purposes:
                wrapper._invalid_call = True
                if len(args) > nargs:
                    # too many arguments, invoke call:
                    return func(*args, **kwargs)
                for name in posargnames[:len(args)]:
                    if name in kwargs:
                        # duplicate argument, invoke call:
                        return func(*args, **kwargs)
                for name in posargnames[len(args):]:
                    if name not in kwargs:
                        # missing argument, invoke call:
                        return func(*args, **kwargs)
                for name in kwargs:
                    if name not in argnames:
                        # unexpected kwarg, invoke call:
                        return func(*args, **kwargs)
                # call is valid, unset test marker:
                wrapper._invalid_call = False
            args = list(args)
            ncoords = []
            all_single = allow_single
            for name in coord_names:
                idx = posargnames.index(name)
                if idx < len(args):
                    args[idx], is_single, ncoord = _check_coords(args[idx],
                                                                 name)
                    all_single &= is_single
                    ncoords.append(ncoord)
                else:
                    kwargs[name], is_single, ncoord = _check_coords(kwargs[name],
                                                                    name)
                    all_single &= is_single
                    ncoords.append(ncoord)
            if check_lengths_match and ncoords:
                if ncoords.count(ncoords[0]) != len(ncoords):
                    raise ValueError("{}(): {} must contain the same number of "
                                     "coordinates, got {}."
                                     "".format(fname, ", ".join(coord_names),
                                               ncoords))
            # If all input coordinate arrays were 1-d, so should be the output:
            if all_single and reduce_result_if_single:
                return func(*args, **kwargs)[0]
            return func(*args, **kwargs)
        return wrapper
    return check_coords_decorator


def check_atomgroup_not_empty(groupmethod):
    """Decorator triggering a ``ValueError`` if the underlying group is empty.

    Avoids downstream errors in computing properties of empty atomgroups. 

    Raises
    ------
    ValueError
        If the input :class:`~MDAnalysis.core.groups.AtomGroup`,
        of a decorated method is empty.


    .. versionadded:: 2.4.0
    """
    @wraps(groupmethod)
    def wrapper(group, *args, **kwargs):
        # Throw error if the group is empty.
        if not group.atoms:
            raise ValueError("AtomGroup is empty.")
        # Proceed as usual if the calling group is not empty.
        else:
            result = groupmethod(group, *args, **kwargs)
        return result
    return wrapper


# ------------------------------------------------------------------
#
# our own deprecate function, derived from numpy (see
# https://github.com/MDAnalysis/mdanalysis/pull/1763#issuecomment-403231136)
#
# From numpy/lib/utils.py 1.14.5 (used under the BSD 3-clause licence,
# https://www.numpy.org/license.html#license) and modified

def _set_function_name(func, name):
    func.__name__ = name
    return func


class _Deprecate(object):
    """
    Decorator class to deprecate old functions.

    Refer to `deprecate` for details.

    See Also
    --------
    deprecate


    .. versionadded:: 0.19.0
    """

    def __init__(self, old_name=None, new_name=None,
                 release=None, remove=None, message=None):
        self.old_name = old_name
        self.new_name = new_name
        if release is None:
            raise ValueError("deprecate: provide release in which "
                             "feature was deprecated.")
        self.release = str(release)
        self.remove = str(remove) if remove is not None else remove
        self.message = message

    def __call__(self, func, *args, **kwargs):
        """
        Decorator call.  Refer to ``decorate``.

        """
        old_name = self.old_name
        new_name = self.new_name
        message = self.message
        release = self.release
        remove = self.remove

        if old_name is None:
            try:
                old_name = func.__name__
            except AttributeError:
                old_name = func.__name__
        if new_name is None:
            depdoc = "`{0}` is deprecated!".format(old_name)
        else:
            depdoc = "`{0}` is deprecated, use `{1}` instead!".format(
                old_name, new_name)

        warn_message = depdoc

        remove_text = ""
        if remove is not None:
            remove_text = "`{0}` will be removed in release {1}.".format(
                old_name, remove)
            warn_message += "\n" + remove_text
        if message is not None:
            warn_message += "\n" + message

        def newfunc(*args, **kwds):
            """This function is deprecated."""
            warnings.warn(warn_message, DeprecationWarning, stacklevel=2)
            return func(*args, **kwds)

        newfunc = _set_function_name(newfunc, old_name)

        # Build the doc string
        # First line: func is deprecated, use newfunc instead!
        # Normal docs follows.
        # Last: .. deprecated::

        # make sure that we do not mess up indentation, otherwise sphinx
        # docs do not build properly
        try:
            doc = dedent_docstring(func.__doc__)
        except TypeError:
            doc = ""

        deprecation_text = dedent_docstring("""\n\n
        .. deprecated:: {0}
           {1}
           {2}
        """.format(release,
                   message if message else depdoc,
                   remove_text))

        doc = "{0}\n\n{1}\n{2}\n".format(depdoc, doc, deprecation_text)

        newfunc.__doc__ = doc
        try:
            d = func.__dict__
        except AttributeError:
            pass
        else:
            newfunc.__dict__.update(d)
        return newfunc


def deprecate(*args, **kwargs):
    r"""Issues a DeprecationWarning, adds warning to `old_name`'s
    docstring, rebinds ``old_name.__name__`` and returns the new
    function object.

    This function may also be used as a decorator.

    It adds a restructured text ``.. deprecated:: release`` block with
    the sphinx deprecated role to the end of the docs. The `message`
    is added under the deprecation block and contains the `release` in
    which the function was deprecated.

    Parameters
    ----------
    func : function
        The function to be deprecated.
    old_name : str, optional
        The name of the function to be deprecated. Default is None, in
        which case the name of `func` is used.
    new_name : str, optional
        The new name for the function. Default is None, in which case the
        deprecation message is that `old_name` is deprecated. If given, the
        deprecation message is that `old_name` is deprecated and `new_name`
        should be used instead.
    release : str
        Release in which the function was deprecated. This is given as
        a keyword argument for technical reasons but is required; a
        :exc:`ValueError` is raised if it is missing.
    remove : str, optional
        Release for which removal of the feature is planned.
    message : str, optional
        Additional explanation of the deprecation.  Displayed in the
        docstring after the warning.

    Returns
    -------
    old_func : function
        The deprecated function.

    Examples
    --------
    When :func:`deprecate` is used as a function as in the following
    example,

    .. code-block:: python

       oldfunc = deprecate(func, release="0.19.0", remove="1.0",
                           message="Do it yourself instead.")

    then ``oldfunc`` will return a value after printing
    :exc:`DeprecationWarning`; ``func`` is still available as it was
    before.

    When used as a decorator, ``func`` will be changed and issue the
    warning and contain the deprecation note in the do string.

    .. code-block:: python

       @deprecate(release="0.19.0", remove="1.0",
                  message="Do it yourself instead.")
       def func():
           \"\"\"Just pass\"\"\"
           pass

    The resulting doc string (``help(func)``) will look like:

    .. code-block:: reST

       `func` is deprecated!

       Just pass.

       .. deprecated:: 0.19.0
          Do it yourself instead.
          `func` will be removed in 1.0.

    (It is possible but confusing to change the name of ``func`` with
    the decorator so it is not recommended to use the `new_func`
    keyword argument with the decorator.)

    .. versionadded:: 0.19.0

    """
    # Deprecate may be run as a function or as a decorator
    # If run as a function, we initialise the decorator class
    # and execute its __call__ method.

    if args:
        fn = args[0]
        args = args[1:]
        return _Deprecate(*args, **kwargs)(fn)
    else:
        return _Deprecate(*args, **kwargs)
#
# ------------------------------------------------------------------


def dedent_docstring(text):
    """Dedent typical python doc string.

    Parameters
    ----------
    text : str
        string, typically something like ``func.__doc__``.

    Returns
    -------
    str
        string with the leading common whitespace removed from each
        line

    See Also
    --------
    textwrap.dedent


    .. versionadded:: 0.19.0
    """
    lines = text.splitlines()
    if len(lines) < 2:
        return text.lstrip()

    # treat first line as special (typically no leading whitespace!) which messes up dedent
    return lines[0].lstrip() + "\n" + textwrap.dedent("\n".join(lines[1:]))


def check_box(box):
    """Take a box input and deduce what type of system it represents based on
    the shape of the array and whether all angles are 90 degrees.

    Parameters
    ----------
    box : array_like
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.

    Returns
    -------
    boxtype : {``'ortho'``, ``'tri_vecs'``}
        String indicating the box type (orthogonal or triclinic).
    checked_box : numpy.ndarray
        Array of dtype ``numpy.float32`` containing box information:
          * If `boxtype` is ``'ortho'``, `cecked_box` will have the shape ``(3,)``
            containing the x-, y-, and z-dimensions of the orthogonal box.
          * If  `boxtype` is ``'tri_vecs'``, `cecked_box` will have the shape
            ``(3, 3)`` containing the triclinic box vectors in a lower triangular
            matrix as returned by
            :meth:`~MDAnalysis.lib.mdamath.triclinic_vectors`.

    Raises
    ------
    ValueError
        If `box` is not of the form ``[lx, ly, lz, alpha, beta, gamma]``
        or contains data that is not convertible to ``numpy.float32``.

    See Also
    --------
    MDAnalysis.lib.mdamath.triclinic_vectors


    .. versionchanged: 0.19.0
       * Enforced correspondence of `box` with specified format.
       * Added automatic conversion of input to :class:`numpy.ndarray` with
         dtype ``numpy.float32``.
       * Now also returns the box in the format expected by low-level functions
         in :mod:`~MDAnalysis.lib.c_distances`.
       * Removed obsolete box types ``tri_box`` and ``tri_vecs_bad``.
    """
    if box is None:
        raise ValueError("Box is None")
    from .mdamath import triclinic_vectors  # avoid circular import
    box = np.asarray(box, dtype=np.float32, order='C')
    if box.shape != (6,):
        raise ValueError("Invalid box information. Must be of the form "
                         "[lx, ly, lz, alpha, beta, gamma].")
    if np.all(box[3:] == 90.):
        return 'ortho', box[:3]
    return 'tri_vecs', triclinic_vectors(box)


def store_init_arguments(func):
    """Decorator to store arguments passed to the init method of a class.

    Arguments are stored as a dictionary in ``cls._kwargs``.

    Notes
    -----
    * Only does a shallow copy, if the arguments are changed
      by the class after passing through the decorator this will be
      reflected in the stored arguments.
    * If not empty, ``args`` is not unpacked and stored as-is in the
      dictionary. If no ``args`` are passed, then no ``arg`` entry will be
      stored in the dictionary.


    .. versionadded:: 2.2.0
    """
    sig = inspect.signature(func)

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if not hasattr(self, "_kwargs"):
            arg_values = sig.bind(self, *args, **kwargs)
            arg_values.apply_defaults()
            self._kwargs = {}
            for key, arg in arg_values.arguments.items():
                if key != "self":
                    if key == "kwargs":
                        for k, v in arg.items():
                            self._kwargs[k] = v
                    elif key == "args":
                        if len(arg) > 0:
                            self._kwargs[key] = arg
                    else:
                        self._kwargs[key] = arg
        return func(self, *args, **kwargs)
    return wrapper


def no_copy_shim():
    if np.lib.NumpyVersion >= "2.0.0rc1":
        copy = None
    else:
        copy = False
    return copy


def atoi(s: str) -> int:
    """Convert the leading number part of a string to an integer.

    Parameters
    ----------
    s : str
        The string to convert to an integer.

    Returns
    -------
    number : int
        The first numeric part of the string converted to an integer.
        If the string does not start with a number, 0 is returned.

    Examples
    --------
    >>> from MDAnalysis.lib.util import atoi
    >>> atoi('34f4')
    34
    >>> atoi('foo')
    0
 

    .. versionadded:: 2.8.0
    """
    try:
        return int(''.join(itertools.takewhile(str.isdigit, s.strip())))
    except ValueError:
        return 0
