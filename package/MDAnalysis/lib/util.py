# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Helper functions --- :mod:`MDAnalysis.lib.util`
====================================================

Small helper functions that don't fit anywhere else.

Files and directories
---------------------

.. autofunction:: filename
.. function:: openany(directory[,mode='r'])

   Context manager to open a compressed (bzip2, gzip) or plain file
   (uses :func:`anyopen`).

.. autofunction:: anyopen
.. autofunction:: greedy_splitext
.. autofunction:: which
.. autofunction:: realpath
.. autofunction:: guess_format

Streams
-------

Many of the readers are not restricted to just reading files. They can
also use gzip-compressed or bzip2-compressed files (through the
internal use of :func:`openany`). It is also possible to provide more
general streams as inputs, such as a :func:`cStringIO.StringIO`
instances (essentially, a memory buffer) by wrapping these instances
into a :class:`NamedStream`. This :class:`NamedStream` can then be
used in place of an ordinary file name (typically, with a
class:`~MDAnalysis.core.AtomGroup.Universe` but it is also possible to
*write* to such a stream using :func:`MDAnalysis.Writer`).

.. rubric: Examples

In the following example, we use a PDB stored as a string ``pdb_s``::

   import MDAnalysis
   from MDAnalysis.lib.util import NamedStream
   import cStringIO

   pdb_s = "TITLE     Lonely Ion\\nATOM      1  NA  NA+     1      81.260  64.982  10.926  1.00  0.00\\n"
   u = MDAnalysis.Universe(NamedStream(cStringIO.StringIO(pdb_s), "ion.pdb"))
   print(u)
   #  <Universe with 1 atoms>
   print(u.atoms.positions)
   # [[ 81.26000214  64.98200226  10.92599964]]

It is important to provide a proper pseudo file name with the correct extension
(".pdb") to :class:`NamedStream` because the file type recognition uses the
extension of the file name to determine the file format or alternatively
provide the ``format="pdb"`` keyword argument to the
:class:`~MDAnalysis.core.AtomGroup.Universe`.

The use of streams becomes more interesting when MDAnalysis is used as glue
between different analysis packages and when one can arrange things so that
intermediate frames (typically in the PDB format) are not written to disk but
remain in memory via e.g. :mod:`cStringIO` buffers.


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

File parsing
------------

.. autoclass:: FORTRANReader
   :members:
.. autodata:: FORTRAN_format_regex


Data manipulation and handling
------------------------------

.. autofunction:: fixedwidth_bins


Strings
-------

.. autofunction:: convert_aa_code
.. autofunction:: parse_residue
.. autofunction:: conv_float


Class decorators
----------------

.. autofunction:: cached

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


.. versionchanged:: 0.11.0
   Moved mathematical functions into lib.mdamath
"""

__docformat__ = "restructuredtext en"

from six.moves import range, map
import six

import os
import os.path
import errno
from contextlib import contextmanager
import bz2
import gzip
import re
import io
import warnings
from functools import wraps
import numpy as np
import functools

from ..exceptions import StreamWarning


# Python 3.0, 3.1 do not have the builtin callable()
try:
    callable(list)
except NameError:
    # http://bugs.python.org/issue10518
    import collections

    def callable(obj):
        return isinstance(obj, collections.Callable)


def filename(name, ext=None, keep=False):
    """Return a new name that has suffix attached; replaces other extensions.

    :Arguments:
      *name*
           filename; extension is replaced unless keep=True;
           *name* can also be a :class:`NamedStream` (and its
           :attr:`NamedStream.name` will be changed accordingly)
      *ext*
           extension
      *keep*
           - ``False``: replace existing extension with *ext*;
           - ``True``: keep old extension if one existed

    .. versionchanged:: 0.9.0
       Also permits :class:`NamedStream` to pass through.
    """
    if ext is not None:
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
def openany(datasource, mode='r', reset=True):
    """Context manager for :func:`anyopen`.

    Open the *datasource* and close it when the context of the :keyword:`with`
    statement exits.

    *datasource* can be a filename or a stream (see :func:`isstream`). A stream
    is reset to its start if possible (via :meth:`~io.IOBase.seek` or
    :meth:`~cString.StringIO.reset`).

    The advantage of this function is that very different input sources
    ("streams") can be used for a "file", ranging from files on disk (including
    compressed files) to open file objects to sockets and strings---as long as
    they have a file-like interface.

    :Arguments:
     *datasource*
        a file or a stream
     *mode*
        'r' or 'w'
     *reset*
        try to read (*mode* 'r') the stream from the start [``True``]

    .. rubric:: Examples

    Open a gzipped file and process it line by line::

       with openany("input.pdb.gz") as pdb:
          for line in pdb:
              if line.startswith('ATOM'): print(line)

    Open a URL and read it::

       import urllib2
       with openany(urllib2.urlopen("http://www.MDAnalysis.org/")) as html:
          print(html.read())

    .. SeeAlso:: :func:`anyopen`
    """
    stream = anyopen(datasource, mode=mode, reset=reset)
    try:
        yield stream
    finally:
        stream.close()


# On python 3, we want to use bz2.open to open and uncompress bz2 files. That
# function allows to specify the type of the uncompressed file (bytes ot text).
# The function does not exist in python 2, thus we must use bz2.BZFile to
# which we cannot tell if the uncompressed file contains bytes or text.
# Therefore, on python 2 we use a proxy function that removes the type of the
# uncompressed file from the `mode` argument.
try:
    bz2.open
except AttributeError:
    # We are on python 2 and bz2.open is not available
    def bz2_open(filename, mode):
        """Open and uncompress a BZ2 file"""
        mode = mode.replace('t', '').replace('b', '')
        return bz2.BZ2File(filename, mode)
else:
    # We are on python 3 so we can use bz2.open
    bz2_open = bz2.open


def anyopen(datasource, mode='r', reset=True):
    """Open datasource (gzipped, bzipped, uncompressed) and return a stream.

    *datasource* can be a filename or a stream (see :func:`isstream`). By
    default, a stream is reset to its start if possible (via
    :meth:`~io.IOBase.seek` or :meth:`~cString.StringIO.reset`).

    If possible, the attribute ``stream.name`` is set to the filename or
    "<stream>" if no filename could be associated with the *datasource*.

    :Arguments:
     *datasource*
        a file (from :class:`file` or :func:`open`) or a stream (e.g. from
        :func:`urllib2.urlopen` or :class:`cStringIO.StringIO`)
     *mode*
        'r' or 'w' or 'a', more complicated modes ('r+', 'w+' are not supported because
        only the first letter is looked at) [``'r'``]
     *reset*
        try to read (*mode* 'r') the stream from the start [``True``]

    :Returns: tuple ``stream`` which is a file-like object

    .. SeeAlso:: :func:`openany` to be used with the :keyword:`with` statement.

    .. versionchanged:: 0.9.0
       Only returns the ``stream`` and tries to set ``stream.name = filename`` instead of the previous
       behavior to return a tuple ``(stream, filename)``.
    """
    handlers = {'bz2': bz2_open, 'gz': gzip.open, '': open}

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
                openfunc = handlers[ext]
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
            openfunc = handlers[ext]
            stream = openfunc(datasource, mode=mode)
            if stream is None:
                raise IOError(errno.EIO, "Cannot open file or stream in mode={mode!r}.".format(**vars()), repr(filename))
    else:
        raise NotImplementedError("Sorry, mode={mode!r} is not implemented for {datasource!r}".format(**vars()))
    try:
        stream.name = filename
    except (AttributeError, TypeError):
        pass  # can't set name (e.g. cStringIO.StringIO)
    return stream


def _get_stream(filename, openfunction=open, mode='r'):
    """Return open stream if *filename* can be opened with *openfunction* or else ``None``."""
    try:
        stream = openfunction(filename, mode=mode)
    except IOError:
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
    p : path, string

    Returns
    -------
    Tuple ``(root, extension)`` where ``root`` is the full path and
    filename with all extensions removed whereas ``extension`` is the
    string of all extensions.

    Example
    -------
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
    """Detect if *obj* is a stream.

    We consider anything a stream that has the methods

    - ``close()``

    and either set of the following

    - ``read()``, ``readline()``, ``readlines()``
    - ``write()``, ``writeline()``, ``writelines()``

    .. SeeAlso:: :mod:`io`

    :Arguments:
      *obj*
          stream or string

    :Returns: ``True`` is *obj* is a stream, ``False`` otherwise

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
    """Determine full path of executable *program* on :envvar:`PATH`.

    (Jay at http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python)
    """

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
class NamedStream(io.IOBase):
    """Stream that also provides a (fake) name.

    By wrapping a stream *stream* in this class, it can be passed to
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

    .. rubric:: Example

    Wrap a :func:`cStringIO.StringIO` instance to write to::

      import cStringIO
      import os.path
      stream = cStringIO.StringIO()
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

    .. Note::

       This class uses its own :meth:`__getitem__` method so if *stream*
       implements :meth:`stream.__getitem__` then that will be masked and this
       class should not be used.

    .. Warning::

       By default, :meth:`NamedStream.close` will **not close the
       stream** but instead :meth:`~NamedStream.reset` it to the
       beginning. [#NamedStreamClose]_ Provide the ``force=True`` keyword
       to :meth:`NamedStream.close` to always close the stream.

    """

    def __init__(self, stream, filename, reset=True, close=False):
        """Initialize the :class:`NamedStream` from a *stream* and give it a *name*.

        The constructor attempts to rewind the stream to the beginning unless
        the keyword *reset* is set to ``False``. If rewinding fails, a
        :class:`MDAnalysis.StreamWarning` is issued.

        .. Note::

           By default, this stream will *not* be closed by :keyword:`with` and
           :meth:`close` (see there) unless the *close* keyword is set to
           ``True``.

        Arguments
        ---------
        stream : stream
               an open stream (e.g. :class:`file` or :func:`cStringIO.StringIO`)
        filename : str
               the filename that should be associated with the stream

        Keywords
        --------
        reset : boolean, default ``True``
               start the stream from the beginning (either :meth:`reset` or :meth:`seek`)
               when the class instance is constructed
        close : booelan, default ``True``
               close the stream when a :keyword:`with` block exits or when
               :meth:`close` is called; note that the default is **not to close
               the stream**

        .. versionadded:: 0.9.0
        """
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
        """Change the stream position to the given byte *offset* .

        *offset* is interpreted relative to the position indicated by
        *whence*. Values for *whence* are:

        - :data:`io.SEEK_SET` or 0 – start of the stream (the default); *offset*
          should be zero or positive
        - :data:`io.SEEK_CUR` or 1 – current stream position; *offset* may be
          negative
        - :data:`io.SEEK_END` or 2 – end of the stream; *offset* is usually
          negative

        :Returns: the new absolute position.
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
        """Truncate the stream's size to *size*.

        The size defaults to the current position (if no *size* argument is
        supplied). The current file position is not changed.
        """
        try:
            return self.stream.truncate(*size)
        except AttributeError:
            return super(NamedStream, self).truncate(*size)

    def seekable(self):
        """Return ``True`` if the stream supports random access.

        If ``False``, :meth:`seek`, :meth:`tell` and :meth:`truncate` will raise :exc:`IOError`.
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
            raise IOError("This NamedStream does not use a file descriptor.")

    # fake the important parts of the string API
    # (other methods such as rfind() are automatically dealt with via __getattr__)
    def __getitem__(self, x):
        return self.name[x]

    def __eq__(self, x):
        return self.name == x

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
    """Join all args and return the real path, rooted at /.

    Expands '~', '~user', and environment variables such as :envvar`$HOME`.

    Returns ``None`` if any of the args is ``None``.
    """
    if None in args:
        return None
    return os.path.realpath(os.path.expanduser(os.path.expandvars(os.path.join(*args))))


def get_ext(filename):
    """Return the lower-cased extension of *filename* without a leading dot.

    :Returns: root, ext
    """
    root, ext = os.path.splitext(filename)
    if ext.startswith(os.extsep):
        ext = ext[1:]
    return root, ext.lower()


def check_compressed_format(root, ext):
    """Check if this is a supported gzipped/bzip2ed file format and return UPPERCASE format."""
    # XYZReader&others are setup to handle both plain and compressed (bzip2, gz) files
    # ..so if the first file extension is bzip2 or gz, look at the one to the left of it
    if ext.lower() in ("bz2", "gz"):
        try:
            root, ext = get_ext(root)
        except:
            raise TypeError("Cannot determine coordinate format for '{0}.{1}'"
                            "".format(root, ext))

    return ext.upper()


def format_from_filename_extension(filename):
    """Guess file format from the file extension"""
    try:
        root, ext = get_ext(filename)
    except:
        raise TypeError(
            "Cannot determine file format for file '{0}'.\n"
            "           You can set the format explicitly with "
            "'Universe(..., format=FORMAT)'.".format(filename))
    format = check_compressed_format(root, ext)
    return format


def guess_format(filename):
    """Return the format of *filename*

    The current heuristic simply looks at the filename extension
    and can work around compressed format extensions

    *filename* can also be a stream, in which case
    *filename.name* is looked at for a hint to the format

    :Raises:
       *ValueError*

    .. versionadded:: 0.11.0
       Moved into lib.util
    """
    if isstream(filename):
        # perhaps StringIO or open stream
        try:
            format = format_from_filename_extension(filename.name)
        except AttributeError:
            # format is None so we need to complain:
            raise ValueError("guess_format requires an explicit format specifier "
                             "for stream {0}".format(filename))
    else:
        # iterator, list, filename: simple extension checking... something more
        # complicated is left for the ambitious.
        # Note: at the moment the upper-case extension *is* the format specifier
        # and list of filenames is handled by ChainReader
        format = (format_from_filename_extension(filename)
                  if not iterable(filename) else 'CHAIN')

    return format.upper()


def iterable(obj):
    """Returns ``True`` if *obj* can be iterated over and is *not* a  string
    nor a :class:`NamedStream`"""
    if isinstance(obj, (six.string_types, NamedStream)):
        return False  # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True  # any iterator will do
    try:
        len(obj)  # anything else that might work
    except TypeError:
        return False
    return True


def asiterable(obj):
    """Returns obj so that it can be iterated over; a string is *not* treated as iterable"""
    if not iterable(obj):
        obj = [obj]
    return obj

#: Regular expresssion (see :mod:`re`) to parse a simple `FORTRAN edit descriptor`_.
#: ``(?P<repeat>\d?)(?P<format>[IFELAX])(?P<numfmt>(?P<length>\d+)(\.(?P<decimals>\d+))?)?``
#:
#: .. _FORTRAN edit descriptor: http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
FORTRAN_format_regex = "(?P<repeat>\d+?)(?P<format>[IFEAX])(?P<numfmt>(?P<length>\d+)(\.(?P<decimals>\d+))?)?"
_FORTRAN_format_pattern = re.compile(FORTRAN_format_regex)


def strip(s):
    """Convert *s* to a string and return it white-space stripped."""
    return str(s).strip()


class FixedcolumnEntry(object):
    """Represent an entry at specific fixed columns.

    Reads from line[start:stop] and converts according to
    typespecifier.
    """
    convertors = {'I': int, 'F': float, 'E': float, 'A': strip}

    def __init__(self, start, stop, typespecifier):
        """
        :Arguments:
         *start*
            first column
         *stop*
            last column + 1
         *typespecifier*
            'I': int, 'F': float, 'E': float, 'A': stripped string

        The start/stop arguments follow standard Python convention in that
        they are 0-based and that the *stop* argument is not included.
        """
        self.start = start
        self.stop = stop
        self.typespecifier = typespecifier
        self.convertor = self.convertors[typespecifier]

    def read(self, line):
        """Read the entry from *line* and convert to appropriate type."""
        try:
            return self.convertor(line[self.start:self.stop])
        except ValueError:
            raise ValueError("{0!r}: Failed to read&convert {1!r}".format(self, line[self.start:self.stop]))

    def __len__(self):
        """Length of the field in columns (stop - start)"""
        return self.stop - self.start

    def __repr__(self):
        return "FixedcolumnEntry({0:d},{1:d},{2!r})".format(self.start, self.stop, self.typespecifier)


class FORTRANReader(object):
    """FORTRANReader provides a method to parse FORTRAN formatted lines in a file.

    Usage::

       atomformat = FORTRANReader('2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10')
       for line in open('coordinates.crd'):
           serial,TotRes,resName,name,x,y,z,chainID,resSeq,tempFactor = atomformat.read(line)

    Fortran format edit descriptors; see `Fortran Formats`_ for the syntax.

    Only simple one-character specifiers supported here: *I F E A X* (see
    :data:`FORTRAN_format_regex`).

    Strings are stripped of leading and trailing white space.

    .. _`Fortran Formats`: http://www.webcitation.org/5xbaWMV2x
    .. _`Fortran Formats (URL)`:
       http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
    """

    def __init__(self, fmt):
        """Set up the reader with the FORTRAN format string.

        The string *fmt* should look like '2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10'.
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
        """Parse *line* according to the format string and return list of values.

        Values are converted to Python types according to the format specifier.

        :Returns: list of entries with appropriate types
        :Raises: :exc:`ValueError` if any of the conversions cannot be made
                 (e.g. space for an int)

        .. SeeAlso:: :meth:`FORTRANReader.number_of_matches`
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

          parse_FORTRAN_format(edit_descriptor) --> dict

        :Returns: dict with totallength (in chars), repeat, length,
                  format, decimals
        :Raises: :exc:`ValueError` if the *edit_descriptor* is not recognized
                 and cannot be parsed

        .. Note::

           Specifiers: *L ES EN T TL TR / r S SP SS BN BZ* are *not*
           supported, and neither are the scientific notation *Ew.dEe*
           forms.
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
    """Return bins of width delta that cover xmin,xmax (or a larger range).

    dict = fixedwidth_bins(delta,xmin,xmax)

    The dict contains 'Nbins', 'delta', 'min', and 'max'.
    """
    if not np.all(xmin < xmax):
        raise ValueError('Boundaries are not sane: should be xmin < xmax.')
    _delta = np.asarray(delta, dtype=np.float_)
    _xmin = np.asarray(xmin, dtype=np.float_)
    _xmax = np.asarray(xmax, dtype=np.float_)
    _length = _xmax - _xmin
    N = np.ceil(_length / _delta).astype(np.int_)  # number of bins
    dx = 0.5 * (N * _delta - _length)  # add half of the excess to each end
    return {'Nbins': N, 'delta': _delta, 'min': _xmin - dx, 'max': _xmax + dx}


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
amino_acid_codes = {one: three for three, one in canonical_inverse_aa_codes.items()}
#: non-default charge state amino acids or special charge state descriptions
#: (Not fully synchronized with :class:`MDAnalysis.core.Selection.ProteinSelection`.)
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

    .. SeeAlso:: Data are defined in :data:`amino_acid_codes` and :data:`inverse_aa_codes`.
    """
    if len(x) == 1:
        d = amino_acid_codes
    else:
        d = inverse_aa_codes

    try:
        return d[x.upper()]
    except KeyError:
        raise ValueError("No conversion for {0} found (1 letter -> 3 letter or 3/4 letter -> 1 letter)".format(x))


#: Regular expression to match and parse a residue-atom selection; will match
#: "LYS300:HZ1" or "K300:HZ1" or "K300" or "4GB300:H6O" or "4GB300" or "YaA300".
RESIDUE = re.compile("""
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

    Examples:
     - "LYS300:HZ1" --> ("LYS", 300, "HZ1")
     - "K300:HZ1" --> ("LYS", 300, "HZ1")
     - "K300" --> ("LYS", 300, None)
     - "4GB300:H6O" --> ("4GB", 300, "H6O")
     - "4GB300" --> ("4GB", 300, None)

    :Argument: The *residue* must contain a 1-letter or 3-letter or
               4-letter residue string, a number (the resid) and
               optionally an atom identifier, which must be separate
               from the residue with a colon (":"). White space is
               allowed in between.

    :Returns: `(3-letter aa string, resid, atomname)`; known 1-letter
              aa codes are converted to 3-letter codes
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
    """Convert an object *s* to float if possible.

    Function to be passed into :func:`map` or a list comprehension. If
    the argument can be interpreted as a float it is converted,
    otherwise the original object is passed back.
    """
    try:
        return float(s)
    except ValueError:
        return s


def cached(key):
    """Cache a property within a class

    Requires the Class to have a cache dict called ``_cache``.

    Usage::

       class A(object):
           def__init__(self):
               self._cache = dict()

           @property
           @cached('keyname')
           def size(self):
               # This code gets ran only if the lookup of keyname fails
               # After this code has been ran once, the result is stored in
               # _cache with the key: 'keyname'
               size = 10.0

    .. versionadded:: 0.9.0
    """

    def cached_lookup(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            try:
                return self._cache[key]
            except KeyError:
                self._cache[key] = ret = func(self, *args, **kwargs)
                return ret

        return wrapper

    return cached_lookup


def blocks_of(a, n, m):
    """Extract a view of (n, m) blocks along the diagonal of the array `a`

    Parameters
    ----------
    a : array_like
        starting array
    n : int
        size of block in first dimension
    m : int
        size of block in second dimension


    Returns
    -------
      (nblocks, n, m) view of the original array.
      Where nblocks is the number of times the miniblock fits in the original.

    Raises
    ------
      ValueError
        If the supplied `n` and `m` don't divide `a` into an integer number
        of blocks.

    Examples
    --------
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
      n, m must divide a into an identical integer number of blocks.

      Uses strides so probably requires that the array is C contiguous

      Returns a view, so editing this modifies the original array

    .. versionadded:: 0.12.0
    """
    # based on:
    # http://stackoverflow.com/a/10862636
    # but generalised to handle non square blocks.

    nblocks = a.shape[0] / n
    nblocks2 = a.shape[1] / m

    if not nblocks == nblocks2:
        raise ValueError("Must divide into same number of blocks in both"
                         " directions.  Got {} by {}"
                         "".format(nblocks, nblocks2))

    new_shape = (nblocks, n, m)
    new_strides = (n * a.strides[0] + m * a.strides[1],
                   a.strides[0], a.strides[1])

    return np.lib.stride_tricks.as_strided(a, new_shape, new_strides)
