# $Id$
"""
Helper functions
================

Small helper functions that don't fit anywhere else.

Files and directories
---------------------

.. autofunction:: filename
.. function:: openany(directory[,mode='r'])

   Context manager to open a compressed (bzip2, gzip) or plain file
   (uses :func:`anyopen`).

.. autofunction:: anyopen


Containers and lists
--------------------

.. autofunction:: iterable
.. autofunction:: asiterable


File parsing
------------

.. autoclass:: FORTRANReader
.. autodata:: FORTRAN_format_regex

"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
from contextlib import contextmanager
import bz2, gzip
import re

def filename(name,ext=None,keep=False):
    """Return a new name that has suffix attached; replaces other extensions.

    name        filename; extension is replaced unless keep=True
    ext         extension
    keep        False: replace existing extension; True: keep if exists
    """ 
    name = str(name)
    if ext is None:
        return name
    if not ext.startswith(os.path.extsep):
        ext = os.path.extsep + ext
    #if name.find(ext) > 0:    # normally >= 0 but if we start with '.' then we keep it
    #    return name
    root, origext = os.path.splitext(name)
    if len(origext) == 0 or not keep:
        return root + ext
    return name

@contextmanager
def openany(datasource, mode='r'):
    """Open the datasource and close it when the context exits."""
    stream, filename = anyopen(datasource, mode=mode)
    try:
        yield stream
    finally:
        stream.close()

def anyopen(datasource, mode='r'):
    """Open datasource (gzipped, bzipped, uncompressed) and return a stream.

    :Arguments:
    - *datasource*: a file or a stream
    - *mode*: 'r' or 'w'
    """
    # TODO: - make this act as ContextManager (and not return filename)
    #       - need to add binary 'b' to mode for compressed files?

    handlers = {'bz2': bz2.BZ2File, 'gz': gzip.open, '': file}
    
    if mode.startswith('r'):
        if hasattr(datasource,'next') or hasattr(datasource,'readline'):
            stream = datasource
            filename = '(%s)' % stream.name  # maybe that does not always work?
        else:
            stream = None
            filename = datasource
            for ext in ('bz2', 'gz', ''):   # file == '' should be last
                openfunc = handlers[ext]
                stream = _get_stream(datasource, openfunc, mode=mode)
                if not stream is None:
                    break
            if stream is None:
                raise IOError("Cannot open %(filename)r in mode=%(mode)r." % vars())
    elif mode.startswith('w'):
        if hasattr(datasource, 'write'):
            stream = datasource
            filename = '(%s)' % stream.name  # maybe that does not always work?
        else:
            stream = None
            filename = datasource
            name, ext = os.path.splitext(filename)
            if ext.startswith('.'):
                ext = ext[1:]
            if not ext in ('bz2', 'gz'):
                ext = ''   # anything else but bz2 or gz is just a normal file
            openfunc = handlers[ext]
            stream = openfunc(datasource, mode=mode)
            if stream is None:
                raise IOError("Cannot open %(filename)r in mode=%(mode)r with type %(ext)r." % vars())
    else:
        raise NotImplementedError("Sorry, mode=%(mode)r is not implemented for %(datasource)r" % vars())

    return stream, filename

def _get_stream(filename, openfunction=file, mode='r'):
    try:
        stream = openfunction(filename, mode=mode)
    except IOError:
        return None

    try:
        stream.readline()
        stream.close()
        stream = openfunction(filename,'r')
    except IOError:
        stream.close()
        stream = None
    return stream


def iterable(obj):
    """Returns ``True`` if *obj* can be iterated over and is *not* a  string."""
    if type(obj) is str:
        return False    # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True    # any iterator will do 
    try: 
        len(obj)       # anything else that might work
    except TypeError: 
        return False
    return True

def asiterable(obj):
    """Returns obj so that it can be iterated over; a string is *not* treated as iterable"""
    if not iterable(obj):
        obj = [obj]
    return obj

#: Regular expresssion (see :mod:`re`) to parse a simple `FORTRAN edit descriptor`_.
#:   ``(?P<repeat>\d?)(?P<format>[IFELAX])(?P<numfmt>(?P<length>\d+)(\.(?P<decimals>\d+))?)?``
#: .. _FORTRAN edit descriptor: http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
FORTRAN_format_regex = "(?P<repeat>\d?)(?P<format>[IFEAX])(?P<numfmt>(?P<length>\d+)(\.(?P<decimals>\d+))?)?"
_FORTRAN_format_pattern = re.compile(FORTRAN_format_regex)

def strip(s):
    """Convert *s* to a string and return it whit-space stripped."""
    return str(s).strip()

class FixedcolumnEntry(object):
    convertors = {'I': int, 'F': float, 'E': float, 'A': strip}
    def __init__(self, start, stop, typespecifier):
        """Read string  line[start:stop] and convert according to typespecifier"""
        self.start = start
        self.stop = stop
        self.typespecifier = typespecifier
        self.convertor = self.convertors[typespecifier]
    def read(self, line):
        return self.convertor(line[self.start:self.stop])
    def __repr__(self):
        return "FixedcolumnEntry(%d,%d,%r)" % (self.start, self.stop, self.typespecifier)

class FORTRANReader(object):
    """FORTRANReader provides a method to parse FORTRAN formatted lines in a file.

    Usage::
       atomformat = FORTRANReader('2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10')
       for line in open('coordinates.crd'):
           serial,TotRes,resName,name,x,y,z,chainID,resSeq,tempFactor = atomformat.read(line)

    Fortran format edit descriptors: See
    http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
    for the syntax.

    Only simple one-character specifiers supported here: *I F E A X* (see
    :data:`FORTRAN_format_regex`).

    Strings are stripped of leading and trailing white space.
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
                    self.entries.append(FixedcolumnEntry(start,stop,d['format']))
                    start = stop
            else:
                start += d['totallength']

    def read(self, line):
        """Parse *line* according to the format string and return list of values.

        Values are converted to Python types according to the format specifier."""
        return [e.read(line) for e in self.entries]

    def parse_FORTRAN_format(self, edit_descriptor):
        """Parse the descriptor.

          parse_FORTRAN_format(edit_descriptor) --> dict

        :Returns: dict with totallength (in chars), repeat, length,
                  format, decimals; if the specifier could not be parsed,
                  an empty dict is returned.
        :Raises: :exc:`ValueError` if the *edit_descriptor* is not recognized.

        .. Note:: Specifiers: *L ES EN T TL TR / r S SP SS BN BZ* are *not*
           supported, and neither are the scientific notation *Ew.dEe* forms.
        """

        m = _FORTRAN_format_pattern.match(edit_descriptor)
        if m is None:
            raise ValueError("unrecognized FORTRAN format %r" % edit_descriptor)
        d = m.groupdict()
        if d['repeat'] == '':
            d['repeat'] = 1
        if d['format'] == 'X':
            d['length'] = 1
        for k in ('repeat','length','decimals'):
            try:
                d[k] = int(d[k])
            except ValueError:   # catches ''
                d[k] = 0
            except TypeError:    # keep None
                pass
        d['totallength'] = d['repeat'] * d['length']
        return d

