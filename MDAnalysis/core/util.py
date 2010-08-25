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


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
from contextlib import contextmanager
import bz2, gzip


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
