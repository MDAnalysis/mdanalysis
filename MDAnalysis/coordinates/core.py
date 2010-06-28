"""
:mod:`MDAnalysis.coordinates.core` --- Common functions for coordinate reading
==============================================================================

Important base classes are collected in :mod:`MDAnalysis.coordinates.base`.

.. function:: get_reader_for
.. function:: guess_format

"""
import os.path
import MDAnalysis.coordinates

def get_reader_for(filename):
    """Return the appropriate topology parser for *filename*."""
    return MDAnalysis.coordinates._trajectory_readers[guess_format(filename)]

def guess_format(filename):
    """Returns the type of coordinate file *filename*.

    The current heuristic simply looks at the filename extension but
    more complicated probes could be implemented here or in the
    individual packages (e.g. as static methods).
    """
    # simple extension checking... something more complicated is left
    # for the ambitious
    root, ext = os.path.splitext(filename)
    try:
        if ext.startswith('.'):
            ext = ext[1:]
        ext = ext.lower()
    except:
        raise TypeError("Cannot determine coordinate format for %r" % filename)
    
    if not ext in MDAnalysis.coordinates._trajectory_readers:
        raise TypeError("Unknown coordinate trajectory extension %r from %r; only %r are implemented in MDAnalysis." % 
                        (ext, filename, MDAnalysis.coordinates._trajectory_readers.keys()))
    return ext



