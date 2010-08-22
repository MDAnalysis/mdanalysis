"""
:mod:`MDAnalysis.coordinates.core` --- Common functions for coordinate reading
==============================================================================

Important base classes are collected in :mod:`MDAnalysis.coordinates.base`.

.. function:: get_reader_for
.. function:: guess_format

"""
import os.path
import MDAnalysis.coordinates
import numpy

try:
    from numpy import rad2deg   # numpy 1.3+
except ImportError:
    def rad2deg(x):             # no need for the numpy out=[] argument 
        return 180.0*x/numpy.pi

def get_reader_for(filename):
    """Return the appropriate trajectory reader for *filename*."""
    return MDAnalysis.coordinates._trajectory_readers[guess_format(filename)]

def get_writer_for(filename=None, format='DCD'):
    """Return an appropriate trajectory or frame writer for *filename*.
    
    The format is determined by the *format* argument or the extension
    of *filename*. The default is to return a dcd writer (*format* = 'dcd').
    """
    if filename:
        root, ext = get_ext(filename)
        format = check_compressed_format(root, ext)
    try:
        return MDAnalysis.coordinates._trajectory_writers[format]
    except KeyError:
        try:
            return MDAnalysis.coordinates._frame_writers[format]
        except KeyError:
            raise TypeError("No trajectory or frame writer for format %r" % format)

def get_ext(filename):
    """Return the lower-cased extension of *filename* without a leading dot.

    :Returns: root, ext
    """
    root, ext = os.path.splitext(filename)
    if ext.startswith('.'):
        ext = ext[1:]
    return root, ext.lower()

def guess_format(filename):
    """Returns the type of coordinate file *filename*.

    The current heuristic simply looks at the filename extension but
    more complicated probes could be implemented here or in the
    individual packages (e.g. as static methods).
    """
    # simple extension checking... something more complicated is left
    # for the ambitious

    # Note: at the moment the upper-case extension *is* the format specifier

    try:
        root, ext = get_ext(filename)
    except:
        raise TypeError("Cannot determine coordinate format for %r" % filename)
    format = ext.upper()
    format = check_compressed_format(root, ext)    
    if not format in MDAnalysis.coordinates._trajectory_readers:
        raise TypeError("Unknown coordinate trajectory extension %r from %r; only %r are implemented in MDAnalysis." % 
                        (ext, filename, MDAnalysis.coordinates._trajectory_readers.keys()))
    return format
    
def check_compressed_format(root, ext):
    """Check if this is a supported gzipped/bzip2ed file format and return format."""
    filename = root + '.' + ext  # only needed for diagnostics
    # XYZReader is setup to handle both plain and compressed (bzip2, gz) files
    # ..so if the first file extension is bzip2 or gz, look at the one to the left of it 
    if ext.lower() in ("bz2","gz"):
        try:
            root, ext = get_ext(root) 
        except:
            raise TypeError("Cannot determine coordinate format for %r" % filename)
        if ext.upper() != "XYZ": 
            # only bzipped xyz files can be parsed right now (might be useful to parse foo.pdb.bz2 ?) 
            raise TypeError("Cannot handle %r in compressed format" % filename)
    return ext.upper()

def _veclength(v):
    """Length of vector *v*."""
    return numpy.sqrt(numpy.dot(v,v))

def _angle(a,b):
    """Angle between two vectors *a* and *b* in degrees."""
    angle = numpy.arccos(numpy.dot(a,b) / (_veclength(a)*_veclength(b)))
    return rad2deg(angle)



