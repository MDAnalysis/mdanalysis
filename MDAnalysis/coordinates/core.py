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
    
    # XYZReader is setup to handle both plain and compressed (bzip2, gz) files
    # ..so if the first file extension is bzip2 or gz, look at the one to the left of it 
    if ext in ("bz2","gz"):
	    root, ext = os.path.splitext(root)
	    try:
		  if ext.startswith('.'):
			ext = ext[1:]
		  ext = ext.lower()
	    except:
		  raise TypeError("Cannot determine coordinate format for %r" % filename)	
	    if ext != "xyz": # only bzipped xyz files can be parsed right now (might be useful to parse foo.pdb.bz2 ?) 
		  raise TypeError("Cannot handle %r in bzipped format" % filename)    		      		    	
   
    if not ext in MDAnalysis.coordinates._trajectory_readers:
        raise TypeError("Unknown coordinate trajectory extension %r from %r; only %r are implemented in MDAnalysis." % 
                        (ext, filename, MDAnalysis.coordinates._trajectory_readers.keys()))
    return ext
    
def _veclength(v):
    """Length of vector *v*."""
    return numpy.sqrt(numpy.dot(v,v))

def _angle(a,b):
    """Angle between two vectors *a* and *b* in degrees."""
    angle = numpy.arccos(numpy.dot(a,b) / (_veclength(a)*_veclength(b)))
    return rad2deg(angle)



