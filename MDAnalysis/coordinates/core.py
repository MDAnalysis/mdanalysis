"""
Common functions for coordinate reading --- :mod:`MDAnalysis.coordinates.core`  
==============================================================================

Important base classes are collected in :mod:`MDAnalysis.coordinates.base`.

.. autofunction:: reader
.. autofunction:: writer

Helper functions:

.. autofunction:: get_reader_for
.. autofunction:: get_writer_for
.. autofunction:: guess_format

"""
import os.path
import MDAnalysis.coordinates
import MDAnalysis.core.util

import numpy

from numpy import sin, cos, sqrt
try:
    from numpy import rad2deg, deg2rad   # numpy 1.3+
except ImportError:
    def rad2deg(x):             # no need for the numpy out=[] argument 
        return 180.0*x/numpy.pi
    def deg2rad(x):             # no need for the numpy out=[] argument 
        return x*numpy.pi/180.0


def get_reader_for(filename, permissive=False):
    """Return the appropriate trajectory reader class for *filename*."""
    format = guess_format(filename)
    if permissive:
        return MDAnalysis.coordinates._trajectory_readers_permissive[format]
    return MDAnalysis.coordinates._trajectory_readers[format]

def reader(filename, **kwargs):
    """Provide a trajectory reader instance for *filename*.
    
    This function guesses the file format from the extension of *filename* and
    it will throw a :exc:`TypeError` if the extension is not recognized.

    In most cases, no special keyword arguments are necessary. For some readers
    (such as PDB) it might be useful to set the *permissive* = ``True`` flag to
    select a simpler but faster reader.

    All other keywords are passed on to the underlying Reader classes; see
    their documentation for details.

    .. SeeAlso:: For trajectory formats: :class:`~DCD.DCDReader`,
       :class:`~XTC.XTCReader`, :class:`~TRR.TRRReader`,
       :class:`~XYZ.XYZReader`.  For single frame formats:
       :class:`~CRD.CRDReader`, :class:`~PDB.PDBReader` and
       :class:`~PDB.PrimitivePDBReader`, :class:`~GRO.GROReader`,

    :Arguments:
       *filename*
          filename of the input trajectory or coordinate file
       *permissive*
          If set to ``True``, a reader is selected that is more tolerant of the
          input (currently only implemented for PDB). [``False``]
       *kwargs*
           Keyword arguments for the selected Reader class.
    """
    Reader = get_reader_for(filename, permissive=kwargs.pop('permissive', False))
    return Reader(filename, **kwargs)

def get_writer_for(filename=None, format='DCD'):
    """Return an appropriate trajectory or frame writer class for *filename*.
    
    The format is determined by the *format* argument or the extension
    of *filename*. The default is to return a dcd writer (*format* = 'dcd').

    :Arguments:
      *filename*
         The filename for the trajectory is examined for its extension and 
         the Writer is chosen accordingly.
      *format*
         If no *filename* is supplied then the format can be explicitly set;
         possible values are "DCD", "XTC", "TRR"; "PDB", "CRD", "GRO".
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

def writer(filename, numatoms=None, **kwargs):
    """Initialize a trajectory writer instance for *filename*.

    :Arguments:
       *filename*
           Output filename of the trajectory; the extension determines the 
           format.
       *numatoms*
            The number of atoms in the output trajectory; can be ommitted 
            for single-frame writers.
       *kwargs*
            Keyword arguments for the writer; all trajectory Writers accept
            at least

               *start*
                   starting time [0]
               *step*
                   step size in frames [1]
               *delta*
                   length of time between two frames, in ps [1.0]

            Some readers accept additional arguments, which need to be looked
            up in the documentation of the reader.

            .. SeeAlso:: :class:`~MDAnalysis.coordinates.DCD.DCDWriter` for DCD
               trajectories or :class:`~MDAnalysis.coordinates.XTC.XTCWriter`
               and :class:`~MDAnalysis.coordinates.TRR.TRRWriter` for Gromacs.
    """
    Writer = get_writer_for(filename, format=kwargs.pop('format',None))
    return Writer(filename, numatoms=numatoms, **kwargs)

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

    if MDAnalysis.core.util.iterable(filename):
        # list of filenames, handled by ChainReader
        format = 'CHAIN'
    else:
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
        if not ext.upper() in MDAnalysis.coordinates._compressed_formats: 
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

def triclinic_box(x,y,z):
    """Convert the three triclinic box vectors to [A,B,C,alpha,beta,gamma].

    Angles are in degrees.
    """
    A, B, C = [_veclength(v) for v in x,y,z]
    alpha =  _angle(x,y)
    beta  =  _angle(x,z)
    gamma =  _angle(y,z)
    return numpy.array([A,B,C,alpha,beta,gamma], dtype=numpy.float32)

def triclinic_vectors(dimensions):
	"""Convert [A,B,C,alpha,beta,gamma] to a triclinic box representation.

	Original code by Tsjerk Wassenaar; see http://www.mail-archive.com/gmx-users@gromacs.org/msg28032.html

        :Arguments:
          *dimensions*
             list of box lengths and angles (in degrees) such as
             [A,B,C,alpha,beta,gamma]

        :Returns: numpy 3x3 array B, with B[0] = first box vector,
                  B[1] = second vector, B[2] third box vector.

        .. note:: The first vector is always pointing along the
                  X-axis i.e. parallel to (1,0,0).
	"""
	B = numpy.zeros((3,3), dtype=numpy.float32)
	x, y, z, a, b, c = dimensions[:6]

	B[0][0] = x
	if a == 90. and b == 90. and c == 90.:
		B[1][1] = y
		B[2][2] = z
	else:
		a = deg2rad(a)
		b = deg2rad(b)
		c = deg2rad(c)
		B[1][0] = y*cos(c)
		B[1][1] = y*sin(c)
		B[2][0] = z*cos(b)
		B[2][1] = z*(cos(a)-cos(b)*cos(c))/sin(c)
		B[2][2] = sqrt(z*z-B[2][0]**2-B[2][1]**2)
        return B



