# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

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

.. autofunction:: triclinic_box
.. autofunction:: triclinic_vectors
.. autofunction:: box_volume

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


def get_reader_for(filename, permissive=False, format=None):
    """Return the appropriate trajectory reader class for *filename*.

    Automatic detection is disabled when an explicit *format* is
    provided.
    """
    format = guess_format(filename, format=format)
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
    if isinstance(filename,tuple):
        Reader = get_reader_for(filename[0],permissive=kwargs.pop('permissive', False),format=filename[1])
        return Reader(filename[0], **kwargs)
    else:
        Reader = get_reader_for(filename, permissive=kwargs.pop('permissive', False))
        return Reader(filename, **kwargs)

def get_writer_for(filename=None, format='DCD', multiframe=None):
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
      *multiframe*
         ``True``: write multiple frames to the trajectory; ``False``: only
         write a single coordinate frame; ``None``: first try trajectory (multi
         frame writers), then the single frame ones. Default is ``None``.

    .. versionchanged:: 0.7.6
       Added *multiframe* keyword; the default ``None`` reflects the previous
       behaviour.
    """
    if filename:
        root, ext = get_ext(filename)
        format = check_compressed_format(root, ext)
    if multiframe is None:
        try:
            return MDAnalysis.coordinates._trajectory_writers[format]
        except KeyError:
            try:
                return MDAnalysis.coordinates._frame_writers[format]
            except KeyError:
                raise TypeError("No trajectory or frame writer for format %r" % format)
    elif multiframe is True:
        try:
            return MDAnalysis.coordinates._trajectory_writers[format]
        except KeyError:
            raise TypeError("No trajectory  writer for format %r" % format)
    elif multiframe is False:
        try:
            return MDAnalysis.coordinates._frame_writers[format]
        except KeyError:
            raise TypeError("No single frame writer for format %r" % format)
    else:
        raise ValueError("Unknown value %r for multiframe, only True, False, None allowed" % multiframe)

def writer(filename, numatoms=None, **kwargs):
    """Initialize a trajectory writer instance for *filename*.

    :Arguments:
       *filename*
            Output filename of the trajectory; the extension determines the
            format.
       *numatoms*
            The number of atoms in the output trajectory; can be ommitted
            for single-frame writers.
       *multiframe*
            ``True``: write a trajectory with multiple frames; ``False``
            only write a single frame snapshot; ``None`` first try to get
            a multiframe writer and then fall back to single frame [``None``]
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

    .. versionchanged:: 0.7.6
       Added *multiframe* keyword. See also :func:`get_writer_for`.
    """
    Writer = get_writer_for(filename, format=kwargs.pop('format',None),
                            multiframe=kwargs.pop('multiframe', None))
    return Writer(filename, numatoms=numatoms, **kwargs)

def get_ext(filename):
    """Return the lower-cased extension of *filename* without a leading dot.

    :Returns: root, ext
    """
    root, ext = os.path.splitext(filename)
    if ext.startswith(os.extsep):
        ext = ext[1:]
    return root, ext.lower()

def guess_format(filename, format=None):
    """Returns the type of coordinate file *filename*.

    The current heuristic simply looks at the filename extension but
    more complicated probes could be implemented here or in the
    individual packages (e.g. as static methods).

    If *format* is supplied then it overrides the auto detection.
    """
    if format is None:
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
    else:
        # internally, formats are all uppercase
        # enable chain reader with explicit format
        if MDAnalysis.core.util.iterable(filename):
            # list of filenames, handled by ChainReader
            format = 'CHAIN'
        else:
            format = str(format).upper()

    # sanity check
    if format != 'CHAIN' and not format in MDAnalysis.coordinates._trajectory_readers:
        raise TypeError("Unknown coordinate trajectory format %r for %r; only %r are implemented in MDAnalysis." %
                        (format, filename, MDAnalysis.coordinates._trajectory_readers.keys()))
    return format

def check_compressed_format(root, ext):
    """Check if this is a supported gzipped/bzip2ed file format and return UPPERCASE format."""
    filename = root + '.' + ext  # only needed for diagnostics
    # XYZReader&others are setup to handle both plain and compressed (bzip2, gz) files
    # ..so if the first file extension is bzip2 or gz, look at the one to the left of it
    if ext.lower() in ("bz2","gz"):
        try:
            root, ext = get_ext(root)
        except:
            raise TypeError("Cannot determine coordinate format for %r" % filename)
        if not ext.upper() in MDAnalysis.coordinates._compressed_formats:
            raise TypeError("Cannot handle %r in compressed format" % filename)
    return ext.upper()

def _veclength(v):
    """Length of vector *v*."""
    # note: this is 3 times faster than numpy.linalg.norm
    return numpy.sqrt(numpy.dot(v,v))

def _angle(a,b):
    """Angle between two vectors *a* and *b* in degrees.

    If one of the lengths is 0 then the angle is returned as 0
    (instead of `nan`).
    """
    angle = numpy.arccos(numpy.dot(a,b) / (_veclength(a)*_veclength(b)))
    if numpy.isnan(angle):
        return 0.0
    return rad2deg(angle)

def triclinic_box(x,y,z):
    """Convert the three triclinic box vectors to [A,B,C,alpha,beta,gamma].

    Angles are in degrees.

    * alpha  = angle(y,z)
    * beta   = angle(x,z)
    * gamma  = angle(x,y)

    .. SeeAlso:: Definition of angles: http://en.wikipedia.org/wiki/Lattice_constant
    """
    A, B, C = [_veclength(v) for v in x,y,z]
    alpha =  _angle(y,z)
    beta  =  _angle(x,z)
    gamma =  _angle(x,y)
    return numpy.array([A,B,C,alpha,beta,gamma], dtype=numpy.float32)

def triclinic_vectors(dimensions):
    """Convert `[A,B,C,alpha,beta,gamma]` to a triclinic box representation.

    Original `code by Tsjerk Wassenaar`_ posted on the Gromacs mailinglist.

    If *dimensions* indicates a non-periodic system (i.e. all lengths
    0) then null-vectors are returned.

    .. _code by Tsjerk Wassenaar:
       http://www.mail-archive.com/gmx-users@gromacs.org/msg28032.html

    :Arguments:
      *dimensions*
        list of box lengths and angles (in degrees) such as
        [A,B,C,alpha,beta,gamma]

    :Returns: numpy 3x3 array B, with B[0] = first box vector,
              B[1] = second vector, B[2] third box vector.

    .. note::

       The first vector is always pointing along the X-axis
       i.e. parallel to (1,0,0).

    .. versionchanged:: 0.7.6
       Null-vectors are returned for non-periodic (or missing) unit cell.

    """
    B = numpy.zeros((3,3), dtype=numpy.float32)
    x, y, z, a, b, c = dimensions[:6]

    if numpy.all(dimensions[:3] == 0):
        return B

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

def box_volume(dimensions):
    """Return the volume of the unitcell described by *dimensions*.

    The volume is computed as `det(x1,x2,x2)` where the xi are the
    triclinic box vectors from :func:`triclinic_vectors`.

    :Arguments:
       *dimensions*
          list of box lengths and angles (in degrees) such as
          [A,B,C,alpha,beta,gamma]

    :Returns: numpy 3x3 array B, with B[0] = first box vector,
              B[1] = second vector, B[2] third box vector.
    """
    return numpy.linalg.det(triclinic_vectors(dimensions))
