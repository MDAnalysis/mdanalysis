# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
import numpy
from numpy import sin, cos, sqrt
from numpy import rad2deg, deg2rad

import MDAnalysis.coordinates
import MDAnalysis.lib.util
from ..lib.mdamath import triclinic_box, triclinic_vectors, box_volume


def get_reader_for(filename, permissive=False, format=None):
    """Return the appropriate trajectory reader class for *filename*.

    Automatic detection is disabled when an explicit *format* is
    provided.
    """
    format = guess_format(filename, format=format)
    if permissive and format == 'PDB':
        return MDAnalysis.coordinates._trajectory_readers['Permissive_PDB']
    return MDAnalysis.coordinates._trajectory_readers[format]


def reader(filename, **kwargs):
    """Provide a trajectory reader instance for *filename*.

    This function guesses the file format from the extension of *filename* and
    it will throw a :exc:`TypeError` if the extension is not recognized.

    In most cases, no special keyword arguments are necessary. For PDB readers
    it might be useful to set the *permissive* = ``True`` flag to
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
    if isinstance(filename, tuple):
        Reader = get_reader_for(filename[0], permissive=kwargs.pop('permissive', False), format=filename[1])
        return Reader(filename[0], **kwargs)
    else:
        Reader = get_reader_for(filename, permissive=kwargs.pop('permissive', False))
        return Reader(filename, **kwargs)


def get_writer_for(filename=None, format='DCD', multiframe=None):
    """Return an appropriate trajectory or frame writer class for *filename*.

    The format is determined by the *format* argument or the extension of
    *filename*. The default is to return a dcd writer (*format* = 'dcd'). If
    the *filename* is not provided or if it is something like a
    :class:`cStringIO.StringIO` instance then the *format* argument must be
    used.

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
    if isinstance(filename, basestring) and filename:
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
    Writer = get_writer_for(filename, format=kwargs.pop('format', None),
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


def format_from_filename_extension(filename):
    """Guess file format from the file extension"""
    try:
        root, ext = get_ext(filename)
    except:
        raise TypeError(
            "Cannot determine coordinate file format for file '{0}'.\n"
            "           You can set the format explicitly with "
            "'Universe(..., format=FORMAT)'.".format(filename))
            #TypeError: ...."
    format = ext.upper()
    format = check_compressed_format(root, ext)
    return format


def guess_format(filename, format=None):
    """Returns the type of file *filename*.

    The current heuristic simply looks at the filename extension but more
    complicated probes could be implemented here or in the individual packages
    (e.g. as static methods). *filename* can also be a stream, in which case
    *filename.name* is looked at for a hint to the format if *format* is not
    provided.

    If *format* is supplied then it overrides the auto detection.
    """
    if format is None:
        if MDAnalysis.lib.util.isstream(filename):
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
            format = format_from_filename_extension(filename) if not MDAnalysis.lib.util.iterable(
                filename) else 'CHAIN'
    else:
        # format was set; but a list of filenames is always handled by ChainReader
        format = format if not MDAnalysis.lib.util.iterable(filename) else 'CHAIN'

    format = str(format).upper()

    # sanity check
    if format != 'CHAIN' and not format in MDAnalysis.coordinates._trajectory_readers:
        raise TypeError(
            "Unknown coordinate trajectory format '{0}' for '{1}'. The FORMATs \n"
            "           {2}\n"
            "           are implemented in MDAnalysis.\n"
            "           See http://docs.mdanalysis.org/documentation_pages/coordinates/init.html#id1\n"
            "           Use the format keyword to explicitly set the format: 'Universe(...,format=FORMAT)'\n"
            "           For missing formats, raise an issue at "
            "http://issues.mdanalysis.org".format(
                format, filename, MDAnalysis.coordinates._trajectory_readers.keys()))
            #TypeError: ...."
    return format


def check_compressed_format(root, ext):
    """Check if this is a supported gzipped/bzip2ed file format and return UPPERCASE format."""
    filename = root + '.' + ext  # only needed for diagnostics
    # XYZReader&others are setup to handle both plain and compressed (bzip2, gz) files
    # ..so if the first file extension is bzip2 or gz, look at the one to the left of it
    if ext.lower() in ("bz2", "gz"):
        try:
            root, ext = get_ext(root)
        except:
            raise TypeError("Cannot determine coordinate format for '{0}'".format(filename))
        if not ext.upper() in MDAnalysis.coordinates._compressed_formats:
            raise TypeError("Cannot handle coordinates '{0}' in compressed format".format(filename))
    return ext.upper()


