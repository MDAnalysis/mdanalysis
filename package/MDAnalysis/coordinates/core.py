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

"""

from ..lib import util
from ..lib.mdamath import triclinic_box, triclinic_vectors, box_volume


def get_reader_for(filename, permissive=False, format=None):
    """Return the appropriate trajectory reader class for *filename*.

    Automatic detection is disabled when an explicit *format* is
    provided.
    """
    from . import _trajectory_readers

    if format is None:
        format = util.guess_format(filename)
    format = format.upper()
    if permissive and format == 'PDB':
        return _trajectory_readers['Permissive_PDB']
    try:
        return _trajectory_readers[format]
    except KeyError:
        raise ValueError(
            "Unknown coordinate trajectory format '{0}' for '{1}'. The FORMATs \n"
            "           {2}\n"
            "           are implemented in MDAnalysis.\n"
            "           See http://docs.mdanalysis.org/documentation_pages/coordinates/init.html#id1\n"
            "           Use the format keyword to explicitly set the format: 'Universe(...,format=FORMAT)'\n"
            "           For missing formats, raise an issue at "
            "http://issues.mdanalysis.org".format(
                format, filename, _trajectory_readers.keys()))


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
    from . import _trajectory_writers, _frame_writers

    if isinstance(filename, basestring) and filename:
        root, ext = util.get_ext(filename)
        format = util.check_compressed_format(root, ext)
    if multiframe is None:
        try:
            return _trajectory_writers[format]
        except KeyError:
            try:
                return _frame_writers[format]
            except KeyError:
                raise TypeError("No trajectory or frame writer for format %r" % format)
    elif multiframe is True:
        try:
            return _trajectory_writers[format]
        except KeyError:
            raise TypeError("No trajectory  writer for format %r" % format)
    elif multiframe is False:
        try:
            return _frame_writers[format]
        except KeyError:
            raise TypeError("No single frame writer for format %r" % format)
    else:
        raise ValueError("Unknown value %r for multiframe, only True, False, None allowed" % multiframe)


def writer(filename, n_atoms=None, **kwargs):
    """Initialize a trajectory writer instance for *filename*.

    :Arguments:
       *filename*
            Output filename of the trajectory; the extension determines the
            format.
       *n_atoms*
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
    return Writer(filename, n_atoms=n_atoms, **kwargs)
