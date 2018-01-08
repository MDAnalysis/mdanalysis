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
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""Functions for fetching Readers

These functions officially live in in topology/core (parsers) and
coordinates/core (all others).  They are declared here to avoid
circular imports.

"""
from __future__ import absolute_import

import copy
import inspect
import mmtf
import numpy as np
from MDAnalysis.lib.util import isstream

from .. import _READERS, _PARSERS, _MULTIFRAME_WRITERS, _SINGLEFRAME_WRITERS
from ..lib import util


def get_reader_for(filename, format=None):
    """Return the appropriate trajectory reader class for `filename`.

    Parameters
    ----------
    filename
        filename of the input trajectory or coordinate file.  Also can
        handle a few special cases, see notes below.
    format : str or :class:`Reader` (optional)
        Define the desired format.  Can be a string to request a given
        Reader.
        If a class is passed, it will be assumed that this is
        a Reader and will be returned.

    Returns
    -------
    :class:`Reader`
        A Reader object


    Raises
    ------
    ValueError
        If no appropriate Reader is found


    Notes
    -----
    There are a number of special cases that can be handled:

    - If `filename` is a numpy array,
      :class:`~MDAnalysis.coordinates.memory.MemoryReader` is returned.
    - If `filename` is an MMTF object,
      :class:`~MDAnalysis.coordinates.MMTF.MMTFReader` is returned.
    - If `filename` is an iterable of filenames,
      :class:`~MDAnalysis.coordinates.chain.ChainReader` is returned.

    Automatic detection is disabled when an explicit `format` is provided,
    unless a list of filenames is given, in which case
    :class:`~MDAnalysis.coordinates.chain.ChainReader` is returned and `format`
    passed to the :class:`~MDAnalysis.coordinates.chain.ChainReader`.

    """
    # check if format is actually a Reader
    if inspect.isclass(format):
        return format

    # ChainReader gets returned even if format is specified
    if not isinstance(filename, np.ndarray) and util.iterable(filename) and not isstream(filename):
        format = 'CHAIN'
    # Only guess if format is not specified
    elif format is None:
        # Checks for specialised formats
        if isinstance(filename, np.ndarray):
            # memoryreader slurps numpy arrays
            format = 'MEMORY'
        elif isinstance(filename, mmtf.MMTFDecoder):
            # mmtf slurps mmtf object
            format = 'MMTF'
        else:
            # else let the guessing begin!
            format = util.guess_format(filename)
    format = format.upper()
    try:
        return _READERS[format]
    except KeyError:
        raise ValueError(
            "Unknown coordinate trajectory format '{0}' for '{1}'. The FORMATs \n"
            "           {2}\n"
            "           are implemented in MDAnalysis.\n"
            "           See https://docs.mdanalysis.org/documentation_pages/coordinates/init.html#id1\n"
            "           Use the format keyword to explicitly set the format: 'Universe(...,format=FORMAT)'\n"
            "           For missing formats, raise an issue at "
            "http://issues.mdanalysis.org".format(
                format, filename, _READERS.keys()))


def get_writer_for(filename, format=None, multiframe=None):
    """Return an appropriate trajectory or frame writer class for `filename`.

    The format is determined by the `format` argument or the extension of
    `filename`. If `format` is provided, it takes precedence over the
    extension of `filename`.

    Parameters
    ----------
    filename : str or ``None``
        If no *format* is supplied, then the filename for the trajectory is
        examined for its extension and the Writer is chosen accordingly.
        If ``None`` is provided, then
        :class:`~MDAnalysis.coordinates.null.NullWriter` is selected (and
        all output is discarded silently).
    format : str (optional)
        Explicitly set a format.
    multiframe : bool (optional)
        ``True``: write multiple frames to the trajectory; ``False``: only
        write a single coordinate frame; ``None``: first try trajectory (multi
        frame writers), then the single frame ones. Default is ``None``.

    Returns
    -------
    :class:`Writer`
        A Writer object

    Raises
    ------
    ValueError:
        The format could not be deduced from `filename` or an unexpected value
        was provided for the `multiframe` argument.
    TypeError:
        No writer was found for the required format or the required `filename`
        argument was omitted.


    .. versionchanged:: 0.7.6
       Added `multiframe` keyword; the default ``None`` reflects the previous
       behaviour.

    .. versionchanged:: 0.14.0
       Removed the default value for the `format` argument. Now, the value
       provided with the `format` parameter takes precedence over the extension
       of `filename`. A :exc:`ValueError` is raised if the format cannot be
       deduced from `filename`.

    .. versionchanged:: 0.16.0
       The `filename` argument has been made mandatory.
    """
    if filename is None:
        format = 'NULL'
    elif format is None:
        try:
            root, ext = util.get_ext(filename)
        except (TypeError, AttributeError):
            # An AttributeError is raised if filename cannot
            # be manipulated as a string.
            # A TypeError is raised in py3.6
            # "TypeError: expected str, bytes or os.PathLike object"
            raise ValueError('File format could not be guessed from "{0}"'
                             .format(filename))
        else:
            format = util.check_compressed_format(root, ext)
    format = format.upper()
    if multiframe is None:
        # Multiframe takes priority, else use singleframe
        options = copy.copy(_SINGLEFRAME_WRITERS)  # do copy to avoid changing in place
        options.update(_MULTIFRAME_WRITERS)  # update overwrites existing entries
        errmsg = "No trajectory or frame writer for format '{0}'"
    elif multiframe is True:
        options = _MULTIFRAME_WRITERS
        errmsg = "No trajectory writer for format '{0}'"
    elif multiframe is False:
        options = _SINGLEFRAME_WRITERS
        errmsg = "No single frame writer for format '{0}'"
    else:
        raise ValueError("Unknown value '{0}' for multiframe,"
                         " only True, False, None allowed"
                         "".format(multiframe))

    try:
        return options[format]
    except KeyError:
        raise TypeError(errmsg.format(format))


def get_parser_for(filename, format=None):
    """Return the appropriate topology parser for `filename`.

    Automatic detection is disabled when an explicit `format` is
    provided.

    Parameters
    ----------
    filename : str or mmtf.MMTFDecoder
        name of the topology file; if this is an instance of
        :class:`mmtf.MMTFDecoder` then directly use the MMTF format.
    format : str
        description of the file format

    Raises
    ------
    ValueError
        If no appropriate parser could be found.

    """
    if inspect.isclass(format):
        return format

    # Only guess if format is not provided
    if format is None:
        if isinstance(filename, mmtf.MMTFDecoder):
            format = 'mmtf'
        else:
            format = util.guess_format(filename)
    format = format.upper()
    try:
        return _PARSERS[format]
    except KeyError:
        try:
            rdr = get_reader_for(filename)
        except ValueError:
            raise ValueError(
                "'{0}' isn't a valid topology format, nor a coordinate format\n"
                "   from which a topology can be minimally inferred.\n"
                "   You can use 'Universe(topology, ..., topology_format=FORMAT)'\n"
                "   to explicitly specify the format and\n"
                "   override automatic detection. Known FORMATs are:\n"
                "   {1}\n"
                "   See https://docs.mdanalysis.org/documentation_pages/topology/init.html#supported-topology-formats\n"
                "   For missing formats, raise an issue at \n"
                "   http://issues.mdanalysis.org".format(format, _PARSERS.keys()))
        else:
            return _PARSERS['MINIMAL']
