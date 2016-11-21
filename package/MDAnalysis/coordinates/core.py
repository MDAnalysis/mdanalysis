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
from __future__ import absolute_import

import six

from . import (
    _READERS,
    _SINGLEFRAME_WRITERS,
    _MULTIFRAME_WRITERS,
)

from ..lib import util
from ..lib.util import get_reader_for, get_writer_for
from ..lib.mdamath import triclinic_box, triclinic_vectors, box_volume


def reader(filename, **kwargs):
    """Provide a trajectory reader instance for *filename*.

    This function guesses the file format from the extension of *filename* and
    it will throw a :exc:`TypeError` if the extension is not recognized.

    In most cases, no special keyword arguments are necessary.

    All other keywords are passed on to the underlying Reader classes; see
    their documentation for details.

    Parameters
    ----------
    filename : str or tuple
        filename (or tuple of filenames) of the input coordinate file
    kwargs
        Keyword arguments for the selected Reader class.

    Returns
    -------
    A Reader object

    .. SeeAlso:: For trajectory formats: :class:`~DCD.DCDReader`,
       :class:`~XTC.XTCReader`, :class:`~TRR.TRRReader`,
       :class:`~XYZ.XYZReader`.  For single frame formats:
       :class:`~CRD.CRDReader`, and
       :class:`~PDB.PDBReader`, :class:`~GRO.GROReader`,

    .. deprecated:: 0.15.0
    The "permissive" flag is not used anymore (and effectively
    defaults to True); it will be completely removed in 0.16.0.
    """
    if isinstance(filename, tuple):
        Reader = get_reader_for(filename[0],
                                format=filename[1])
        return Reader(filename[0], **kwargs)
    else:
        Reader = get_reader_for(filename)
        return Reader(filename, **kwargs)


def writer(filename, n_atoms=None, **kwargs):
    """Initialize a trajectory writer instance for *filename*.

    Parameters:
    -----------
    filename : str
        Output filename of the trajectory; the extension determines the
        format.
    n_atoms : int, optional
        The number of atoms in the output trajectory; can be ommitted
        for single-frame writers.
    multiframe : bool, optional
        ``True``: write a trajectory with multiple frames; ``False``
        only write a single frame snapshot; ``None`` first try to get
        a multiframe writer and then fall back to single frame [``None``]
    kwargs : optional
        Keyword arguments for the writer; all trajectory Writers accept
        at least
            *start*
                starting time [0]
            *step*
                step size in frames [1]
            *dt*
                length of time between two frames, in ps [1.0]
       Some readers accept additional arguments, which need to be looked
       up in the documentation of the reader.

    Returns
    -------
    A Writer object

    .. SeeAlso:: :class:`~MDAnalysis.coordinates.DCD.DCDWriter` for DCD
                 trajectories or :class:`~MDAnalysis.coordinates.XTC.XTCWriter`
                 and :class:`~MDAnalysis.coordinates.TRR.TRRWriter` for Gromacs.

    .. versionchanged:: 0.7.6
       Added *multiframe* keyword. See also :func:`get_writer_for`.
    """
    Writer = get_writer_for(filename, format=kwargs.pop('format', None),
                            multiframe=kwargs.pop('multiframe', None))
    return Writer(filename, n_atoms=n_atoms, **kwargs)
