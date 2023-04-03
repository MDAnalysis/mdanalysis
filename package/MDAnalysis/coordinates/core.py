# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
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

from ..core._get_readers import get_reader_for, get_writer_for


def reader(filename, format=None, **kwargs):
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
    :class:`~base.Reader`
        A trajectory Reader instance

    See Also
    --------
    :ref:`Supported coordinate formats`


    """
    if isinstance(filename, tuple):
        Reader = get_reader_for(filename[0],
                                format=filename[1])
        filename = filename[0]
    else:
        Reader = get_reader_for(filename, format=format)
    try:
        return Reader(filename, **kwargs)
    except ValueError:
        errmsg = f'Unable to read {filename} with {Reader}.'
        raise TypeError(errmsg) from None


def writer(filename, n_atoms=None, **kwargs):
    """Initialize a trajectory writer instance for *filename*.

    Parameters
    ----------
    filename : str
        Output filename of the trajectory; the extension determines the
        format.
    n_atoms : int (optional)
        The number of atoms in the output trajectory; can be ommitted
        for single-frame writers.
    multiframe : bool (optional)
        ``True``: write a trajectory with multiple frames; ``False``
        only write a single frame snapshot; ``None`` first try to get
        a multiframe writer and then fall back to single frame [``None``]
    kwargs : optional
        Keyword arguments for the writer; all trajectory Writers accept
        ``start``: starting time [0], ``step``: step size in frames [1],
        ``dt``: length of time between two frames, in ps [1.0] Some readers
        accept additional arguments, which need to be looked up in the
        documentation of the reader.

    Returns
    -------
    :class:`~base.Writer`
        A trajectory Writer instance

    See Also
    --------
    :ref:`Supported coordinate formats`


    .. versionchanged:: 0.7.6
       Added `multiframe` keyword. See also :func:`get_writer_for`.

    """
    Writer = get_writer_for(filename, format=kwargs.pop('format', None),
                            multiframe=kwargs.pop('multiframe', None))
    return Writer(filename, n_atoms=n_atoms, **kwargs)
