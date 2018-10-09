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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Common functions for auxiliary reading --- :mod:`MDAnalysis.auxiliary.core`
===========================================================================

.. autofunction:: get_auxreader_for
.. autofunction:: auxreader
"""
from __future__ import absolute_import

from six import string_types

from . import _AUXREADERS
from ..lib import util

def get_auxreader_for(auxdata=None, format=None):
    """Return the appropriate auxiliary reader class for *auxdata*/*format*.

    If *format* is provided, will attempt to find an AuxReader corresponding
    to that format. If *auxdata* is provided, the format will first be guessed.

    Parameters
    ----------
    auxdata
        (Optional) The auxiliary data (e.g. filename of file containing
        auxiliary data).
    format
        (Optional). Known format of *auxdata*.

    Returns
    -------
    :class:`~MDAnalysis.auxiliary.base.AuxReader`
        AuxReader class corresponding to the supplied/guessed format.

    Raises
    ------
    ValueError
        If an AuxReader for the format (provided or guessed from *auxdata*)
        cannot be found.

    """
    if not auxdata and not format:
        raise ValueError('Must provide either auxdata or format')

    if format is None:
        if isinstance(auxdata, string_types):
            ## assume it's a filename?
            format = util.guess_format(auxdata)
        else:
            ## TBA if add non-file-format readers
            pass
        format = format.upper()
        try:
            return _AUXREADERS[format]
        except KeyError:
            raise ValueError("Unknown auxiliary data format for auxdata: "
                             "{0}".format(auxdata))
    else:
        try:
            return _AUXREADERS[format]
        except KeyError:
            raise ValueError("Unknown auxiliary data format {0}".format(format))

def auxreader(auxdata, format=None, **kwargs):
    """ Return an auxiliary reader instance for *auxdata*.

    An appropriate reader class is first obtained using 
    :func:`get_auxreader_for`, and an auxiliary reader instance for *auxdata*
    then created and returned.

    Parameters
    ----------
    auxdata
        Auxiliary data (e.g. filename of file containing auxiliary data).
    format
        (Optional). The format of *auxdata*, if known.
    **kwargs
         Additional AuxReader options.

    Returns
    -------
    :class:`~MDAnalysis.auxiliary.base.AuxReader` instance
        Appropriate auxiliary reader instance for *auxdata*.
    """
    reader = get_auxreader_for(auxdata, format=format)
    return reader(auxdata, **kwargs)
