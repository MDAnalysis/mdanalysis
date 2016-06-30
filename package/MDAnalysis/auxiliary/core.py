# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
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
Common functions for auxiliary reading --- :mod:`MDAnalysis.auxiliary.core`
===========================================================================

.. autofunction:: get_auxreader_for
"""

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
        (Optional) The auxiliary data (e.g. array) or filename of file containing
        auxiliary data.
    format
        (Optional). Known format of *auxdata*; if not provided, will be guessed 
        from datatype/file extension.

    Returns
    -------
    An AuxReader object

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
            ## arrays etc
            pass
        format = format.upper()
        try:
            return _AUXREADERS[format]
        except KeyError:
            raise ValueError("Unknown auxiliary data format for {0}".format(auxdata))
    else:
        try:
            return _AUXREADERS[format]
        except KeyError:
            raise ValueError("Unknown auxiliary data format {0}".format(format))
