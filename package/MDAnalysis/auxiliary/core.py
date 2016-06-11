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

def get_auxreader_for(auxdata, format=None):
    """Return the appropriate auxiliary reader class for *auxdata*.

    Paramaters
    ----------
    auxdata
        The auxiliary data (e.g. array) or filename of file containing
        auxiliary data.
    format
        (Optional). Known format of *auxdata*; will be guessed from datatype/
        file extension if not provided.

    Returns
    -------
    An AuxReader object

    Raises
    ------
    ValueError
        If an AuxReader for the guessed format cannot be found.

    """

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
