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

"""
Custom exceptions and warnings --- :mod:`MDAnalysis.exceptions`
===============================================================

"""

class SelectionError(Exception):
    """Raised when a atom selection failed."""


class FinishTimeException(Exception):
    """For Issue 188."""


class NoDataError(ValueError):
    """Raised when empty input is not allowed or required data are missing."""


class ApplicationError(OSError):
    """Raised when an external application failed.

    The error code is specific for the application.

    .. versionadded:: 0.7.7
    """


class SelectionWarning(Warning):
    """Warning indicating a possible problem with a selection."""


class MissingDataWarning(Warning):
    """Warning indicating is that required data are missing."""


class ConversionWarning(Warning):
    """Warning indicating a problem with converting between units."""


class FileFormatWarning(Warning):
    """Warning indicating possible problems with a file format."""


class StreamWarning(Warning):
    """Warning indicating a possible problem with a stream.

    :exc:`StreamWarning` is used when streams are substituted for simple access
    by filename (see in particular
    :class:`~MDAnalysis.lib.util.NamedStream`). This does not work everywhere
    in MDAnalysis (yet).
    """
