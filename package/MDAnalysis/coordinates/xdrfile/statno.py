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
Status codes and symbols --- :mod:`~MDAnalysis.coordinates.xdrfile.statno`
==========================================================================

The module makes all the status symbols available that are used in
:mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2` (and which are
described there). The value of each symbol is an integer (as defined
in ``xdrfile.h``).

The dictionary :data:`xdrfile.errno.errorcode` maps numeric codes to
symbol names.

.. data:: ERRORCODE
.. data:: ERRORSYMBOLS

"""

from . import libxdrfile2

#: List of all error symbols ``exdr*`` extracted from :mod:`libxdrfile2`.
ERRORSYMBOLS = [k for k in libxdrfile2.__dict__.keys() if k[:4] == 'exdr']

#: Dictionary that maps error codes to symbol names.
ERRORCODE = dict(((libxdrfile2.__dict__[k], k) for k in ERRORSYMBOLS))

globals().update(dict((ERRORCODE[n], n) for n in ERRORCODE))
