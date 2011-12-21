# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Status codes and symbols
========================

The module makes all the status symbols available that are used in
:mod:`libxdrfile`. The value of each symbol is an integer (as defined
in ``xdrfile.h``).

The dictionary :data:`xdrfile.errno.errorcode` maps numeric codes to
symbol names.

.. data:: errorcode
.. data:: errorsymbols
"""

import libxdrfile

#: List of all error symbols ``exdr*`` extracted from :mod:`libxdrfile`.
errorsymbols = [k for k in libxdrfile.__dict__.keys() if k[:4] == 'exdr']

#: Dictionary that maps error codes to symbol names.
errorcode = dict(((libxdrfile.__dict__[k], k) for k in errorsymbols))

globals().update(dict((errorcode[n], n) for n in errorcode))

