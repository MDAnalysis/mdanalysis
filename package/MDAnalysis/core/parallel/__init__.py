# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2012 Naveen Michaud-Agrawal,
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
Parallel algorithms for MDAnalysis
==================================

The :mod:`MDAnalysis.core.parallel` module contains implementations of
standard functions that can make use of parallelization available on
modern multi-core processors.

.. Warning::

   Using parallel code is under active development in MDAnalysis and
   it is possible that the parallel code has some bugs or
   incompatibilities or less features than the serial code.

.. automodule:: MDAnalysis.core.parallel.distances

"""

import distances

