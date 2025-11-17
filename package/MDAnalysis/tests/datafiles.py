# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
Location of data files for the MDAnalysis unit tests
====================================================

Real MD simulation data, used for examples and the unit tests::

  from MDAnalysis.tests.datafiles import *

.. Note::

   The actual files are contained in the separate
   :mod:`MDAnalysisTests` package which must be downloaded from
   http://pypi.python.org/pypi/MDAnalysisTests and installed.
"""

try:
    from MDAnalysisTests.datafiles import *
except ImportError:
    print("*** ERROR ***")
    print("In order to run the MDAnalysis test cases you must install the")
    print("MDAnalysisTestData package (which has been separated from the ")
    print("library code itself since release 0.7.4). Go to ")
    print()
    print("     http://pypi.python.org/pypi/MDAnalysisTests")
    print()
    print("and download and install the `MDAnalysisTests-x.y.z.tar.gz'")
    print("that matches your MDAnalysis release.")
    raise ImportError("MDAnalysisTests package not installed.") from None
