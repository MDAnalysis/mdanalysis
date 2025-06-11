# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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

"""Homogeneous Transformation Matrices and Quaternions --- :mod:`MDAnalysis.lib.transformations`
==============================================================================================

.. deprecated:: 2.10.0
   This module will be removed in 3.0.0

The ``lib.transformations`` module makes the functions from the
`transformations package`_ available.

.. SeeAlso::

   `transformations package`_ by Christoph Gohlke.

   The function :func:`MDAnalysis.lib.transformations.rotaxis` was
   moved to :func:`MDAnalysis.lib.mdamath.rotaxis`.

.. _`transformations package`:
   https://github.com/cgohlke/transformations?tab=readme-ov-file#homogeneous-transformation-matrices-and-quaternions

"""

from transformations import *

from .mdamath import rotaxis

from .util import deprecate

rotaxis = deprecate(
    rotaxis,
    old_name="MDAnalysis.lib.transformations.rotaxis",
    new_name="MDAnalysis.lib.mdamath.rotaxis",
    release="2.10.0",
    remove="3.0.0",
)
