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
Null output --- :mod:`MDAnalysis.coordinates.null`
==================================================

The :class:`NullWriter` provides a Writer instance that behaves like
any other writer but effectively ignores its input and does not write
any output files. This is similar to writing to ``/dev/null``.

This class exists to allow developers writing generic code and tests.

"""
from __future__ import absolute_import

from . import base


class NullWriter(base.Writer):
    """A trajectory Writer that does not do anything.

    The NullWriter is the equivalent to ``/dev/null``: it behaves like
    a Writer but ignores all input. It can be used in order to
    suppress output.

    """
    format = 'NULL'
    units = {'time': 'ps', 'length': 'Angstrom'}

    def __init__(self, filename, **kwargs):
        pass

    def write_next_timestep(self, ts=None):
        pass
