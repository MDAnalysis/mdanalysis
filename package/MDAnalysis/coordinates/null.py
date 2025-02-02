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
Null output --- :mod:`MDAnalysis.coordinates.null`
==================================================

The :class:`NullWriter` provides a Writer instance that behaves like
any other writer but effectively ignores its input and does not write
any output files. This is similar to writing to ``/dev/null``.

This class exists to allow developers writing generic code and tests.

Classes
-------

.. autoclass:: NullWriter
   :members:

"""
from . import base


class NullWriter(base.WriterBase):
    """A trajectory Writer that does not do anything.

    The NullWriter is the equivalent to ``/dev/null``: it behaves like
    a Writer but ignores all input. It can be used in order to
    suppress output.
    """

    format = "NULL"
    multiframe = True
    units = {"time": "ps", "length": "Angstrom"}

    def __init__(self, filename, **kwargs):
        pass

    def _write_next_frame(self, obj):
        try:
            atoms = obj.atoms
        except AttributeError:
            errmsg = "Input obj is neither an AtomGroup or Universe"
            raise TypeError(errmsg) from None
