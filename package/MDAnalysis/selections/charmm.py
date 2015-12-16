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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
CHARMM selections
=================

Write :class:`MDAnalysis.core.AtomGroup.AtomGroup` selection to a `str` file
that defines a `CHARMM selection`_. To be used in CHARMM_ like this::

  stream macro.str

The selection is named *mdanalysis001*.

.. autoclass:: SelectionWriter
   :inherited-members:

.. _CHARMM: http://www.charmm.org
.. _CHARMM selection: http://www.charmm.org/documentation/c34b1/select.html
"""
from __future__ import absolute_import

from . import base

class SelectionWriter(base.SelectionWriter):
    format = "CHARMM"
    ext = "str"
    continuation = '-'
    commentfmt = "! %s"
    default_numterms = 4  # be conservative because CHARMM only reads 72 columns

    def _translate(self, atoms, **kwargs):
        # CHARMM index is 1-based
        def _index(atom):
            return "BYNUM {0:d}".format((atom.index + 1))

        return base.join(atoms, ' .or.', _index)

    def _write_head(self, out, **kwargs):
        out.write(self.comment("MDAnalysis CHARMM selection"))
        out.write("DEFINE {name!s} SELECT ".format(**kwargs) + self.continuation + '\n')

    def _write_tail(self, out, **kwargs):
        out.write("END")
