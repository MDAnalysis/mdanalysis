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
PyMOL selections
=================

Write :class:`MDAnalysis.core.AtomGroup.AtomGroup` selection to a
script `pml`_ file that defines PyMOL_ atomselect macros. To be used
in PyMOL like this::

  @macros.pml

The selections should appear in the user interface.

.. _PyMOL: http://www.pymol.org
.. _pml: http://pymol.sourceforge.net/newman/user/S0210start_cmds.html#6_5_1

.. autoclass:: SelectionWriter
   :inherited-members:
"""
from __future__ import absolute_import

from . import base

class SelectionWriter(base.SelectionWriter):
    format = "PyMol"
    ext = "pml"
    continuation = '\\'  # quoted backslash!
    commentfmt = "# %s"
    default_numterms = 6

    def _translate(self, atoms, **kwargs):
        # PyMol index is 1-based
        def _index(atom):
            return "index %d" % (atom.index + 1)

        return base.join(atoms, ' |', _index)

    def _write_head(self, out, **kwargs):
        out.write(self.comment("MDAnalysis PyMol selection"))
        out.write("select %(name)s, " % kwargs + self.continuation + '\n')
