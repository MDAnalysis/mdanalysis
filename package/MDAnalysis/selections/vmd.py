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
VMD selections
==============

Write :class:`MDAnalysis.core.AtomGroup.AtomGroup` selection to a VMD_ `tcl` file
that defines `atomselect macros`_. To be used in VMD like this::

  source macros.vmd
  set sel [atomselect top mdanalysis001]  # use macro 001

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _atomselect macros: http://www.ks.uiuc.edu/Research/vmd/current/ug/node120.html

.. autoclass:: SelectionWriter
   :inherited-members:
"""
from __future__ import absolute_import

from . import base

class SelectionWriter(base.SelectionWriter):
    format = "VMD"
    ext = "vmd"
    continuation = '\\'
    commentfmt = "# %s"

    def _write_head(self, out, **kwargs):
        out.write(self.comment("MDAnalysis VMD selection"))
        out.write("atomselect macro %(name)s {index " % kwargs)

    def _translate(self, atoms, **kwargs):
        # VMD index is 0-based (as is MDAnalysis)
        return [str(atom.index) for atom in atoms]

    def _write_tail(self, out, **kwargs):
        out.write("}")
