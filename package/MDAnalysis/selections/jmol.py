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
Jmol selections
=================

Write :class:`MDAnalysis.core.AtomGroup.AtomGroup` selection to a `str` file
that defines a `Jmol selection`_. To be used in Jmol_ like this::

  script macro.spt
  select ~selection

The selection is named *mdanalysis001*.TODO

.. autoclass:: SelectionWriter
   :inherited-members:

.. _Jmol: http://wiki.jmol.org/index.php/Main_Page
.. _Jmol selection: http://chemapps.stolaf.edu/jmol/docs/#define
"""
from __future__ import absolute_import

from . import base

class SelectionWriter(base.SelectionWriter):
    format = "Jmol"
    ext = "spt"
    default_numterms = None
    commentfmt = '#'

    def _translate(self, atoms, **kwargs):
        # Jmol indexing is 0 based when using atom bitsets
        def _index(atom):
            return str(atom.index)

        return base.join(atoms, ' ', _index)

    def _write_head(self, out, **kwargs):
        out.write("@~{name!s} ({{".format(**kwargs))

    def _write_tail(self, out, **kwargs):
        out.write("});")
