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
VMD selections
==============

Write :class:`MDAnalysis.core.groups.AtomGroup` selection to a VMD_ `tcl` file
(for example, "macros.vmd") that defines `atomselect macros`_. To be used in
VMD like this (assuming the default macro name "mdanalysis001"):

.. code-block:: tcl

   source macros.vmd
   set sel [atomselect top mdanalysis001]

In the VMD_ GUI the macro "mdanalysis001" appears in the
:menuselection:`Graphics --> Representations` window under
:guilabel:`Selections: Singlewords`.


.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _atomselect macros: http://www.ks.uiuc.edu/Research/vmd/current/ug/node120.html

.. autoclass:: SelectionWriter
   :inherited-members:

"""
from . import base


class SelectionWriter(base.SelectionWriterBase):
    format = "VMD"
    ext = "vmd"
    continuation = "\\"
    commentfmt = "# %s"

    def _write_head(self, out, **kwargs):
        out.write(self.comment("MDAnalysis VMD selection"))
        out.write("atomselect macro {name!s} {{index ".format(**kwargs))

    def _translate(self, atoms, **kwargs):
        # VMD index is 0-based (as is MDAnalysis)
        return [str(atom.index) for atom in atoms]

    def _write_tail(self, out, **kwargs):
        out.write("}")
