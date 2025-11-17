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
Gromacs selections
==================

Write :class:`MDAnalysis.core.groups.AtomGroup` selection to a `ndx`_ file
that defines a Gromacs_ index group. To be used in Gromacs like this::

  <GROMACS_COMMAND> -n macro.ndx

The index groups are named *mdanalysis001*, *mdanalysis002*, etc.

.. _Gromacs: http://www.gromacs.org
.. _ndx: https://manual.gromacs.org/current/reference-manual/file-formats.html#ndx

.. autoclass:: SelectionWriter
   :inherited-members:
"""
from . import base


class SelectionWriter(base.SelectionWriterBase):
    format = ["Gromacs", "ndx"]
    ext = "ndx"
    default_numterms = 12

    def _translate(self, atoms, **kwargs):
        # Gromacs index is 1-based; MDAnalysis is 0-based
        return [str(atom.index + 1) for atom in atoms]

    def _write_head(self, out, **kwargs):
        out.write("[ {name!s} ]\n".format(**kwargs))
