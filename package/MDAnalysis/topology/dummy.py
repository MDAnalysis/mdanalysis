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
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
Dummy Topology parser class --- :mod:`MDAnalysis.topology.dummy`
================================================================

A rather simple Parser which just creates a number of blank
:class:`MDAnalysis.core.AtomGroup.Atom` instances.
These will have no unique properties except their index.
Designed to act as a base for building a system from scratch.


.. versionadded:: 0.13.0
"""
from .base import TopologyReader
from ..core.AtomGroup import Atom


class DummyParser(TopologyReader):
    def parse(self):
        # the 'filename' of this Parser is an int of the number
        # of atoms to pass back.
        # Yay for dynamic typing
        atoms = [Atom(i, 'NAME', 'TYPE', 'RESNAME', 1, 'SEGID',
                      0.0, 0.0, universe=self._u)
                 for i in range(self.filename)]
        return {"atoms": atoms}
