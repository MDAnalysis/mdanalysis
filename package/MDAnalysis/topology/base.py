# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
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
Base topology reader classes --- :mod:`MDAnalysis.topology.base`
================================================================

Derive topology reader classes from the base class in this module.


"""

from ..coordinates.base import IObase


class TopologyReader(IObase):
    """Base class for topology readers

    All topology readers must:
      * Be initialised with a filename
      * Return a struct dict by calling their parse function

    :Raises:
       IOError upon failing to read a topology file
       ValueError upon failing to make sense of the read data

    .. versionadded:: 0.9
    """
    def __init__(self, filename, guess_bonds_mode=False, **kwargs):
        self.filename = filename
        self.guess_bonds_mode = guess_bonds_mode
        self.kwargs = kwargs

    def parse(self):
        raise NotImplementedError("Override this in each subclass")
