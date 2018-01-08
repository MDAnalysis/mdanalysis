# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Minimal topology parser
=======================

Parses only the number of atoms from a coordinate/trajectory file.
This minimal topology can be used to access coordinate data without the need
for a full-fledged topology file.



Classes
-------

.. autoclass:: MinimalParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

from ..core._get_readers import get_reader_for
from .base import TopologyReaderBase


class MinimalParser(TopologyReaderBase):
    """Parses a coordinate file/trajectory. Reads only the number of atoms.

    Produces a minimal topology from only the number of atoms via
    :meth:`MDAnalysis.core.dummy.make_topology`.
    """
    format = 'MINIMAL'

    def parse(self, **kwargs):
        """Return the minimal *Topology* object for the coordinate file"""
        # Import here to prevent circular imports.
        from ..core.dummy import make_topology

        if 'size' not in kwargs:
            n_atoms = kwargs.get('n_atoms')
            if n_atoms is None:
                reader = get_reader_for(self.filename)
                n_atoms = reader(self.filename, **kwargs).n_atoms
            kwargs['size'] = (n_atoms, 1, 1)
        return make_topology(**kwargs)

