# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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

from ..core._get_readers import get_reader_for
from ..core.topology import Topology
from .base import TopologyReaderBase


class MinimalParser(TopologyReaderBase):
    """Produces a minimal topology from only the number of atoms.

    This requires that the number of atoms be given in one of two ways:
     - The number of atoms can be given as the 'n_atoms' keyword argument.
     - If this is not given, then a Reader object for the filename will be
       created and the `parse_n_atoms` method on this Reader will be called,
       (requiring that the Reader has this capability).

    This requires that the coordinate format has
    """

    format = "MINIMAL"

    def parse(self, **kwargs):
        """Return the minimal *Topology* object"""
        try:
            n_atoms = kwargs["n_atoms"]
        except KeyError:
            reader = get_reader_for(self.filename)
            n_atoms = reader.parse_n_atoms(self.filename, **kwargs)

        return Topology(n_atoms, 1, 1)
