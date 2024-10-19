# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
XYZ Topology Parser
===================

.. versionadded:: 0.9.1

Reads an xyz file and pulls the atom information from it.  Because
xyz only has atom name information, all information about residues
and segments won't be populated.

Classes
-------

.. autoclass:: XYZParser
   :members:

"""
import numpy as np

from ..lib.util import openany
from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    Resids,
    Resnums,
    Segids,
    Elements,
)


class XYZParser(TopologyReaderBase):
    """Parse a list of atoms from an XYZ file.

    Creates the following attributes:
     - Atomnames


    .. versionadded:: 0.9.1

    .. versionchanged: 1.0.0
       Store elements attribute, based on XYZ atom names
     .. versionchanged:: 2.8.0
        Removed type and mass guessing (attributes guessing takes place now
        through universe.guess_TopologyAttrs() API).

    """
    format = 'XYZ'

    def parse(self, **kwargs):
        """Read the file and return the structure.

        Returns
        -------
        MDAnalysis Topology object
        """
        with openany(self.filename) as inf:
            natoms = int(inf.readline().strip())
            inf.readline()

            names = np.zeros(natoms, dtype=object)

            # Can't infinitely read as XYZ files can be multiframe
            for i in range(natoms):
                name = inf.readline().split()[0]
                names[i] = name


        attrs = [Atomnames(names),
                 Atomids(np.arange(natoms) + 1),
                 Resids(np.array([1])),
                 Resnums(np.array([1])),
                 Segids(np.array(['SYSTEM'], dtype=object)),
                 Elements(names)]

        top = Topology(natoms, 1, 1,
                       attrs=attrs)

        return top
