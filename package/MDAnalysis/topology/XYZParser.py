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
from __future__ import absolute_import

from six.moves import range
import numpy as np

from . import guessers
from ..lib.util import openany
from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    Atomtypes,
    Masses,
    Resids,
    Resnums,
    Segids,
)


class XYZParser(TopologyReaderBase):
    """Parse a list of atoms from an XYZ file.

    Creates the following attributes:
     - Atomnames

    Guesses the following attributes:
     - Atomtypes
     - Masses

    .. versionadded:: 0.9.1
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

        # Guessing time
        atomtypes = guessers.guess_types(names)
        masses = guessers.guess_masses(atomtypes)

        attrs = [Atomnames(names),
                 Atomids(np.arange(natoms) + 1),
                 Atomtypes(atomtypes, guessed=True),
                 Masses(masses, guessed=True),
                 Resids(np.array([1])),
                 Resnums(np.array([1])),
                 Segids(np.array(['SYSTEM'], dtype=object)),
                 ]

        top = Topology(natoms, 1, 1,
                       attrs=attrs)

        return top
