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
"""TXYZ topology parser
====================

Tinker_ topology parser: reads information from .txyz and .arc files.
Atom types are read from column 6, while bond connectivity is read from column 7
onwards.

.. _Tinker: https://dasher.wustl.edu/tinker/

See Also
--------
MDAnalysis.coordinates.TXYZ : further documentation on the Tinker format


Classes
-------

.. autoclass:: TXYZParser
   :members:
   :inherited-members:

"""

import itertools
import numpy as np
import warnings

from . import guessers
from .tables import SYMB2Z
from ..lib.util import openany
from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    Atomtypes,
    Bonds,
    Masses,
    Resids,
    Resnums,
    Segids,
    Elements,
)


class TXYZParser(TopologyReaderBase):
    """Parse a list of atoms from a Tinker XYZ file.

    Creates the following attributes:

    - Atomnames
    - Atomtypes
    - Elements (if all atom names are element symbols)

    .. versionadded:: 0.17.0
    .. versionchanged:: 2.4.0
       Adding the `Element` attribute if all names are valid element symbols.
    """
    format = ['TXYZ', 'ARC']

    def parse(self, **kwargs):
        """Read the file and return the structure.

        Returns
        -------
        MDAnalysis Topology object
        """
        with openany(self.filename) as inf:
            #header
            natoms = int(inf.readline().split()[0])

            atomids = np.zeros(natoms, dtype=int)
            names = np.zeros(natoms, dtype=object)
            types = np.zeros(natoms, dtype=object)
            bonds = []
            # Find first atom line, maybe there's box information
            fline = inf.readline()
            try:
                # If a box second value will be a float
                # If an atom, second value will be a string
                float(fline.split()[1])
            except ValueError:
                # If float conversion failed, we have first atom line
                pass
            else:
                # If previous try succeeded it was a box
                # so read another line to find the first atom line
                fline = inf.readline()
            # Can't infinitely read as XYZ files can be multiframe
            for i, line in zip(range(natoms), itertools.chain([fline], inf)):
                line = line.split()
                atomids[i]= line[0]
                names[i] = line[1]
                types[i] = line[5]
                bonded_atoms = line[6:]
                for other_atom in bonded_atoms:
                    other_atom = int(other_atom) - 1
                    if i < other_atom:
                        bonds.append((i, other_atom))

        # Guessing time
        masses = guessers.guess_masses(names)

        attrs = [Atomnames(names),
                 Atomids(atomids),
                 Atomtypes(types),
                 Bonds(tuple(bonds)),
                 Masses(masses, guessed=True),
                 Resids(np.array([1])),
                 Resnums(np.array([1])),
                 Segids(np.array(['SYSTEM'], dtype=object)),
                 ]
        if all(n.capitalize() in SYMB2Z for n in names):
            attrs.append(Elements(np.array(names, dtype=object)))
            
        else:
            warnings.warn("Element information is missing, elements attribute "
                          "will not be populated. If needed these can be "
                          "guessed using MDAnalysis.topology.guessers.")
 
        top = Topology(natoms, 1, 1,
                       attrs=attrs)

        return top
