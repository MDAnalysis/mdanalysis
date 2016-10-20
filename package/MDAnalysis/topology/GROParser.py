# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
GRO topology parser
===================

Read a list of atoms from a GROMOS/Gromacs GRO coordinate file to
build a basic topology.

Atom types and masses are guessed.

.. SeeAlso:: :mod:`MDAnalysis.coordinates.GRO`

Classes
-------

.. autoclass:: GROParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

import numpy as np
from six.moves import range

from ..lib.util import openany
from ..core.topologyattrs import (
    Atomnames,
    Atomtypes,
    Atomids,
    Elements,
    Masses,
    Resids,
    Resnames,
    Resnums,
    Segids,
)
from ..core.topology import Topology
from .base import TopologyReader, squash_by
from . import guessers


class GROParser(TopologyReader):
    """Reads a Gromacs GRO file

    Reads the following attributes:
      - resids
      - resnames
      - atomids
      - atomnames

    Guesses the following attributes
      - elements
      - masses
    """
    format = 'GRO'

    def parse(self):
        """Return the *Topology* object for this file"""
        # Gro has the following columns
        # resid, resname, name, index, (x,y,z)
        with openany(self.filename, 'rt') as inf:
            inf.readline()
            n_atoms = int(inf.readline())

            # Allocate shizznizz
            resids = np.zeros(n_atoms, dtype=np.int32)
            resnames = np.zeros(n_atoms, dtype=object)
            names = np.zeros(n_atoms, dtype=object)
            indices = np.zeros(n_atoms, dtype=np.int32)

            for i in range(n_atoms):
                line = inf.readline()
                try:
                    resids[i] = int(line[:5])
                    resnames[i] = line[5:10].strip()
                    names[i] = line[10:15].strip()
                    indices[i] = int(line[15:20])
                except (ValueError, TypeError):
                    raise IOError(
                        "Couldn't read the following line of the .gro file:\n"
                        "{0}".format(line))
        # Check all lines had names
        if not np.all(names):
            missing = np.where(names == '')
            raise IOError("Missing atom name on line: {0}"
                          "".format(missing[0][0] + 3))  # 2 header, 1 based

        # Guess types and masses
        elements = guessers.guess_types(names)
        masses = guessers.guess_masses(elements)

        residx, new_resids, (new_resnames,) = squash_by(resids, resnames)

        # new_resids is len(residues)
        # so resindex 0 has resid new_resids[0]
        attrs = [
            Atomnames(names),
            Atomids(indices),
            Elements(elements, guessed=True),
            Resids(new_resids),
            Resnums(new_resids.copy()),
            Resnames(new_resnames),
            Masses(masses, guessed=True),
            Segids(np.array(['SYSTEM'], dtype=object))
        ]

        top = Topology(n_atoms=n_atoms, n_res=len(new_resids), n_seg=1,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=None)

        return top
