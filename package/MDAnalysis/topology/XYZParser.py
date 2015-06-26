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
   :inherited-members:

"""
from __future__ import absolute_import

from ..core.AtomGroup import Atom
from ..lib.util import openany
from .core import get_atom_mass, guess_atom_charge, guess_atom_element
from .base import TopologyReader


class XYZParser(TopologyReader):
    """Parse a list of atoms from an XYZ file.

    .. versionadded:: 0.9.1
    """

    def parse(self):
        """Read the file and return the structure.

        :Returns: MDAnalysis internal *structure* dict.
        """
        with openany(self.filename, 'r') as inf:
            natoms = int(inf.readline().strip())
            inf.readline()

            segid = "SYSTEM"
            resid = 1
            resname = "SYSTEM"

            atoms = []
            # Can't infinitely read as XYZ files can be multiframe
            for i in range(natoms):
                name = inf.readline().split()[0]

                elem = guess_atom_element(name)
                mass = get_atom_mass(elem)
                charge = guess_atom_charge(name)

                at = Atom(i, name, elem, resname, resid,
                          segid, mass, charge, universe=self._u)

                atoms.append(at)

        struc = {"atoms": atoms}

        return struc
