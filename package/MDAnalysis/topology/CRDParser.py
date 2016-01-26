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
CRD topology parser
===================

Read a list of atoms from a CHARMM CARD coordinate file (CRD) to build a basic topology.

Atom types, charges and masses are guessed.

Classes
-------

.. autoclass:: CRDParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

from ..core.AtomGroup import Atom
from ..lib.util import openany, FORTRANReader
from .core import guess_atom_type, guess_atom_mass, guess_atom_charge
from .base import TopologyReader


class CRDParser(TopologyReader):
    """Parse a CHARMM CARD coordinate file for topology information."""
    format = 'CRD'

    def parse(self):
        """Parse CRD file *filename* and return the dict `structure`.

        Only reads the list of atoms.

        :Returns: MDAnalysis internal *structure* dict

        .. SeeAlso:: The *structure* dict is defined in
                     `MDAnalysis.topology`
        """
        extformat = FORTRANReader('2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10')
        stdformat = FORTRANReader('2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5')

        atoms = []
        atom_serial = 0
        with openany(self.filename) as crd:
            for linenum, line in enumerate(crd):
                # reading header
                if line.split()[0] == '*':
                    continue
                elif line.split()[-1] == 'EXT' and bool(int(line.split()[0])) is True:
                    r = extformat
                    continue
                elif line.split()[0] == line.split()[-1] and line.split()[0] != '*':
                    r = stdformat
                    continue
                # anything else should be an atom
                try:
                    serial, TotRes, resName, name, x, y, z, chainID, resSeq, tempFactor = r.read(line)
                except:
                    raise ValueError("Check CRD format at line {0}: {1}"
                                     "".format(linenum, line.rstrip()))

                atomtype = guess_atom_type(name)
                mass = guess_atom_mass(name)
                charge = guess_atom_charge(name)
                atoms.append(Atom(atom_serial, name, atomtype, resName,
                                  TotRes, chainID, mass, charge, universe=self._u))
                atom_serial += 1

        structure = {}
        structure['atoms'] = atoms

        return structure
