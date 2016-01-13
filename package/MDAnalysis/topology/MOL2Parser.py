# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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
MOL2 file format --- :mod:`MDAnalysis.coordinates.MOL2`
========================================================

Classes to read Tripos_ molecule structure format (MOL2_) coordinate
and topology files. Used by the DOCK_ docking code.

.. _MOL2: http://www.tripos.com/data/support/mol2.pdf
.. _Tripos: http://www.tripos.com/
.. _DOCK: http://dock.compbio.ucsf.edu/


Classes
-------

.. autoclass:: MOL2Parser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

import os

from ..lib.util import openany
from ..core.AtomGroup import Atom
from .core import guess_atom_type, guess_atom_mass, guess_atom_charge
from .base import TopologyReader


class MOL2Parser(TopologyReader):
    """Read topology from a Tripos_ MOL2_ file.

    .. versionchanged:: 0.9
       Now subclasses TopologyReader
    """
    format = 'MOL2'

    def parse(self, filename=None):
        """Parse MOL2 file *filename* and return the dict `structure`.

        Only reads the list of atoms.

        :Returns: MDAnalysis internal *structure* dict

        .. SeeAlso:: The *structure* dict is defined in
                     :func:`MDAnalysis.topology.PSFParser.PSFParser`.

        """
        if not filename:
            filename = self.filename

        blocks = []

        with openany(filename) as f:
            for i, line in enumerate(f):
                # found new molecules
                if "@<TRIPOS>MOLECULE" in line:
                    if len(blocks):
                        break
                    blocks.append({"start_line": i, "lines": []})
                blocks[-1]["lines"].append(line)

        if not len(blocks):
            raise ValueError("The mol2 file '{0}' needs to have at least one"
                             " @<TRIPOS>MOLECULE block".format(filename))
        block = blocks[0]

        sections = {}
        cursor = None

        for line in block["lines"]:
            if "@<TRIPOS>" in line:
                cursor = line.split("@<TRIPOS>")[1].strip().lower()
                sections[cursor] = []
                continue
            elif line.startswith("#") or line == "\n":
                continue
            sections[cursor].append(line)

        atom_lines, bond_lines = sections["atom"], sections["bond"]

        if not len(atom_lines):
            raise ValueError("The mol2 block ({0}:{1}) has no atoms".format(
                os.path.basename(filename), block["start_line"]))
        if not len(bond_lines):
            raise ValueError("The mol2 block ({0}:{1}) has no bonds".format(
                os.path.basename(filename), block["start_line"]))

        atoms = []
        for a in atom_lines:
            aid, name, x, y, z, atom_type, resid, resname, charge = a.split()
            aid = int(aid) - 1
            #x, y, z = float(x), float(y), float(z)
            resid = int(resid)
            charge = float(charge)
            element = guess_atom_type(name)
            mass = guess_atom_mass(element)
            # atom type is sybl atom type
            atoms.append(Atom(aid, name, atom_type, resname,
                              resid, "X", mass, charge, universe=self._u))
            #guess_atom_type(a.split()[1]

        bonds = []
        bondorder = {}
        for b in bond_lines:
            # bond_type can be: 1, 2, am, ar
            bid, a0, a1, bond_type = b.split()
            a0, a1 = int(a0) - 1, int(a1) - 1
            bond = tuple(sorted([a0, a1]))
            bondorder[bond] = bond_type
            bonds.append(bond)
        structure = {"atoms": atoms,
                     "bonds": bonds,
                     "bondorder": bondorder}
        return structure
