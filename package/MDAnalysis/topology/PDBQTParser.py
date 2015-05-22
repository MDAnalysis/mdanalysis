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
PDBQT topology parser
=====================

Use a PDBQT_ file to build a minimum internal structure representation (list of
atoms), including AutoDock_ atom types (stored as :attr:`Atom.type`) and
partial charges (:attr:`Atom.charge`).

* Reads a PDBQT file line by line and does not require sequential atom numbering.
* Multi-model PDBQT files are not supported.

.. Note:: Only reads atoms and their names; connectivity is not
          deduced. Masses are guessed and set to 0 if unknown.

.. SeeAlso:: `MDAnalysis.coordinates.PDBQT`

.. _PDBQT:
   http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
.. _AutoDock:
   http://autodock.scripps.edu/

Classes
-------

.. autoclass:: PDBQTParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

from ..core.AtomGroup import Atom
from .core import guess_atom_type, guess_atom_mass, guess_atom_charge
from .base import TopologyReader


class PDBQTParser(TopologyReader):
    """Read topology from a PDBQT file."""

    def parse(self):
        """Parse atom information from PDBQT file *filename*.

        :Returns: MDAnalysis internal *structure* dict

        .. SeeAlso:: The *structure* dict is defined in
                     :func:`MDAnalysis.topology.PSFParser.PSFParser` and the file
                     is read with
                     :class:`MDAnalysis.coordinates.PDBQT.PDBQTReader`.
        """
        atoms = self._parseatoms()

        structure = {'atoms': atoms}

        return structure

    def _parseatoms(self):
        from ..coordinates.PDBQT import PDBQTReader
        pdb = PDBQTReader(self.filename)

        atoms = []
        # translate list of atoms to MDAnalysis Atom.
        for iatom, atom in enumerate(pdb._atoms):
            atomname = atom.name
            atomtype = atom.type        # always set in PDBQT
            resname = atom.resName
            resid = int(atom.resSeq)
            chain = atom.chainID.strip()
            segid = chain or "SYSTEM"   # no empty segids (or Universe throws IndexError)
            mass = guess_atom_mass(atomname)
            charge = float(atom.partialCharge)  # always set in PDBQT
            bfactor = atom.tempFactor
            # occupancy = atom.occupancy
            atoms.append(Atom(iatom, atomname, atomtype, resname, resid, segid,
                              mass, charge, bfactor=bfactor, universe=self._u))
        return atoms
