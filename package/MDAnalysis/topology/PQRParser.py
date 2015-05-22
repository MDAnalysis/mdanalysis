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
PQR topology parser
===================

Read atoms with charges from a PQR_ file (as written by PDB2PQR_). No
connectivity is deduced.

.. SeeAlso:: The file format is described in :mod:`MDAnalysis.coordinates.PQR`.

.. _PQR:     http://www.poissonboltzmann.org/file-formats/biomolecular-structurw/pqr
.. _APBS:    http://www.poissonboltzmann.org/apbs
.. _PDB2PQR: http://www.poissonboltzmann.org/pdb2pqr
.. _PDB:     http://www.rcsb.org/pdb/info.html#File_Formats_and_Standards

Classes
-------

.. autoclass:: PQRParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

from ..core.AtomGroup import Atom
from .core import guess_atom_type, guess_atom_mass
from .base import TopologyReader


class PQRParser(TopologyReader):
    """Parse atom information from PQR file *filename*.

    Only reads the list of atoms. Reads the charges and radii from the
    PQR file and populates the
    :attr:`MDAnalysis.core.AtomGroup.Atom.charge` and
    :attr:`MDAnalysis.core.AtomGroup.Atom.radius` attribute.

    .. versionchanged:: 0.9.0
       Read chainID from a PQR file and use it as segid (before we always used
       'SYSTEM' as the new segid).
    """
    def parse(self):
        """Parse atom information from PQR file *filename*.

        :Returns: MDAnalysis internal *structure* dict

        .. SeeAlso:: The *structure* dict is defined in
                     `MDAnalysis.topology` and
                     :class:`~MDAnalysis.coordinates.PQR.PQRReader` is used to read
                     the PQR file.
        """
        from ..coordinates.PQR import PQRReader
        pqr = PQRReader(self.filename)

        atoms = self._parseatoms(pqr)
        structure = {'atoms': atoms}

        return structure

    def _parseatoms(self, pqr):
        atoms = []

        # translate list of atoms to MDAnalysis Atom.
        for iatom, atom in enumerate(pqr._atoms):
            atomname = atom.name
            atomtype = guess_atom_type(atomname)
            resname = atom.resName
            resid = int(atom.resSeq)
            chain = atom.chainID.strip()
            # no empty segids (or Universe throws IndexError)
            segid = atom.segID.strip() or chain or "SYSTEM"
            mass = guess_atom_mass(atomname)
            charge = float(atom.charge)
            radius = atom.radius

            atoms.append(Atom(iatom, atomname, atomtype, resname, resid,
                              segid, mass, charge, radius=radius, universe=self._u))
        return atoms
