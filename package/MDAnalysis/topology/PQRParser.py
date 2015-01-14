# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2014 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see AUTHORS for the full list)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
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
"""

import MDAnalysis.coordinates.PQR
from MDAnalysis.topology.core import guess_atom_type, guess_atom_mass

class PQRParser(object):
    """Parse atom information from PQR file *filename*.

    Only reads the list of atoms. Reads the charges and radii from the
    PQR file and populates the
    :attr:`MDAnalysis.core.AtomGroup.Atom.charge` and
    :attr:`MDAnalysis.core.AtomGroup.Atom.radius` attribute.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and
                 :class:`~MDAnalysis.coordinates.PQR.PQRReader` is used to read
                 the PQR file.

    .. versionchanged:: 0.8.2
       Read chainID from a PQR file and use it as segid (before we always used
       'SYSTEM' as the new segid).
    """

    def __init__(self, filename, guess_bonds_mode=False):
        self.PQRReader = MDAnalysis.coordinates.PQR.PQRReader
        self.filename = filename
        self.guess_bonds_mode = guess_bonds_mode

    def parse(self):
        self.structure = {}
        pqr =  self.PQRReader(self.filename)

        self.__parseatoms_(pqr)
        # TODO: reconstruct bonds from CONECT or guess from distance search
        #       (e.g. like VMD)
        return self.structure

    def __parseatoms_(self, pqr):
        from MDAnalysis.core.AtomGroup import Atom
        attr = "_atoms"  # name of the atoms section
        atoms = []       # list of Atom objects

        # translate list of atoms to MDAnalysis Atom.
        for iatom,atom in enumerate(pqr._atoms):
            atomname = atom.name
            atomtype = guess_atom_type(atomname)
            resname = atom.resName
            resid = atom.resSeq
            chain = atom.chainID.strip()
            segid = atom.segID.strip() or chain or "SYSTEM"  # no empty segids (or Universe throws IndexError)
            mass = guess_atom_mass(atomname)
            charge = atom.charge
            radius = atom.radius

            atoms.append(Atom(iatom,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge),
                              radius=radius))
        self.structure[attr] = atoms

# function to keep compatible with the current API; should be cleaned up...
def parse(filename):
    """Parse atom information from PQR file *filename*.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                 :class:`MDAnalysis.coordinates.PQR.PQRReader`.

    """
    return PQRParser(filename).parse()

