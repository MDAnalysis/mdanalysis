# Primitive PQR parser
# -*- coding: utf-8 -*-
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

class PQRParseError(Exception):
    pass

def parse(filename):
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
    """
    structure = {}
    pqr =  MDAnalysis.coordinates.PQR.PQRReader(filename)

    __parseatoms_(pqr, structure)
    # TODO: reconstruct bonds from CONECT or guess from distance search
    #       (e.g. like VMD)
    return structure

def __parseatoms_(pqr, structure):
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
        segid = atom.segID.strip() or "SYSTEM"  # no empty segids (or Universe throws IndexError)
        mass = guess_atom_mass(atomname)
        charge = atom.charge
        radius = atom.radius

        atoms.append(Atom(iatom,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge),
                          radius=radius))

    structure[attr] = atoms
