
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
PDB Topology Parser
=========================================================================

This topology parser uses a standard PDB file to build a minimum
internal structure representation (list of atoms).

The topology reader reads a PDB file line by line and ignores atom
numbers but only reads residue numbers up to 9,999 correctly. If you
have systems containing at least 10,000 residues then you need to use
a different file format (e.g. the "extended" PDB, *XPDB* format, see
:mod:`~MDAnalysis.topology.ExtendedPDBParser`) that can handle residue
numbers up to 99,999.

.. Note::

   The parser processes atoms and their names. Masses are guessed and set to 0
   if unknown. Partial charges are not set.

.. SeeAlso::

   * :mod:`MDAnalysis.topology.ExtendedPDBParser`
   * :class:`MDAnalysis.coordinates.PDB.PDBReader`
   * :class:`MDAnalysis.core.AtomGroup.Universe`

Classes
-------

.. autoclass:: PDBParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import, print_function

import numpy as np
import warnings

from ..core.AtomGroup import Atom
from .core import get_atom_mass, guess_atom_element
from ..lib.util import openany
from .base import TopologyReader


class PDBParser(TopologyReader):
    """Parser that obtains a list of atoms from a standard PDB file.

    See Also
    --------
    :class:`MDAnalysis.coordinates.PDB.PDBReader`

    .. versionadded:: 0.8
    """
    format = ['PDB','ENT']

    def parse(self):
        """Parse atom information from PDB file *filename*.

        Returns
        -------
        MDAnalysis internal *structure* dict

        See Also
        --------
        The *structure* dict is defined in `MDAnalysis.topology` and the file
        is read with :class:`MDAnalysis.coordinates.PDB.PDBReader`.

        """
        structure = {}

        atoms = self._parseatoms()
        structure['atoms'] = atoms

        bonds = self._parsebonds(atoms)
        structure['bonds'] = bonds

        return structure

    def _parseatoms(self):
        iatom = 0
        atoms = []

        with openany(self.filename, 'rt') as f:
            resid_prev = 0  # resid looping hack
            for i, line in enumerate(f):
                line = line.strip()  # Remove extra spaces
                if not line:  # Skip line if empty
                    continue
                record = line[:6].strip()

                if record.startswith('END'):
                    break
                elif line[:6] in ('ATOM  ', 'HETATM'):
                    try:
                        serial = int(line[6:11])
                    except ValueError:
                        # serial can become '***' when they get too high
                        self._wrapped_serials = True
                        serial = None
                    name = line[12:16].strip()
                    altLoc = line[16:17].strip()
                    resName = line[17:21].strip()
                    # empty chainID is a single space ' '!
                    chainID = line[21:22].strip()
                    if self.format == "XPDB":  # fugly but keeps code DRY
                        # extended non-standard format used by VMD
                        resSeq = int(line[22:27])
                        resid = resSeq
                    else:
                        resSeq = int(line[22:26])
                        resid = resSeq

                        while resid - resid_prev < -5000:
                            resid += 10000
                        resid_prev = resid
                        # insertCode = _c(27, 27, str)  # not used
                        # occupancy = float(line[54:60])
                    try:
                        tempFactor = float(line[60:66])
                    except ValueError:
                        tempFactor = 0.0
                    segID = line[66:76].strip()
                    element = line[76:78].strip()

                    segid = segID.strip() or chainID.strip() or "SYSTEM"

                    elem = guess_atom_element(name)

                    atomtype = element or elem
                    mass = get_atom_mass(elem)
                    # charge = guess_atom_charge(name)
                    charge = 0.0

                    atom = Atom(iatom, name, atomtype, resName, resid,
                                segid, mass, charge,
                                bfactor=tempFactor, serial=serial,
                                altLoc=altLoc, universe=self._u,
                                resnum=resSeq)
                    iatom += 1
                    atoms.append(atom)

        return atoms

    def _parsebonds(self, atoms):
        # Could optimise this by saving lines in the main loop
        # then doing post processing after all Atoms have been read
        # ie do one pass through the file only
        # Problem is that in multiframe PDB, the CONECT is at end of file,
        # so the "break" call happens before bonds are reached.

        # If the serials wrapped, this won't work
        if hasattr(self, '_wrapped_serials'):
            warnings.warn("Invalid atom serials were present, bonds will not"
                          " be parsed")
            return tuple([])

        # Mapping between the atom array indicies a.index and atom ids
        # (serial) in the original PDB file
        mapping = dict((a.serial, a.index) for a in atoms)

        bonds = set()
        with openany(self.filename, "rt") as f:
            lines = (line for line in f if line[:6] == "CONECT")
            for line in lines:
                atom, atoms = _parse_conect(line.strip())
                for a in atoms:
                    bond = tuple([mapping[atom], mapping[a]])
                    bonds.add(bond)

        bonds = tuple(bonds)

        return bonds


def _parse_conect(conect):
    """parse a CONECT record from pdbs

    Parameters
    ----------
    conect : str
        white space striped CONECT record

    Returns
    -------
    atom_id : int
        atom index of bond
    bonds : set
        atom ids of bonded atoms

    Raises
    ------
    RuntimeError
        Raised if ``conect`` is not a valid CONECT record
    """
    atom_id = np.int(conect[6:11])
    n_bond_atoms = len(conect[11:]) // 5
    if len(conect[11:]) % n_bond_atoms != 0:
        raise RuntimeError("Bond atoms aren't aligned proberly for CONECT "
                           "record: {}".format(conect))
    bond_atoms = (int(conect[11 + i * 5: 16 + i * 5]) for i in
                  range(n_bond_atoms))
    return atom_id, bond_atoms
