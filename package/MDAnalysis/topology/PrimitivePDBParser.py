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
Primitive PDB topology parser
=============================

This topology parser uses a standard PDB file to build a minimum
internal structure representation (list of atoms).

The topology reader reads a PDB file line by line and ignores atom
numbers but only reads residue numbers up to 9,999 correctly. If you
have systems containing at least 10,000 residues then you need to use
a different file format (e.g. the "extended" PDB, *XPDB* format, see
:mod:`~MDAnalysis.topology.ExtendedPDBParser`) that can handle residue
numbers up to 99,999.

.. SeeAlso::

   * :mod:`MDAnalysis.topology.ExtendedPDBParser`
   * :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`
   * :class:`MDAnalysis.core.AtomGroup.Universe`

Classes
-------

.. autoclass:: PrimitivePDBParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import, print_function

import numpy as np
import warnings
import numpy as np

from ..lib.util import openany
from .base import TopologyReader, squash_by
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomtypes,
    Atomids,
    AltLocs,
    Tempfactors,
    ICodes,
    ChainIDs,
    Occupancies,
    Resids,
    Resnames,
    Resnums,
    Segids,
    Bonds,
    AtomAttr,
)


class PrimitivePDBParser(TopologyReader):
    """Parser that obtains a list of atoms from a standard PDB file.

    Creates the following per atom Attributes:
    - serials [optional]
    ` - names
    - types
    - altLocs
    - insertion_codes (non XPDB)
    - chainids
    - tempfactors [optional]
    - occupancies [optional]
    per residue Attributes:
    - resids
    - resnums (non XPDB)
    - resnames
    per segment Attributes:
    - segids

    .. seealso:: :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`

    See Also
    --------
    :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`

    .. versionadded:: 0.8
    """
    format = 'Permissive_PDB'

    def parse(self):
        """Parse atom information from PDB file

        Returns
        -------
        MDAnalysis Topology object
        """
        top = self._parseatoms()

        try:
            bonds = self._parsebonds(top.ids.values)
        except AttributeError:
            warnings.warn("Invalid atom serials were present, "
                          "bonds will not be parsed")
        else:
            top.add_TopologyAttr(bonds)

        return top

    def _parseatoms(self):
        """Create the initial Topology object"""
        resid_prev = 0  # resid looping hack

        serials = []
        names = []
        altlocs = []
        chainids = []
        icodes = []
        tempfactors = []
        occupancies = []
        atomtypes = []

        resids = []
        resnums = []
        resnames = []

        segids = []

        with openany(self.filename) as f:
            for line in f:
                line = line.strip()  # Remove extra spaces
                if not line:  # Skip line if empty
                    continue
                if line.startswith('END'):
                    break
                if not line.startswith(('ATOM', 'HETATM')):
                    continue

                try:
                    serial = int(line[6:11])
                except ValueError:
                    # serial can become '***' when they get too high
                    self._wrapped_serials = True
                    serial = None
                serials.append(serial)

                names.append(line[12:16].strip())
                altlocs.append(line[16:17].strip())
                resnames.append(line[17:21].strip())
                # empty chainID is a single space ' '!
                chainid = line[21:22].strip()
                if not chainid:
                    chainid = None
                chainids.append(chainid)

                if self.format == "XPDB":  # fugly but keeps code DRY
                    # extended non-standard format used by VMD
                    resids.append(int(line[22:27]))
                else:
                    resid = int(line[22:26])
                    icodes.append(line[26:27].strip())

                    resnums.append(resid)

                    while resid - resid_prev < -5000:
                        resid += 10000
                    resid_prev = resid
                    resids.append(resid)

                try:
                    occupancy = float(line[54:60])
                except ValueError:
                    occupancy = None
                occupancies.append(occupancy)

                try:
                    # bfactor
                    tempFactor = float(line[60:66])
                except ValueError:
                    tempFactor = None
                tempfactors.append(tempFactor)

                segids.append(line[66:76].strip())
                atomtypes.append(line[76:78].strip())

        # If segids not present, try to use chainids
        if not any(segids):
            segids, chainids = chainids, None
        if not icodes:
            icodes = None

        n_atoms = len(serials)
        attrs = []
        attrs.append(Atomnames(np.array(names, dtype=object)))
        # Optional attributes
        for vals, Attr, dtype in (
                (icodes, ICodes, object),
                (atomtypes, Atomtypes, object),
                (altlocs, AltLocs, object),
                (chainids, ChainIDs, object),
                (serials, Atomids, np.int32),
                (tempfactors, Tempfactors, np.float32),
                (occupancies, Occupancies, np.float32),
        ):
            # Skip if:
            # - vals is None
            # - any entries are None
            # - all entries are empty
            if vals is None:
                continue
            if any(v == None for v in vals):
                continue
            attrs.append(Attr(np.array(vals, dtype=dtype)))

        # Residue level stuff from here
        resnames = np.array(resnames, dtype=object)
        resids = np.array(resids, dtype=np.int32)
        segids = np.array(segids, dtype=object)

        if resnums:
            resnums = np.array(resnums, dtype=np.int32)
            residx, resids, (resnames, resnums, segids) = squash_by(
                resids, resnames, resnums, segids)
            attrs.append(Resnums(resnums))
        else:
            residx, resids, (resnames, segids) = squash_by(
                resids, resnames, segids)
        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnames(resnames))

        if any(segids) and not any(val == None for val in segids):
            segidx, segids = squash_by(segids)[:2]
            n_segments = len(segids)
            attrs.append(Segids(segids))
        else:
            n_segments = 1
            segidx = None

        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top

    def _parsebonds(self, serials):
        # Could optimise this by saving lines in the main loop
        # then doing post processing after all Atoms have been read
        # ie do one pass through the file only
        # Problem is that in multiframe PDB, the CONECT is at end of file,
        # so the "break" call happens before bonds are reached.

        # Mapping between the atom array indicies a.index and atom ids
        # (serial) in the original PDB file
        mapping = dict((s, i) for i, s in enumerate(serials))

        bonds = set()
        with openany(self.filename, "r") as f:
            lines = (line for line in f if line[:6] == "CONECT")
            for line in lines:
                atom, atoms = _parse_conect(line.strip())
                for a in atoms:
                    bond = tuple([mapping[atom], mapping[a]])
                    bonds.add(bond)

        bonds = tuple(bonds)

        return Bonds(bonds)


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

