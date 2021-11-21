# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
RDKit topology parser --- :mod:`MDAnalysis.converters.RDKitParser`
==================================================================

Converts an `RDKit <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol>`_ :class:`rdkit.Chem.rdchem.Mol` into a :class:`MDAnalysis.core.Topology`.


See Also
--------
:mod:`MDAnalysis.converters.RDKit`


Classes
-------

.. autoclass:: RDKitParser
   :members:
   :inherited-members:

"""

import logging
import warnings
import numpy as np

from ..topology.base import TopologyReaderBase, change_squash
from ..topology import guessers
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Elements,
    Masses,
    Charges,
    Aromaticities,
    Bonds,
    Resids,
    Resnums,
    Resnames,
    RSChirality,
    Segids,
    AltLocs,
    ChainIDs,
    ICodes,
    Occupancies,
    Tempfactors,
)
from ..core.topology import Topology

logger = logging.getLogger("MDAnalysis.converters.RDKitParser")


def _rdkit_atom_to_RS(atom):
    """Fetches RDKit chiral tags"""
    try:
        return atom.GetProp("_CIPCode")
    except KeyError:
        return ""


class RDKitParser(TopologyReaderBase):
    """
    For RDKit structures

    Creates the following Attributes:
     - Atomids
     - Atomnames
     - Aromaticities
     - Elements
     - Masses
     - Bonds
     - Resids
     - Resnums
     - RSChirality
     - Segids

    Guesses the following:
     - Atomtypes

    Depending on RDKit's input, the following Attributes might be present:
     - Charges
     - Resnames
     - AltLocs
     - ChainIDs
     - ICodes
     - Occupancies
     - Tempfactors

    Attributes table:

    +---------------------------------------------+-------------------------+
    | RDKit attribute                             | MDAnalysis equivalent   |
    +=============================================+=========================+
    | atom.GetMonomerInfo().GetAltLoc()           | altLocs                 |
    +---------------------------------------------+-------------------------+
    | atom.GetIsAromatic()                        | aromaticities           |
    +---------------------------------------------+-------------------------+
    | atom.GetMonomerInfo().GetChainId()          | chainIDs                |
    +---------------------------------------------+-------------------------+
    | atom.GetDoubleProp('_GasteigerCharge')      | charges                 |
    | atom.GetDoubleProp('_TriposPartialCharge')  |                         |
    +---------------------------------------------+-------------------------+
    | atom.GetSymbol()                            | elements                |
    +---------------------------------------------+-------------------------+
    | atom.GetMonomerInfo().GetInsertionCode()    | icodes                  |
    +---------------------------------------------+-------------------------+
    | atom.GetIdx()                               | indices                 |
    +---------------------------------------------+-------------------------+
    | atom.GetMass()                              | masses                  |
    +---------------------------------------------+-------------------------+
    | atom.GetMonomerInfo().GetName()             | names                   |
    | atom.GetProp('_TriposAtomName')             |                         |
    +---------------------------------------------+-------------------------+
    | atom.GetProp('_CIPCode')                    | chiralities             |
    +---------------------------------------------+-------------------------+
    | atom.GetMonomerInfo().GetOccupancy()        | occupancies             |
    +---------------------------------------------+-------------------------+
    | atom.GetMonomerInfo().GetResidueName()      | resnames                |
    +---------------------------------------------+-------------------------+
    | atom.GetMonomerInfo().GetResidueNumber()    | resnums                 |
    +---------------------------------------------+-------------------------+
    | atom.GetMonomerInfo().GetTempFactor()       | tempfactors             |
    +---------------------------------------------+-------------------------+
    | atom.GetProp('_TriposAtomType')             | types                   |
    +---------------------------------------------+-------------------------+

    Raises
    ------
    ValueError
        If only part of the atoms have MonomerInfo available


    .. versionadded:: 2.0.0
    .. versionchanged:: 2.1.0
       Added R/S chirality support
    """
    format = 'RDKIT'

    @staticmethod
    def _format_hint(thing):
        """Can this Parser read object *thing*?"""
        try:
            from rdkit import Chem
        except ImportError:  # if no rdkit, probably not rdkit
            return False
        else:
            return isinstance(thing, Chem.Mol)

    def parse(self, **kwargs):
        """Parse RDKit into Topology

        Returns
        -------
        MDAnalysis Topology object
        """

        mol = self.filename

        # Atoms
        names = []
        chiralities = []
        resnums = []
        resnames = []
        elements = []
        masses = []
        charges = []
        aromatics = []
        ids = []
        atomtypes = []
        segids = []
        altlocs = []
        chainids = []
        icodes = []
        occupancies = []
        tempfactors = []

        try:
            atom = mol.GetAtomWithIdx(0)
        except RuntimeError:
            top = Topology(n_atoms=0, n_res=0, n_seg=0,
                           attrs=None,
                           atom_resindex=None,
                           residue_segindex=None)
            return top

        # check if multiple charges present
        if atom.HasProp('_GasteigerCharge') and (
        atom.HasProp('_TriposPartialCharge')
        ):
            warnings.warn(
                'Both _GasteigerCharge and _TriposPartialCharge properties '
                'are present. Using Gasteiger charges by default.')

        for atom in mol.GetAtoms():
            ids.append(atom.GetIdx())
            elements.append(atom.GetSymbol())
            masses.append(atom.GetMass())
            aromatics.append(atom.GetIsAromatic())
            chiralities.append(_rdkit_atom_to_RS(atom))
            mi = atom.GetMonomerInfo()
            if mi: # atom name and residue info are present
                names.append(mi.GetName().strip())
                resnums.append(mi.GetResidueNumber())
                resnames.append(mi.GetResidueName())
                segids.append(mi.GetSegmentNumber())
                altlocs.append(mi.GetAltLoc().strip())
                chainids.append(mi.GetChainId().strip())
                icodes.append(mi.GetInsertionCode().strip())
                occupancies.append(mi.GetOccupancy())
                tempfactors.append(mi.GetTempFactor())
            else:
                # atom name (MOL2 only)
                try:
                    names.append(atom.GetProp('_TriposAtomName'))
                except KeyError:
                    pass
                # atom type (MOL2 only)
                try:
                    atomtypes.append(atom.GetProp('_TriposAtomType'))
                except KeyError:
                    pass
                # gasteiger charge (computed):
                # if the user took the time to compute them, make it a priority
                # over charges read from a MOL2 file
                try:
                    charges.append(atom.GetDoubleProp('_GasteigerCharge'))
                except KeyError:
                    # partial charge (MOL2 only)
                    try:
                        charges.append(atom.GetDoubleProp('_TriposPartialCharge'))
                    except KeyError:
                        pass

        # make Topology attributes
        attrs = []
        n_atoms = len(ids)

        if resnums and (len(resnums) != n_atoms):
            raise ValueError(
                "ResidueInfo is only partially available in the molecule. "
                "If you have added hydrogens to the input RDKit molecule with "
                "`Chem.AddHs(mol)`, consider using "
                "`Chem.AddHs(mol, addResidueInfo=True)` instead"
            )

        # * Attributes always present *

        # Atom attributes
        for vals, Attr, dtype in (
            (ids, Atomids, np.int32),
            (elements, Elements, object),
            (masses, Masses, np.float32),
            (aromatics, Aromaticities, bool),
            (chiralities, RSChirality, 'U1'),
        ):
            attrs.append(Attr(np.array(vals, dtype=dtype)))

        # Bonds
        bonds = []
        bond_types = []
        bond_orders = []
        for bond in mol.GetBonds():
            bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            bond_orders.append(bond.GetBondTypeAsDouble())
            bond_types.append(str(bond.GetBondType()))
        attrs.append(Bonds(bonds, types=bond_types, order=bond_orders))

        # * Optional attributes *

        # Atom name
        if names:
            attrs.append(Atomnames(np.array(names, dtype=object)))
        else:
            for atom in mol.GetAtoms():
                name = "%s%d" % (atom.GetSymbol(), atom.GetIdx())
                names.append(name)
            attrs.append(Atomnames(np.array(names, dtype=object)))

        # Atom type
        if atomtypes:
            attrs.append(Atomtypes(np.array(atomtypes, dtype=object)))
        else:
            atomtypes = guessers.guess_types(names)
            attrs.append(Atomtypes(atomtypes, guessed=True))

        # Partial charges
        if charges:
            attrs.append(Charges(np.array(charges, dtype=np.float32)))
        else:
            pass # no guesser yet

        # PDB only
        for vals, Attr, dtype in (
            (altlocs, AltLocs, object),
            (chainids, ChainIDs, object),
            (occupancies, Occupancies, np.float32),
            (tempfactors, Tempfactors, np.float32),
        ):
            if vals:
                attrs.append(Attr(np.array(vals, dtype=dtype)))

        # Residue
        if any(resnums) and not any(val is None for val in resnums):
            resnums = np.array(resnums, dtype=np.int32)
            resnames = np.array(resnames, dtype=object)
            segids = np.array(segids, dtype=object)
            icodes = np.array(icodes, dtype=object)
            residx, (resnums, resnames, icodes, segids) = change_squash(
                (resnums, resnames, icodes, segids),
                (resnums, resnames, icodes, segids))
            n_residues = len(resnums)
            for vals, Attr, dtype in (
                (resnums, Resids, np.int32),
                (resnums.copy(), Resnums, np.int32),
                (resnames, Resnames, object),
                (icodes, ICodes, object),
            ):
                attrs.append(Attr(np.array(vals, dtype=dtype)))
        else:
            attrs.append(Resids(np.array([1])))
            attrs.append(Resnums(np.array([1])))
            residx = None
            n_residues = 1

        # Segment
        if any(segids) and not any(val is None for val in segids):
            segidx, (segids,) = change_squash((segids,), (segids,))
            n_segments = len(segids)
            attrs.append(Segids(segids))
        else:
            n_segments = 1
            attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))
            segidx = None

        # create topology
        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top
