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
RDKit topology parser
=====================

Converts a `RDKit <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol>`_ 
:class:`rdkit.Chem.rdchem.Mol` into a :class:`MDAnalysis.core.Topology`.


Example
-------
TODO

Classes
-------

.. autoclass:: RDKitParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import, division

import logging
import numpy as np

from .base import TopologyReaderBase, change_squash
from . import guessers
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Elements,
    Masses,
    Charges,
    Bonds,
    Resids,
    Resnums,
    Resnames,
    Segids,
    AltLocs,
    ChainIDs,
    Occupancies,
    Tempfactors,
)
from ..core.topology import Topology

logger = logging.getLogger("MDAnalysis.topology.RDKitParser")


class RDKitParser(TopologyReaderBase):
    """
    For RDKit structures
    """
    format = 'RDKIT'

    @staticmethod
    def _format_hint(thing):
        """Can this Parser read object *thing*?

        .. versionadded:: 1.0.0
        """
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
        MDAnalysis *Topology* object
        """

        mol = self.filename

        # Atoms
        names = []
        resnums = []
        resnames = []
        elements = []
        masses = []
        charges = []
        ids = []
        atomtypes = []
        segids = []
        altlocs = []
        chainids = []
        occupancies = []
        tempfactors = []
        for atom in mol.GetAtoms():
            ids.append(atom.GetIdx())
            elements.append(atom.GetSymbol())
            masses.append(atom.GetMass())
            charges.append(atom.GetFormalCharge())
            mi = atom.GetMonomerInfo()
            if mi: # atom name and residue info are present
                names.append(mi.GetName().strip())
                resnums.append(mi.GetResidueNumber())
                resnames.append(mi.GetResidueName())
                segids.append(mi.GetSegmentNumber())
                altlocs.append(mi.GetAltLoc().strip())
                chainids.append(mi.GetChainId().strip())
                occupancies.append(mi.GetOccupancy())
                tempfactors.append(mi.GetTempFactor())
            else:
                for prop, value in atom.GetPropsAsDict(True).items():
                    if 'atomname' in prop.lower(): # usually _TriposAtomName
                        names.append(value)
                    elif 'atomtype' in prop.lower(): # usually _TriposAtomType
                        atomtypes.append(value)
                
        # make Topology attributes
        attrs = []
        n_atoms = len(ids)

        # * Attributes always present *

        # Atom attributes
        for vals, Attr, dtype in (
            (ids, Atomids, np.int32),
            (elements, Elements, object),
            (masses, Masses, np.float32),
            (charges, Charges, np.float32),
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
            residx, (resnums, resnames, segids) = change_squash(
                (resnums, resnames, segids), 
                (resnums, resnames, segids))
            n_residues = len(resnums)
            for vals, Attr, dtype in (
                (resnums, Resids, np.int32),
                (resnums.copy(), Resnums, np.int32),
                (resnames, Resnames, object)
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