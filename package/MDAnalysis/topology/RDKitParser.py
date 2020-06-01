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
from six.moves import range
from six import raise_from

import logging
import numpy as np

from .base import TopologyReaderBase, squash_by
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Elements,
    Masses,
    Bonds,
    Resids,
    Resnums,
    Resnames,
    Segids,
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
        residue_numbers = []
        residue_names = []
        elements = []
        masses = []
        ids = []
        for atom in mol.GetAtoms():
            ids.append(atom.GetIdx())
            elements.append(atom.GetSymbol())
            masses.append(atom.GetMass())
            mi = atom.GetMonomerInfo()
            if mi: # atom name and residue info are present
                names.append(mi.GetName().strip())
                residue_numbers.append(mi.GetResidueNumber())
                residue_names.append(mi.GetResidueName())
            else:
                pass

        # make Topology attributes
        residx, residue_numbers, (residue_names,) = squash_by(
            residue_numbers, np.array(residue_names))
        attrs = []
        n_atoms = len(names)
        n_residues = len(residue_numbers)

        for vals, Attr, dtype in (
            (ids, Atomids, np.int32),
            (elements, Elements, object),
            (names, Atomnames, object),
            (masses, Masses, np.float32),
            (residue_numbers, Resids, np.int32),
            (residue_numbers.copy(), Resnums, np.int32),
            (residue_names, Resnames, object)
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

        attrs.append(Bonds(bonds, types=bond_types, guessed=False, 
                           order=bond_orders))

        # Others
        #TODO
        attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))

        top = Topology(n_atoms, n_residues, 1,
                       attrs=attrs,
                       atom_resindex=residx)

        return top