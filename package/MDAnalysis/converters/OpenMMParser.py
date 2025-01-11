# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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

"""OpenMM topology parser :mod:`MDAnalysis.converters.OpenMMParser`
===================================================================

.. versionadded:: 2.0.0
.. versionchanged:: 2.8.0
   Removed type and mass guessing (attributes guessing takes place
   now through universe.guess_TopologyAttrs() API)


Converts an
`OpenMM topology <http://docs.openmm.org/latest/api-python/generated/openmm.app.topology.Topology.html#openmm.app.topology.Topology>`_
:class:`openmm.app.topology.Topology` into a :class:`MDAnalysis.core.Topology`.

Also converts some objects within the
`OpenMM Application layer <http://docs.openmm.org/latest/api-python/app.html>`_

- `openmm.app.pdbfile.PDBFile <http://docs.openmm.org/latest/api-python/generated/openmm.app.pdbfile.PDBFile.html#openmm.app.pdbfile.PDBFile>`_
- `openmm.app.simulation.Simulation <http://docs.openmm.org/latest/api-python/generated/openmm.app.simulation.Simulation.html#openmm.app.simulation.Simulation>`_
- `openmm.app.modeller.Modeller <http://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller>`_
- `openmm.app.pdbxfile.PDBxFile <http://docs.openmm.org/latest/api-python/generated/openmm.app.pdbxfile.PDBxFile.html#openmm.app.pdbxfile.PDBxFile>`_

The :class:`OpenMMTopologyParser` generates a topology from an OpenMM Topology object.


Classes
-------

.. autoclass:: OpenMMTopologyParser
   :members:
   :inherited-members:

.. autoclass:: OpenMMAppTopologyParser
   :members:
   :inherited-members:

"""

import numpy as np
import warnings

from ..topology.base import TopologyReaderBase
from ..guesser.tables import SYMB2Z
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Bonds,
    ChainIDs,
    Elements,
    Masses,
    Resids,
    Resnums,
    Resnames,
    Segids,
)


class OpenMMTopologyParser(TopologyReaderBase):
    format = "OPENMMTOPOLOGY"

    @staticmethod
    def _format_hint(thing):
        """Can this Parser read object *thing*?"""
        try:
            from openmm import app
        except ImportError:
            try:  # pragma: no cover
                from simtk.openmm import app
            except ImportError:
                return False
        else:
            return isinstance(thing, app.Topology)

    def _mda_topology_from_omm_topology(self, omm_topology):
        """Construct mda topology from omm topology

        Can be used for any openmm object that contains a topology object

        Parameters
        ----------
        omm_topology: openmm.Topology

        Returns
        -------
        top : MDAnalysis.core.topology.Topology


        When partial elements are present, values from available elements
        are used whereas the absent elements are assigned an empty string
        with their atomtype set to 'X' and mass set to 0.0.

        For elements with invalid and unreal symbols, the symbol is used
        as it is for atomtypes but an empty string is used for elements.

        .. versionchanged:: 2.2.0
           The parser now works when element attribute is missing in some or
           all the atoms.

        """

        try:
            from openmm.unit import daltons
        except ImportError:
            try:
                from simtk.unit import daltons
            except ImportError:
                msg = (
                    "OpenMM is required for the OpenMMParser but "
                    "it's not installed. Try installing it with \n"
                    "conda install -c conda-forge openmm"
                )
                raise ImportError(msg)

        atom_resindex = [a.residue.index for a in omm_topology.atoms()]
        residue_segindex = [r.chain.index for r in omm_topology.residues()]
        atomids = [a.id for a in omm_topology.atoms()]
        atomnames = [a.name for a in omm_topology.atoms()]
        chainids = [a.residue.chain.id for a in omm_topology.atoms()]
        resnames = [r.name for r in omm_topology.residues()]
        resids = [r.index + 1 for r in omm_topology.residues()]
        resnums = resids.copy()
        segids = [str(c.index) for c in omm_topology.chains()]
        bonds = [(b.atom1.index, b.atom2.index) for b in omm_topology.bonds()]
        bond_orders = [b.order for b in omm_topology.bonds()]
        bond_types = [b.type for b in omm_topology.bonds()]

        attrs = [
            Atomids(np.array(atomids, dtype=np.int32)),
            Atomnames(np.array(atomnames, dtype=object)),
            Bonds(bonds, types=bond_types, order=bond_orders, guessed=False),
            ChainIDs(np.array(chainids, dtype=object)),
            Resids(resids),
            Resnums(resnums),
            Resnames(resnames),
            Segids(segids),
        ]

        validated_elements = []
        masses = []
        atomtypes = []
        for a in omm_topology.atoms():
            elem = a.element
            if elem is not None:
                if elem.symbol.capitalize() in SYMB2Z:
                    validated_elements.append(elem.symbol)
                else:
                    validated_elements.append("")
                atomtypes.append(elem.symbol)
                masses.append(elem.mass.value_in_unit(daltons))
            else:
                validated_elements.append("")
                masses.append(0.0)
                atomtypes.append("X")

        if not all(validated_elements):
            if any(validated_elements):
                warnings.warn(
                    "Element information missing for some atoms. "
                    "These have been given an empty element record "
                )
                if any(i == "X" for i in atomtypes):
                    warnings.warn(
                        "For absent elements, atomtype has been  "
                        "set to 'X' and mass has been set to 0.0. "
                        "If needed these can be guessed using "
                        "universe.guess_TopologyAttrs("
                        "to_guess=['masses', 'types']). "
                        "(for MDAnalysis version 2.x "
                        "this is done automatically,"
                        " but it will be removed in 3.0)."
                    )

                attrs.append(
                    Elements(np.array(validated_elements, dtype=object))
                )

            else:
                wmsg = (
                    "Element information is missing for all the atoms. "
                    "Elements attribute will not be populated. "
                    "Atomtype attribute will be guessed using atom "
                    "name and mass will be guessed using atomtype."
                    "For MDAnalysis version 2.x this is done automatically, "
                    "but it will be removed in MDAnalysis v3.0. "
                    "These can be guessed using "
                    "universe.guess_TopologyAttrs("
                    "to_guess=['masses', 'types']) "
                    "See MDAnalysis.guessers."
                )

                warnings.warn(wmsg)
        else:
            attrs.append(Elements(np.array(validated_elements, dtype=object)))
        attrs.append(Atomtypes(np.array(atomtypes, dtype=object)))
        attrs.append(Masses(np.array(masses)))

        n_atoms = len(atomids)
        n_residues = len(resids)
        n_segments = len(segids)
        top = Topology(
            n_atoms,
            n_residues,
            n_segments,
            attrs=attrs,
            atom_resindex=atom_resindex,
            residue_segindex=residue_segindex,
        )

        return top

    def parse(self, **kwargs):
        omm_topology = self.filename
        top = self._mda_topology_from_omm_topology(omm_topology)

        return top


class OpenMMAppTopologyParser(OpenMMTopologyParser):
    format = "OPENMMAPP"

    @staticmethod
    def _format_hint(thing):
        """Can this Parser read object *thing*?"""
        try:
            from openmm import app
        except ImportError:
            try:  # pragma: no cover
                from simtk.openmm import app
            except ImportError:
                return False
        else:
            return isinstance(
                thing,
                (app.PDBFile, app.Modeller, app.Simulation, app.PDBxFile),
            )

    def parse(self, **kwargs):
        try:
            omm_topology = self.filename.getTopology()
        except AttributeError:
            omm_topology = self.filename.topology
        top = self._mda_topology_from_omm_topology(omm_topology)

        return top
