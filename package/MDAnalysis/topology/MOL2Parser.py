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
MOL2 file format --- :mod:`MDAnalysis.coordinates.MOL2`
========================================================

Classes to read Tripos_ molecule structure format (MOL2_) coordinate
and topology files. Used by the DOCK_ docking code.

.. _MOL2: http://www.tripos.com/data/support/mol2.pdf
.. _Tripos: http://www.tripos.com/
.. _DOCK: http://dock.compbio.ucsf.edu/


Classes
-------

.. autoclass:: MOL2Parser
   :members:
   :inherited-members:

"""
import os
import numpy as np

from ..lib.util import openany
from .base import TopologyReaderBase, squash_by
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Bonds,
    Charges,
    Elements,
    Resids,
    Resnums,
    Resnames,
    Segids,
)
from ..core.topology import Topology
from ..guesser.tables import SYBYL2SYMB

import warnings


class MOL2Parser(TopologyReaderBase):
    """Reads topology from a Tripos_ MOL2_ file.

    Creates the following Attributes:
     - Atomids
     - Atomnames
     - Atomtypes
     - Charges
     - Resids
     - Resnames
     - Bonds
     - Elements

    
    .. note::

        By default, masses will be guessed on Universe creation.
        This may change in release 3.0.
        See :ref:`Guessers`_ for more information.


    Notes
    -----
    Elements are obtained directly from the SYBYL atom types. If some atoms have
    unknown atom types, they will be assigned an empty element record. If all
    atoms have unknown atom types, the elements attribute will not be set.

    Dealing with optional fields:

      1. ``resid`` will set to 1 when not provided.

      2. If no atoms have ``resname`` field, resnames attribute will not be set;
         If some atoms have ``resname`` while some do not,
         :exc:`ValueError` will occur.

      3. If "NO_CHARGES" shows up in "@<TRIPOS>MOLECULE" section
         and no atoms have the ``charge`` field, charges attribute will not be set;
         If "NO_CHARGES" shows up while ``charge`` field appears,
         :exc:`ValueError` will occur;
         If a charge model is specified while some atoms don't have ``charge``,
         :exc:`ValueError` will occur as well.

    Raises
    ------
    ValueError
      If some atoms have the optional field ``resname`` (aka ``subst_name``)
      while some do not.

    ValueError
      If "NO_CHARGES" shows up in "@<TRIPOS>MOLECULE" section while
      some atoms have the optional field ``charge``.

    ValueError
      If a charge model is specified in "@<TRIPOS>MOLECULE" section while
      some atoms do not have the optional field ``charge``.


    .. versionchanged:: 0.9
       Now subclasses TopologyReaderBase
    .. versionchanged:: 0.20.0
       Allows for comments at the top of the file
       Ignores status bit strings
    .. versionchanged:: 2.0.0
       Bonds attribute is not added if no bonds are present in MOL2 file
    .. versionchanged:: 2.0.0
       Parse elements from atom types.
    .. versionchanged:: 2.2.0
       Read MOL2 files with optional columns omitted.
    .. versionchanged:: 2.8.0
        Removed mass guessing (attributes guessing takes place now
        through universe.guess_TopologyAttrs() API).

    """
    format = 'MOL2'

    def parse(self, **kwargs):
        """Parse MOL2 file *filename* and return the dict `structure`.

        Returns
        -------
        A MDAnalysis Topology object
        """
        blocks = []

        with openany(self.filename) as f:
            for i, line in enumerate(f):
                # found new molecules
                if "@<TRIPOS>MOLECULE" in line:
                    if len(blocks):
                        break
                    blocks.append({"start_line": i, "lines": []})
                if len(blocks):
                    blocks[-1]["lines"].append(line)

        if not len(blocks):
            raise ValueError("The mol2 file '{0}' needs to have at least one"
                             " @<TRIPOS>MOLECULE block".format(self.filename))
        block = blocks[0]

        sections = {}
        cursor = None

        for line in block["lines"]:
            if "@<TRIPOS>" in line:
                cursor = line.split("@<TRIPOS>")[1].strip().lower()
                sections[cursor] = []
                continue
            elif line.startswith("#") or line == "\n":
                continue
            sections[cursor].append(line)

        atom_lines, bond_lines = sections["atom"], sections.get("bond")

        if not len(atom_lines):
            raise ValueError("The mol2 block ({0}:{1}) has no atoms".format(
                os.path.basename(self.filename), block["start_line"]))

        ids = []
        names = []
        types = []
        resids = []
        resnames = []
        charges = []
        has_charges = sections['molecule'][3].strip() != 'NO_CHARGES'
        for a in atom_lines:
            columns = a.split()
            if len(columns) >= 9:
                aid, name, x, y, z, atom_type, \
                    resid, resname, charge = columns[:9]
            elif len(columns) < 6:
                raise ValueError(f"The @<TRIPOS>ATOM block in mol2 file"
                                 f" {os.path.basename(self.filename)}"
                                 f" should have at least 6 fields to be"
                                 f" unpacked: atom_id atom_name x y z"
                                 f" atom_type [subst_id[subst_name"
                                 f" [charge [status_bit]]]]")
            else:
                aid, name, x, y, z, atom_type = columns[:6]
                id_name_charge = [1, None, None]
                for i in range(6, len(columns)):
                    id_name_charge[i-6] = columns[i]
                resid, resname, charge = id_name_charge
            if has_charges:
                if charge is None:
                    raise ValueError(f"The mol2 file {self.filename}"
                                     f" indicates a charge model"
                                     f"{sections['molecule'][3]}, but"
                                     f" no charge provided in line: {a}")
            else:
                if charge is not None:
                    raise ValueError(f"The mol2 file {self.filename}"
                                     f" indicates no charges, but charge"
                                     f" {charge} provided in line: {a}.")

            ids.append(aid)
            names.append(name)
            types.append(atom_type)
            resids.append(resid)
            resnames.append(resname)
            charges.append(charge)

        n_atoms = len(ids)

        validated_elements = np.empty(n_atoms, dtype="U3")
        invalid_elements = set()
        for i, at in enumerate(types):
            if at in SYBYL2SYMB:
                validated_elements[i] = SYBYL2SYMB[at]
            else:
                invalid_elements.add(at)
                validated_elements[i] = ''

        # Print single warning for all unknown elements, if any
        if invalid_elements:
            warnings.warn("Unknown elements found for some "
                          f"atoms: {invalid_elements}. "
                          "These have been given an empty element record.")

        attrs = []
        attrs.append(Atomids(np.array(ids, dtype=np.int32)))
        attrs.append(Atomnames(np.array(names, dtype=object)))
        attrs.append(Atomtypes(np.array(types, dtype=object)))
        if has_charges:
            attrs.append(Charges(np.array(charges, dtype=np.float32)))

        if not np.all(validated_elements == ''):
            attrs.append(Elements(validated_elements))

        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)

        if np.all(resnames):
            residx, resids, (resnames,) = squash_by(
                resids, resnames)
            n_residues = len(resids)
            attrs.append(Resids(resids))
            attrs.append(Resnums(resids.copy()))
            attrs.append(Resnames(resnames))
        elif not np.any(resnames):
            residx, resids, _ = squash_by(resids,)
            n_residues = len(resids)
            attrs.append(Resids(resids))
            attrs.append(Resnums(resids.copy()))
        else:
            raise ValueError(f"Some atoms in the mol2 file {self.filename}"
                             f" have subst_name while some do not.")

        attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))

        # don't add Bonds if there are none (Issue #3057)
        if bond_lines:
            bonds = []
            bondorder = []
            for b in bond_lines:
                # bond_type can be: 1, 2, am, ar
                bid, a0, a1, bond_type = b.split()[:4]

                a0, a1 = int(a0) - 1, int(a1) - 1
                bond = tuple(sorted([a0, a1]))
                bondorder.append(bond_type)
                bonds.append(bond)
            attrs.append(Bonds(bonds, order=bondorder))

        top = Topology(n_atoms, n_residues, 1,
                       attrs=attrs,
                       atom_resindex=residx)

        return top
