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

from . import guessers
from ..lib.util import openany
from .base import TopologyReaderBase, squash_by
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Bonds,
    Charges,
    Masses,
    Resids,
    Resnums,
    Resnames,
    Segids,
)
from ..core.topology import Topology


class MOL2Parser(TopologyReaderBase):
    """Reads topology from a Tripos_ MOL2_ file.

    Creates the following Attributes:
     - Atomids
     - Atomnames
     - Atomtypes
     - Charges
     - Resids,
     - Resnames
     - Bonds

    Guesses the following:
     - masses

    .. versionchanged:: 0.9
       Now subclasses TopologyReaderBase
    .. versionchanged:: 0.20.0
       Allows for comments at the top of the file
       Ignores status bit strings
    .. versionchanged:: 2.0.0
       Bonds attribute is not added if no bonds are present in MOL2 file
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

        for a in atom_lines:
            aid, name, x, y, z, atom_type, resid, resname, charge = a.split()[:9]

            ids.append(aid)
            names.append(name)
            types.append(atom_type)
            resids.append(resid)
            resnames.append(resname)
            charges.append(charge)

        n_atoms = len(ids)

        masses = guessers.guess_masses(types)

        attrs = []
        attrs.append(Atomids(np.array(ids, dtype=np.int32)))
        attrs.append(Atomnames(np.array(names, dtype=object)))
        attrs.append(Atomtypes(np.array(types, dtype=object)))
        attrs.append(Charges(np.array(charges, dtype=np.float32)))
        attrs.append(Masses(masses, guessed=True))

        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)

        residx, resids, (resnames,) = squash_by(
            resids, resnames)
        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(Resnames(resnames))

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
