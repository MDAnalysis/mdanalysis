# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
PDBQT topology parser
=====================

Use a PDBQT_ file to build a minimum internal structure representation (list of
atoms), including AutoDock_ atom types (stored as :attr:`Atom.type`) and
partial charges (:attr:`Atom.charge`).

* Reads a PDBQT file line by line and does not require sequential atom numbering.
* Multi-model PDBQT files are not supported.

.. Note:: Only reads atoms and their names; connectivity is not
          deduced. Masses are guessed and set to 0 if unknown.

.. SeeAlso:: `MDAnalysis.coordinates.PDBQT`

.. _PDBQT:
   http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
.. _AutoDock:
   http://autodock.scripps.edu/

Classes
-------

.. autoclass:: PDBQTParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

import numpy as np

from . import guessers
from ..lib import util
from .base import TopologyReader, change_squash
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    AltLocs,
    Atomtypes,
    Charges,
    Elements,
    Masses,
    Occupancies,
    Resids,
    Resnums,
    Resnames,
    Segids,
    Tempfactors,
)


class PDBQTParser(TopologyReader):
    """Read topology from a PDBQT file.

    Creates the following Attributes:
     - atom ids (serial)
     - atom types
     - atom names
     - altLocs
     - resnames
     - chainIDs (becomes segid)
     - resids
     - icodes
     - occupancies
     - tempfactors
     - charges

    Guesses the following:
     - elements
     - masses
    """
    format = 'PDBQT'

    def parse(self):
        """Parse atom information from PDBQT file *filename*.

        Returns
        -------
        MDAnalysis Topology object
        """
        serials = []
        names = []
        altlocs = []
        resnames = []
        chainids = []
        resids = []
        icodes = []
        occupancies = []
        tempfactors = []
        charges = []
        atomtypes = []

        with util.openany(self.filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith(('ATOM', 'HETATM')):
                    continue
                serials.append(int(line[6:11]))
                names.append(line[12:16].strip())
                altlocs.append(line[16:17].strip())
                resnames.append(line[17:21].strip())
                chainids.append(line[21:22].strip())
                resids.append(int(line[22:26]))
                icodes.append(line[26:27].strip())
                occupancies.append(float(line[54:60]))
                tempfactors.append(float(line[60:66]))
                charges.append(float(line[66:76]))
                atomtypes.append(line[77:80].strip())

        n_atoms = len(serials)

        elements = guessers.guess_types(names)
        masses = guessers.guess_masses(elements)

        attrs = []
        for attrlist, Attr, dtype in (
                (serials, Atomids, np.int32),
                (names, Atomnames, object),
                (altlocs, AltLocs, object),
                (occupancies, Occupancies, np.float32),
                (tempfactors, Tempfactors, np.float32),
                (charges, Charges, np.float32),
                (atomtypes, Atomtypes, object),
        ):
            attrs.append(Attr(np.array(attrlist, dtype=dtype)))
        attrs.append(Elements(elements, guessed=True))
        attrs.append(Masses(masses, guessed=True))

        resids = np.array(resids, dtype=np.int32)
        icodes = np.array(icodes, dtype=object)
        resnames = np.array(resnames, dtype=object)
        chainids = np.array(chainids, dtype=object)

        residx, (resids, icodes, resnames, chainids) = change_squash(
            (resids, icodes), (resids, icodes, resnames, chainids))
        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(Resnames(resnames))

        segidx, (segids,) = change_squash((chainids,), (chainids,))
        n_segments = len(segids)
        attrs.append(Segids(segids))

        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top
