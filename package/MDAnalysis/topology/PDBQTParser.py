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
PDBQT topology parser
=====================

Use a PDBQT_ file to build a minimum internal structure representation (list of
atoms), including AutoDock_ atom types (stored as :attr:`Atom.type`) and
partial charges (:attr:`Atom.charge`).

* Reads a PDBQT file line by line and does not require sequential atom numbering.
* Multi-model PDBQT files are not supported.

Notes
-----
Only reads atoms and their names; connectivity is not
deduced.


See Also
--------
:class:`MDAnalysis.coordinates.PDBQT`


Classes
-------
.. autoclass:: PDBQTParser
   :members:
   :inherited-members:


.. _PDBQT:
   http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
.. _AutoDock:
   http://autodock.scripps.edu/
"""
import numpy as np

from ..lib import util
from .base import TopologyReaderBase, change_squash
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    AltLocs,
    Atomtypes,
    Charges,
    ICodes,
    Occupancies,
    RecordTypes,
    Resids,
    Resnums,
    Resnames,
    Segids,
    ChainIDs,
    Tempfactors,
)


class PDBQTParser(TopologyReaderBase):
    """Read topology from a PDBQT file.

    Creates the following Attributes:
     - atom ids (serial)
     - atom types
     - atom names
     - altLocs
     - resnames
     - chainIDs (assigned to segid as well)
     - resids
     - record_types (ATOM/HETATM)
     - icodes
     - occupancies
     - tempfactors
     - charges


    .. versionchanged:: 0.18.0
       Added parsing of Record types
    .. versionchanged:: 2.7.0
       Columns 67 - 70 in ATOM records, corresponding to the field *footnote*,
       are now ignored. See Autodock's `reference`_.
    .. versionchanged:: 2.8.0
        Removed mass guessing (attributes guessing takes place now
        through universe.guess_TopologyAttrs() API).

       .. _reference: 
          https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf
    """
    format = 'PDBQT'

    def parse(self, **kwargs):
        """Parse atom information from PDBQT file *filename*.

        Returns
        -------
        MDAnalysis Topology object
        """
        record_types = []
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

        with util.openany(self.filename) as f:
            for line in f:
                line = line.strip()
                if not line.startswith(('ATOM', 'HETATM')):
                    continue
                record_types.append(line[:6].strip())
                serials.append(int(line[6:11]))
                names.append(line[12:16].strip())
                altlocs.append(line[16:17].strip())
                resnames.append(line[17:21].strip())
                chainids.append(line[21:22].strip())
                resids.append(int(line[22:26]))
                icodes.append(line[26:27].strip())
                occupancies.append(float(line[54:60]))
                tempfactors.append(float(line[60:66]))
                charges.append(float(line[70:76]))
                atomtypes.append(line[77:80].strip())

        n_atoms = len(serials)

        attrs = []
        for attrlist, Attr, dtype in (
                (record_types, RecordTypes, object),
                (serials, Atomids, np.int32),
                (names, Atomnames, object),
                (altlocs, AltLocs, object),
                (occupancies, Occupancies, np.float32),
                (tempfactors, Tempfactors, np.float32),
                (charges, Charges, np.float32),
                (atomtypes, Atomtypes, object),
        ):
            attrs.append(Attr(np.array(attrlist, dtype=dtype)))

        resids = np.array(resids, dtype=np.int32)
        icodes = np.array(icodes, dtype=object)
        resnames = np.array(resnames, dtype=object)
        chainids = np.array(chainids, dtype=object)

        attrs.append(ChainIDs(chainids))

        residx, (resids, icodes, resnames, chainids) = change_squash(
            (resids, icodes), (resids, icodes, resnames, chainids))
        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(ICodes(icodes))
        attrs.append(Resnames(resnames))

        segidx, (segids,) = change_squash((chainids,), (chainids,))
        n_segments = len(segids)
        attrs.append(Segids(segids))

        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top
