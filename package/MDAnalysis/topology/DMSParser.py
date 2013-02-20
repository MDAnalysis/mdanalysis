# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
DESRES Molecular Structure file format topology parser
======================================================

Classes to read a topology from a DESRES_ Molecular Structure file
format (DMS_) coordinate files (as used by the Desmond_ MD package).

.. _DESRES: http://www.deshawresearch.com
.. _Desmond: http://www.deshawresearch.com/resources_desmond.html
.. _DMS: http://www.deshawresearch.com/Desmond_Users_Guide-0.7.pdf
"""

from __future__ import with_statement

from MDAnalysis.core.AtomGroup import Atom
from MDAnalysis.topology.core import guess_atom_type, Bond
import numpy, sqlite3

class DMSParseError(Exception):
        pass

def parse(filename):
        """Parse DMS file *filename* and return the dict `structure`.

        Only reads the list of atoms.

        :Returns: MDAnalysis internal *structure* dict, which contains
                  Atom and Bond objects

        .. SeeAlso:: The *structure* dict is defined in
                     :func:`MDAnalysis.topology.PSFParser.parse`.

        """
        con = sqlite3.connect(filename)

        def dict_factory(cursor, row):
            """
            Fetch SQL records as dictionaries, rather than the default tuples.
            """
            d = {}
            for idx, col in enumerate(cursor.description):
                d[col[0]] = row[idx]
            return d

        atoms = None
        with con:
            # This will return dictionaries instead of tuples, when calling cur.fetch() or fetchall()
            con.row_factory = dict_factory
            cur = con.cursor()
            cur.execute('SELECT * FROM particle')
            particles = cur.fetchall()

            # p["anum"] contains the atomic number
            atoms = [ (p["id"], Atom(p["id"], p["name"].strip(), guess_atom_type(p["name"].strip()), p["resname"].strip(), p["resid"], p["segid"].strip(), p["mass"], p["charge"]))  for p in particles]

            atoms_dictionary = dict(atoms)

            cur.execute('SELECT * FROM bond')
            bonds = cur.fetchall()

            bonds = [ Bond(atoms_dictionary[b["p0"]],
                           atoms_dictionary[b["p1"]],
                           b["order"] )  for b in bonds]

        # All the records below besides donors and acceptors can be contained in a DMS file.
        # In addition to the coordinates and bonds, DMS may contain the entire force-field information (terms+parameters),
        return {"_atoms": [ atom[1] for atom in atoms],
                "_bonds": bonds,
                "_angles": [],
                "_dihe": [],
                "_impr": [],
                "_donors": [],
                "_acceptors": [],
                }
