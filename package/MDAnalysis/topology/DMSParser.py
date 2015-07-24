# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
DESRES Molecular Structure file format topology parser
======================================================

Classes to read a topology from a DESRES_ Molecular Structure file
format (DMS_) coordinate files (as used by the Desmond_ MD package).

.. _DESRES: http://www.deshawresearch.com
.. _Desmond: http://www.deshawresearch.com/resources_desmond.html
.. _DMS: http://www.deshawresearch.com/Desmond_Users_Guide-0.7.pdf

Classes
-------

.. autoclass:: DMSParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

import sqlite3
import os

from ..core.AtomGroup import Atom
from .core import guess_atom_type
from .base import TopologyReader


class DMSParser(TopologyReader):
    """Read a topology from a DESRES_ Molecular Structure file.

    Format (DMS_) coordinate files (as used by the Desmond_ MD package).

    .. _DESRES: http://www.deshawresearch.com
    .. _Desmond: http://www.deshawresearch.com/resources_desmond.html
    .. _DMS: http://www.deshawresearch.com/Desmond_Users_Guide-0.7.pdf
    """

    def parse(self):
        """Parse DMS file *filename* and return the dict `structure`.

        Only reads the list of atoms.

        :Returns: MDAnalysis internal *structure* dict, which contains
                  Atom and Bond objects

        .. SeeAlso:: The *structure* dict is defined in
                     :mod:`MDAnalysis.topology`.

        """
        # Fix by SB: Needed because sqlite3.connect does not raise anything if file is not there
        if not os.path.isfile(self.filename):
            raise IOError("No such file: {0}".format(self.filename))

        def dict_factory(cursor, row):
            """
            Fetch SQL records as dictionaries, rather than the default tuples.
            """
            d = {}
            for idx, col in enumerate(cursor.description):
                d[col[0]] = row[idx]
            return d

        with sqlite3.connect(self.filename) as con:
            try:
                # This will return dictionaries instead of tuples,
                # when calling cur.fetch() or fetchall()
                con.row_factory = dict_factory
                cur = con.cursor()
                cur.execute('SELECT * FROM particle')
                particles = cur.fetchall()
            except sqlite3.DatabaseError:
                raise IOError("Failed reading the atoms from DMS Database")
            else:
                # p["anum"] contains the atomic number
                try:
                    atoms = [Atom(p["id"], p["name"].strip(),
                                  guess_atom_type(p["name"].strip()),
                                  p["resname"].strip(), p["resid"],
                                  p["segid"].strip(), p["mass"], p["charge"],
                                  universe=self._u)
                             for p in particles]
                except KeyError:
                    raise ValueError("Failed reading atom information")
            try:
                cur.execute('SELECT * FROM bond')
                bonds = cur.fetchall()
            except sqlite3.DatabaseError:
                raise IOError("Failed reading the bonds from DMS Database")
            else:
                bondlist = []
                bondorder = {}
                for b in bonds:
                    desc = tuple(sorted([b['p0'], b['p1']]))
                    bondlist.append(desc)
                    bondorder[desc] = b['order']

        # All the records below besides donors and acceptors can be contained in a DMS file.
        # In addition to the coordinates and bonds, DMS may contain the entire force-field
        # information (terms+parameters),
        structure = {"atoms": atoms,
                     "bonds": tuple(bondlist),
                     "bondorder": bondorder}

        return structure
