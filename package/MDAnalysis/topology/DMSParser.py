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

import numpy as np
import sqlite3
import os

from .base import TopologyReader, squash_by
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Bonds,
    Charges,
    Masses,
    Resids,
    Resnames,
    Segids,
    AtomAttr,  # for custom Attributes
    ChainIDs,
)


class Atomnums(AtomAttr):
    """The number for each Atom"""
    attrname = 'atomnums'
    singular = 'atomnum'


class DMSParser(TopologyReader):
    """Read a topology from a DESRES_ Molecular Structure file.

    Format (DMS_) coordinate files (as used by the Desmond_ MD package).

    Reads the following attributes:
      Atom:
        - atomid
        - atomnum
        - atomname
        - mass
        - charge
        - chainids
      Residue:
        - resname
        - resid
      Segment:
        - segid

    .. _DESRES: http://www.deshawresearch.com
    .. _Desmond: http://www.deshawresearch.com/resources_desmond.html
    .. _DMS: http://www.deshawresearch.com/Desmond_Users_Guide-0.7.pdf
    """

    def parse(self):
        """Parse DMS file *filename* and return the Topology object"""
        # Fix by SB: Needed because sqlite3.connect does not raise anything
        # if file is not there
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

        attrs = {}

        # Row factories for different data types
        facs = {np.int32: lambda c, r: r[0],
                np.float32: lambda c, r: r[0],
                object: lambda c, r: str(r[0].strip())}

        with sqlite3.connect(self.filename) as con:
            # Selecting single column, so just strip tuple
            for attrname, dt in [
                    ('id', np.int32),
                    ('anum', np.int32),
                    ('mass', np.float32),
                    ('charge', np.float32),
                    ('name', object),
                    ('resname', object),
                    ('resid', np.int32),
                    ('chain', object),
                    ('segid', object),
            ]:
                try:
                    cur = con.cursor()
                    cur.row_factory = facs[dt]
                    cur.execute('SELECT {} FROM particle'
                                ''.format(attrname))
                    vals = cur.fetchall()
                except sqlite3.DatabaseError:
                    raise IOError(
                        "Failed reading the atoms from DMS Database")
                else:
                    attrs[attrname] = np.array(vals, dtype=dt)

            try:
                cur.row_factory = dict_factory
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
                attrs['bond'] = bondlist
                attrs['bondorder'] = bondorder

        topattrs = []
        # Bundle in Atom level objects
        for attr, cls in [
                ('id', Atomids),
                ('anum', Atomnums),
                ('mass', Masses),
                ('charge', Charges),
                ('name', Atomnames),
                ('chain', ChainIDs),
        ]:
            topattrs.append(cls(attrs[attr]))

        # Residues
        atom_residx, res_resids, (res_resnames, res_segids) = squash_by(
            attrs['resid'], attrs['resname'], attrs['segid'])
        topattrs.append(Resids(res_resids))
        topattrs.append(Resnames(res_resnames))

        # Segments
        res_segidx, seg_segids = squash_by(
            res_segids)[:2]
        topattrs.append(Segids(seg_segids))

        # Bonds
        topattrs.append(Bonds(attrs['bond']))

        top = Topology(len(attrs['id']), len(res_resids), len(seg_segids),
                       attrs=topattrs,
                       atom_resindex=atom_residx,
                       residue_segindex=res_segidx)

        return top
