# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
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
import numpy as np
import sqlite3
import os

from .base import TopologyReaderBase, change_squash
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Bonds,
    Charges,
    ChainIDs,
    Masses,
    Resids,
    Resnums,
    Resnames,
    Segids,
    AtomAttr,  # for custom Attributes
)


class Atomnums(AtomAttr):
    """The number for each Atom"""
    attrname = 'atomnums'
    singular = 'atomnum'


class DMSParser(TopologyReaderBase):
    """Read a topology from a DESRES_ Molecular Structure file.

    Format (DMS_) coordinate files (as used by the Desmond_ MD package).

    Reads the following attributes:
      Atom:
        - Atomids
        - Atomnums
        - Atomnames
        - Masses
        - Charges
        - Chainids
      Residue:
        - Resnames
        - Resids
      Segment:
        - Segids

    .. note::

        By default, atomtypes will be guessed on Universe creation.
        This may change in release 3.0.
        See :ref:`Guessers` for more information.

    .. _DESRES: http://www.deshawresearch.com
    .. _Desmond: http://www.deshawresearch.com/resources_desmond.html
    .. _DMS: http://www.deshawresearch.com/Desmond_Users_Guide-0.7.pdf
    .. versionchanged:: 2.8.0
        Removed type guessing (attributes guessing takes place now
        through universe.guess_TopologyAttrs() API).

    """
    format = 'DMS'

    def parse(self, **kwargs):
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
                    errmsg = "Failed reading the atoms from DMS Database"
                    raise IOError(errmsg) from None
                else:
                    attrs[attrname] = np.array(vals, dtype=dt)

            try:
                cur.row_factory = dict_factory
                cur.execute('SELECT * FROM bond')
                bonds = cur.fetchall()
            except sqlite3.DatabaseError:
                errmsg = "Failed reading the bonds from DMS Database"
                raise IOError(errmsg) from None
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
        atom_residx, (res_resids,
                      res_resnums,
                      res_resnames,
                      res_segids) = change_squash(
            (attrs['resid'], attrs['resname'], attrs['segid']),
            (attrs['resid'],
             attrs['resid'].copy(),
             attrs['resname'],
             attrs['segid']),
            )

        n_residues = len(res_resids)
        topattrs.append(Resids(res_resids))
        topattrs.append(Resnums(res_resnums))
        topattrs.append(Resnames(res_resnames))

        if any(res_segids) and not any(val is None for val in res_segids):
            res_segidx, (res_segids,) = change_squash((res_segids,),
                                                      (res_segids,))

            uniq_seg = np.unique(res_segids)
            idx2seg = {idx: res_segids[idx] for idx in res_segidx}
            res_segids = uniq_seg
            nidx = {segid: nidx for nidx, segid in enumerate(uniq_seg)}

            res_segidx = np.array([nidx[idx2seg[idx]] for idx in res_segidx])

            n_segments = len(res_segids)
            topattrs.append(Segids(res_segids))
        else:
            n_segments = 1
            topattrs.append(Segids(np.array(['SYSTEM'], dtype=object)))
            res_segidx = None

        topattrs.append(Bonds(attrs['bond']))

        top = Topology(len(attrs['id']), n_residues, n_segments,
                       attrs=topattrs,
                       atom_resindex=atom_residx,
                       residue_segindex=res_segidx)

        return top
