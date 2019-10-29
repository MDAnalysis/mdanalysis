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
ITP topology parser
===================

Reads a GROMACS ITP_ file to build the system. The topology will
contain atom IDs, segids, residue IDs, residue names, atom names, atom types,
charges, chargegroups, masses (guessed if not found), and moltypes. 
Bonds, angles, dihedrals and impropers are also read from
the file.

.. _ITP: http://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#molecule-itp-file

Classes
-------

.. autoclass:: ITPParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import, division
from six.moves import range

import logging
import functools
import re
from math import ceil
import numpy as np

from ..lib.util import openany
from . import guessers
from .base import TopologyReaderBase, squash_by, change_squash
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Masses,
    Moltypes,
    Charges,
    Resids,
    Resnums,
    Resnames,
    Segids,
    Bonds,
    Angles,
    Dihedrals,
    Impropers,
    AtomAttr,
)
from ..core.topology import Topology

logger = logging.getLogger("MDAnalysis.topology.ITP")


class Chargegroups(AtomAttr):
    """The charge group for each Atom"""
    attrname = 'chargegroups'
    singular = 'chargegroup'

class ITPParser(TopologyReaderBase):
    """Read topology information from a GROMACS ITP_ file.

    Creates a Topology with the following Attributes:
    - ids
    - names
    - types
    - masses
    - charges
    - resids
    - resnames
    - segids
    - moltypes
    - bonds
    - angles
    - dihedrals
    - impropers

    .. _ITP: http://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#molecule-itp-file
    """
    format = 'ITP'

    def parse(self, **kwargs):
        """Parse ITP file into Topology

        Returns
        -------
        MDAnalysis *Topology* object
        """
        ids = []
        types = []
        resids = []
        resnames = []
        names = []
        chargegroups = []
        charges = []
        masses = []

        atom_lists = [ids, types, resids, resnames, names, chargegroups, charges, masses]

        bonds = []
        bondtypes = []
        angles = []
        angletypes = []
        dihedrals = []
        dihedraltypes = []
        impropers = []
        impropertypes = []

        segid = 'SYSTEM'
        molnum = 1

        section = None

        # Open and check itp validity
        with openany(self.filename) as itpfile:
            for line in itpfile:
                line = line.split(';')[0].strip()  # ; is for comments
                if not line:  # Skip line if empty
                    continue
                if line.startswith('#'): #ignore include, ifdefs, etc
                    continue
                
                if '[' in line and ']' in line:
                    section = line.split('[')[1].split(']')[0].strip()
                    continue

                if not section:
                    continue
                
                if section == 'moleculetype':
                    segid = line.split()[0]

                # atoms
                elif section == 'atoms':
                    values = line.split()
                    for lst in atom_lists:
                        try:
                            lst.append(values.pop(0))
                        except IndexError: # ran out of values
                            lst.append('')
                
                elif section == 'bonds':
                    values = line.split()
                    bonds.append(values[:2])
                    bondtypes.append(values[2])
                
                elif section == 'angles':
                    values = line.split()
                    angles.append(values[:3])
                    angletypes.append(values[3])
                
                elif section == 'dihedrals':
                    # funct == 1: proper
                    # funct == 2: improper
                    values = line.split()
                    if values[4] == '1':
                        dihedrals.append(values[:4])
                        dihedraltypes.append(values[4])
                    elif values[4] == '2':
                        impropers.append(values[:4])
                        impropertypes.append(values[4])
        
        attrs = []
        
        # atom stuff
        for vals, Attr, dtype in (
            (ids, Atomids, np.int32),
            (types, Atomtypes, object),
            (names, Atomnames, object),
            (chargegroups, Chargegroups, np.int32),
            (charges, Charges, np.float32),
        ):
            if all(vals):
                attrs.append(Attr(np.array(vals, dtype=dtype)))
        
        if not all(masses):
            masses = guess_masses(types)
            attrs.append(Masses(masses, guessed=True))
        else:
            attrs.append(Masses(np.array(masses, dtype=np.float64)))

        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)
        residx, (resids, resnames) = change_squash((resids,), (resids, resnames))
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(Resnames(resnames))
        
        n_atoms = len(ids)
        n_residues = len(resids)
        n_segments = 1
        attrs.append(Segids(np.array([segid], dtype=object)))
        attrs.append(Moltypes(np.array([segid]*n_residues)))
        segidx = None

        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)
        
        # connectivity
        for vals, Attr, attrname, atype in (
            (bonds, Bonds, 'bonds', bondtypes),
            (angles, Angles, 'angles', angletypes),
            (dihedrals, Dihedrals, 'dihedrals', dihedraltypes),
            (impropers, Impropers, 'impropers', impropertypes)
        ):
            tattr = Attr(set(tuple(map(ids.index, x)) for x in vals),
                         types=atype)
            top.add_TopologyAttr(tattr)

        return top


