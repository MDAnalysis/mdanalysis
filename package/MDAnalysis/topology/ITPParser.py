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
from __future__ import absolute_import
from collections import defaultdict

import logging
import numpy as np

from ..lib.util import openany
from . import guessers
from .base import TopologyReaderBase, change_squash, squash_identical
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

        bonds = defaultdict(list)
        angles = defaultdict(list)
        dihedrals = defaultdict(list)
        impropers = defaultdict(list)

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
                    funct = int(values[2])
                    if funct in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10):
                        bonds[tuple(values[:2])].append(funct)
                
                elif section == 'constraints':
                    values = line.split()
                    funct = int(values[2])
                    if funct in (1, 2):
                        bonds[tuple(values[:2])].append(funct)

                elif section == 'settles':
                    values = line.split()
                    funct = int(values[2])
                    if funct in (1,):
                        bonds[tuple(values[:2])].append(funct)
                
                elif section == 'angles':
                    values = line.split()
                    funct = int(values[3])
                    if funct in (1, 2, 3, 4, 5, 6, 8, 10):
                        angles[tuple(values[:3])].append(funct)

                
                elif section == 'dihedrals':
                    # funct == 1: proper
                    # funct == 2: improper
                    values = line.split()
                    funct = int(values[4])
                    if funct in (1, 3, 5, 8, 9, 10, 11):
                        dihedrals[tuple(values[:4])].append(funct)
                    elif funct in (2, 4):
                        impropers[tuple(values[:4])].append(funct)
        
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
            masses = guessers.guess_masses(types)
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
        for dct, Attr, attrname in (
            (bonds, Bonds, 'bonds'),
            (angles, Angles, 'angles'),
            (dihedrals, Dihedrals, 'dihedrals'),
            (impropers, Impropers, 'impropers')
        ):
            vals, types = zip(*list(dct.items()))

            indices = tuple(map(ids.index, x) for x in vals)
            types = [squash_identical(t) for t in types]
            
            tattr = Attr(indices, types=types)
            top.add_TopologyAttr(tattr)

        return top


