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


Preprocessor variables
----------------------

ITP files are often defined with lines that depend on 
whether a keyword flag is given. For example:

    #ifndef FLEXIBLE

    [ settles ]
    ; i     funct   doh     dhh
    1       1       0.09572 0.15139

    #else

    [ bonds ]
    ; i     j       funct   length  force.c.
    1       2       1       0.09572 502416.0 0.09572        502416.0
    1       3       1       0.09572 502416.0 0.09572        502416.0

    [ angles ]
    ; i     j       k       funct   angle   force.c.
    2       1       3       1       104.52  628.02  104.52  628.02

    #endif

Define these preprocessor variables by passing keyword arguments::

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import ITP_tip5p

    u = mda.Universe(ITP_tip5p, FLEXIBLE=True)

These keyword variables are **case-sensitive**. Currently, 
MDAnalysis only supports boolean variables in if/else conditions.

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
from .base import TopologyReaderBase, change_squash, reduce_singular
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
        self.ids = []
        types = []
        resids = []
        resnames = []
        names = []
        chargegroups = []
        charges = []
        masses = []

        atom_lists = [self.ids, types, resids, resnames, names, chargegroups, charges, masses]

        bonds = defaultdict(list)
        angles = defaultdict(list)
        dihedrals = defaultdict(list)
        impropers = defaultdict(list)

        segid = 'SYSTEM'
        molnum = 1

        section = None
        if_queue = []

        # Open and check itp validity
        with openany(self.filename) as self.itpfile:
            for line in self.itpfile:
                line = line.split(';')[0].strip()  # ; is for comments
                if not line:  # Skip line if empty
                    continue

                if line.startswith('#ifdef'):
                    kw = line.split()[1]
                    if kwargs.get(kw):
                        continue
                    else:
                        self.skip_if()

                if line.startswith('#ifndef'):
                    kw = line.split()[1]
                    if not kwargs.get(kw):
                        continue
                    else:
                        self.skip_if()
                
                # Only gets triggered when if condition is true
                if line.startswith('#else'):
                    self.skip_else()

                if line.startswith('#'): #ignore include, endif, etc
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
                    self.add_param(line, bonds, n_funct=2, funct_values=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
                
                elif section == 'constraints':
                    self.add_param(line, bonds, n_funct=2, funct_values=(1, 2))

                elif section == 'settles':
                    # [ settles ] is a triangular constraint for 
                    # water molecules.
                    # In ITP files this is defined with only the 
                    # oxygen atom index. The next two atoms are 
                    # assumed to be hydrogens. Unlike TPRParser,  
                    # the manual only lists this format (as of 2019).
                    # These are treated as 2 bonds and 1 angle.
                    oxygen, funct, doh, dhh = line.split()
                    try:
                        base = self.index_ids([oxygen])[0]
                    except ValueError:
                        pass
                    else:
                        bonds[(base, base+1)].append("settles")
                        bonds[(base, base+2)].append("settles")
                        angles[(base+1, base, base+2)].append("settles")
                
                elif section == 'angles':
                    self.add_param(line, angles, n_funct=3, funct_values=(1, 2, 3, 4, 5, 6, 8, 10))
                
                elif section == 'dihedrals':
                    # funct == 1: proper
                    # funct == 2: improper

                    self.add_param(line, dihedrals, n_funct=4, funct_values=(1, 3, 5, 8, 9, 10, 11))
                    self.add_param(line, impropers, n_funct=4, funct_values=(2, 4))
                    
        
        attrs = []
        
        # atom stuff
        for vals, Attr, dtype in (
            (self.ids, Atomids, np.int32),
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
        
        n_atoms = len(self.ids)
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
            if dct:
                indices, types = zip(*list(dct.items()))
            else:
                indices, types = [], []

            types = [reduce_singular(t) for t in types]
            
            tattr = Attr(indices, types=types)
            top.add_TopologyAttr(tattr)

        return top
    
    def skip_if(self):
        """
        Skip lines in if condition, *inclusive* of the #endif or #else.
        """
        line = next(self.itpfile)
        while not line.startswith('#endif') and not line.startswith('#else'):
            line = next(self.itpfile)
        return line

    def skip_else(self):
        """
        Skip lines in else condition, *inclusive* of #endif.
        """
        line = next(self.itpfile)
        while not line.startswith('#endif'):
            line = next(self.itpfile)
        return line

    def index_ids(self, values):
        """
        Get indices of atom ids (list of strings)
        """
        return tuple(map(self.ids.index, values))
    
    def add_param(self, line, container, n_funct=2, funct_values=[]):
        values = line.split()
        funct = int(values[n_funct])
        if funct in funct_values:
            try:
                container[self.index_ids(values[:n_funct])].append(funct)
            except ValueError:
                pass


