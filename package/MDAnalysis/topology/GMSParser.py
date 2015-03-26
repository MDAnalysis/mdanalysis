# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
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
GMS Topology Parser
===================

Reads an GAMESS output file and pulls element information from it.
Symmetrical assembly is readed (not symmetry element!). 
Atom names becomes as read from GAMESS section.
All information about residues and segments won't be populated.

.. versionadded:: 0.9.1
"""
from __future__ import absolute_import

from ..core.AtomGroup import Atom
from ..core.util import openany
from .core import get_atom_mass, guess_atom_charge, guess_atom_element
from .base import TopologyReader
import re


class GMSParser(TopologyReader):

    def parse(self):
        with openany(self.filename, 'r') as inf:
            natoms = 0
            while True:
                line = inf.readline()
                if not line:
                    raise EOFError
                if re.match(r'^\s+ATOM\s+ATOMIC\s+COORDINATES\s*\(BOHR\).*',\
                        line):
                    break
            line = inf.readline() # skip
            atoms = []
            i = 0
            segid = "SYSTEM"
            resid = 1
            resname = "SYSTEM"
            while True:
                line = inf.readline()
                _m = re.match(\
r'^\s*([A-Za-z_][A-Za-z_0-9]*)\s+([0-9]+\.[0-9]+)\s+(\-?[0-9]+\.[0-9]+)\s+(\-?[0-9]+\.[0-9]+)\s+(\-?[0-9]+\.[0-9]+).*',
                        line)
                if _m == None:
                    break
                name = _m.group(1)
                elem = int(float(_m.group(2)))
                charge = guess_atom_charge(name)
                mass = get_atom_mass(elem)
                #TODO: may be use coordinates info from _m.group(3-5) ??
                at = Atom(i, name, elem, resname, resid,
                          segid, mass, charge)
                atoms.append(at)
                i += 1

        struc = {"_atoms": atoms}

        return struc
