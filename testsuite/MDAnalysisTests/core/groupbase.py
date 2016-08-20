# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import numpy as np

import MDAnalysis as mda
from MDAnalysis.core.topology import Topology
from MDAnalysis.core import groups

def make_Universe():
    """Make a dummy reference Universe"""
    return mda.Universe(make_topology())

def make_topology():
    """Reference topology system

    125 atoms
    25 residue
    5 segments
    """
    return Topology(125, 25, 5,
                    atom_resindex=np.repeat(np.arange(25), 5),
                    residue_segindex=np.repeat(np.arange(5), 5))
