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
    u = mda.Universe()

    u._topology = make_topology()

    # generate Universe version of each class
    # AG, RG, SG, A, R, S
    u._classes = groups.make_classes()

    # Put Group level stuff from topology into class
    for attr in u._topology.attrs:
        u._process_attr(attr)

    # Generate atoms, residues and segments
    u.atoms = u._classes['atomgroup'](
        np.arange(u._topology.n_atoms), u)
    u.residues = u._classes['residuegroup'](
        np.arange( u._topology.n_residues), u)
    u.segments = u._classes['segmentgroup'](np.arange(
        u._topology.n_segments), u)

    # Update Universe namespace with segids
    for seg in u.segments:
        if hasattr(seg, 'segid'):
            if seg.segid[0].isdigit():
                name = 's' + seg.segid
            else:
                name = seg.segid
            u.__dict__[name] = seg

    return u

def make_topology():
    """Reference topology system

    125 atoms
    25 residue
    5 segments
    """
    return Topology(125, 25, 5,
                    atom_resindex=np.repeat(np.arange(25), 5),
                    residue_segindex=np.repeat(np.arange(5), 5))
