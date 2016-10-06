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
"""Mock Universe and Topology for testing purposes



"""
from __future__ import division

from six.moves import range

import numpy as np
import string
import itertools

import MDAnalysis as mda
from MDAnalysis.core.topology import Topology
from MDAnalysis.coordinates.base import SingleFrameReader
from MDAnalysis.core import groups
import MDAnalysis.core.topologyattrs as ta

# Dimensions of the standard mock Universe
# 5 atoms per residues, 5 residues per segment
_N_ATOMS = 125
_N_RESIDUES = 25
_N_SEGMENTS = 5
_ATOMS_PER_RES = _N_ATOMS // _N_RESIDUES
_RESIDUES_PER_SEG = _N_RESIDUES // _N_SEGMENTS

def make_Universe(*extras):
    """Make a dummy reference Universe"""
    return mda.Universe(make_topology(*extras))


def make_topology(*extras):
    """Reference topology system

    125 atoms
    25 residue
    5 segments
    """
    attrs = [_menu[extra]() for extra in extras]

    return Topology(_N_ATOMS, _N_RESIDUES, _N_SEGMENTS,
                    attrs = attrs,
                    atom_resindex=np.repeat(
                        np.arange(_N_RESIDUES), _ATOMS_PER_RES),
                    residue_segindex=np.repeat(
                        np.arange(_N_SEGMENTS), _RESIDUES_PER_SEG))

def make_altLocs(size=None):
    """AltLocs cycling through A B C D E"""
    if size is None:
        size = _N_ATOMS
    alts = itertools.cycle(('A', 'B', 'C', 'D', 'E'))
    return ta.AltLocs(np.array(['{}'.format(next(alts)) for _ in range(size)],
                               dtype=object))

def make_bfactors(size=None):
    if size is None:
        size = _N_ATOMS
    return ta.Bfactors(np.tile(np.array([1.0, 2, 3, 4, 5]), 25))

def make_charges(size=None):
    """Atom charges (-1.5, -0.5, 0.0, 0.5, 1.5) repeated"""
    if size is None:
        size = _N_ATOMS
    charges = itertools.cycle([-1.5, -0.5, 0.0, 0.5, 1.5])
    return ta.Charges(np.array([next(charges)
                                for _ in range(size)]))

def make_resnames(size=None):
    """Creates residues named RsA RsB ... """
    if size is None:
        size = _N_RESIDUES
    return ta.Resnames(np.array(['Rs{}'.format(string.uppercase[i])
                                 for i in range(size)], dtype=object))

def make_segids(size=None):
    """Segids SegA -> SegY"""
    if size is None:
        size = _N_SEGMENTS
    return ta.Segids(np.array(['Seg{}'.format(string.uppercase[i])
                               for i in range(size)], dtype=object))

def make_types(size=None):
    """Atoms are given types TypeA -> TypeE on a loop"""
    if size is None:
        size = _N_ATOMS
    types = itertools.cycle(string.uppercase[:5])
    return ta.Atomtypes(np.array(
        ['Type{}'.format(next(types))
         for _ in range(size)], dtype=object))

def make_names(size=None):
    """Atom names NameAAA -> NameZZZ (all unique)"""
    if size is None:
        size = _N_ATOMS
    # produces, AAA, AAB, AAC, ABA etc
    names = itertools.product(*[string.uppercase] * 3)
    return ta.Atomnames(np.array(
        ['Name{}'.format(''.join(next(names)))
         for _ in range(size)], dtype=object))

def make_occupancies(size=None):
    if size is None:
        size = _N_ATOMS
    return ta.Occupancies(np.tile(np.array([1.0, 2, 3, 4, 5]), 25))

def make_radii(size=None):
    if size is None:
        size = _N_ATOMS
    return ta.Radii(np.tile(np.array([1.0, 2, 3, 4, 5]), 25))

def make_serials(size=None):
    """Serials go from 10 to size+10"""
    if size is None:
        size = _N_ATOMS
    return ta.Atomids(np.arange(size) + 10)

def make_masses(size=None):
    """Atom masses (5.1, 4.2, 3.3, 1.5, 0.5) repeated"""
    if size is None:
        size = _N_ATOMS
    masses = itertools.cycle([5.1, 4.2, 3.3, 1.5, 0.5])
    return ta.Masses(np.array([next(masses)
                               for _ in range(size)]))

def make_resnums(size=None):
    """Resnums 1 and upwards"""
    if size is None:
        size = _N_RESIDUES
    return ta.Resnums(np.arange(size, dtype=np.int64) + 1)

def make_resids(size=None):
    """Resids 1 and upwards"""
    if size is None:
        size = _N_RESIDUES
    return ta.Resids(np.arange(size, dtype=np.int64) + 1)

# Available extra TopologyAttrs to a dummy Universe
_menu = {
    # Atoms
    'altLocs': make_altLocs,
    'bfactors': make_bfactors,
    'charges': make_charges,
    'names': make_names,
    'occupancies': make_occupancies,
    'radii': make_radii,
    'serials': make_serials,
    'types': make_types,
    'masses': make_masses,
    # Residues
    'resnames': make_resnames,
    'resnums': make_resnums,
    'resids': make_resids,
    # Segments
    'segids': make_segids,
}

class FakeReader(SingleFrameReader):
    def __init__(self, n_atoms=None):
        self.n_atoms = n_atoms if not n_atoms is None else _N_ATOMS
        self.filename = 'FakeReader'
        self.n_frames = 1
        self._read_first_frame()

    def _read_first_frame(self):
        self.ts = self._Timestep(self.n_atoms)
        self.ts.positions = np.arange(3 * self.n_atoms).reshape(self.n_atoms, 3)
        self.ts.frame = 0

