# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""Mock Universe and Topology for testing purposes

Avoids import of MDAnalysis at base level because of Issue #344

"""
from __future__ import absolute_import
from __future__ import division

from six.moves import range

import numpy as np
import string
import itertools

# Dimensions of the standard mock Universe
# 5 atoms per residues, 5 residues per segment
_N_ATOMS = 125
_N_RESIDUES = 25
_N_SEGMENTS = 5
_ATOMS_PER_RES = _N_ATOMS // _N_RESIDUES
_RESIDUES_PER_SEG = _N_RESIDUES // _N_SEGMENTS


def make_Universe(extras=None, size=None, trajectory=False, velocities=False, forces=False):
    """Make a dummy reference Universe

    u = make_Universe(('masses', 'charges'))

    Creates a lightweight Universe with only masses and charges.

    Preferable for testing core components because:
      * minimises dependencies within the package
      * very fast compared to a "real" Universe

    Parameters
    ----------
    extras : list of strings
      extra attributes to add to Universe
    size : tuple of int
      dimensions of the Universe (n_atoms, n_residues, n_segments)
    trajectory : bool
      create a fake Reader object attached to Universe
    velocities : bool
      if the fake Reader provides velocities
    force : bool
      if the fake Reader provides forces

    Returns
    -------
    MDAnalysis.core.universe.Universe object

    """
    import MDAnalysis as mda

    u = mda.Universe(make_topology(extras, size=size))
    if trajectory:
        u.trajectory = make_FakeReader(len(u.atoms), velocities, forces)
    return u


def make_topology(extras=None, size=None):
    """Reference topology system

    Parameters
    ----------
    extras : list of strings
      attributes to add to the Universe
    size : tuple of natoms, nres, nseg
      by default: (125 atoms, 25 residue, 5 segments)

    Returns
    -------
    MDAnalysis.core.topology.Topology object
    """
    from MDAnalysis.core.topology import Topology

    if extras is None:
        extras = []
    if size is None:
        size = _N_ATOMS, _N_RESIDUES, _N_SEGMENTS
    attrs = [_MENU[extra](size) for extra in extras]

    return Topology(size[0], size[1], size[2],
                    attrs=attrs,
                    atom_resindex=np.repeat(
                        np.arange(size[1]), size[0] // size[1]),
                    residue_segindex=np.repeat(
                        np.arange(size[2]), size[1] // size[2]))

def make_altLocs(size):
    """AltLocs cycling through A B C D E"""
    import MDAnalysis.core.topologyattrs as ta

    na, nr, ns = size
    alts = itertools.cycle(('A', 'B', 'C', 'D', 'E'))
    return ta.AltLocs(np.array(['{}'.format(next(alts)) for _ in range(na)],
                               dtype=object))

def make_bfactors(size):
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    return ta.Bfactors(np.tile(np.array([1.0, 2, 3, 4, 5]), nr))

def make_tempfactors(size):
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    return ta.Tempfactors(np.tile(np.array([1.0, 2, 3, 4, 5]), nr))

def make_charges(size):
    """Atom charges (-1.5, -0.5, 0.0, 0.5, 1.5) repeated"""
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    charges = itertools.cycle([-1.5, -0.5, 0.0, 0.5, 1.5])
    return ta.Charges(np.array([next(charges)
                                for _ in range(na)]))

def make_resnames(size):
    """Creates residues named RsA RsB ... """
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    return ta.Resnames(np.array(['Rs{}'.format(string.ascii_uppercase[i])
                                 for i in range(nr)], dtype=object))

def make_segids(size):
    """Segids SegA -> SegY"""
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    return ta.Segids(np.array(['Seg{}'.format(string.ascii_uppercase[i])
                               for i in range(ns)], dtype=object))

def make_types(size):
    """Atoms are given types TypeA -> TypeE on a loop"""
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    types = itertools.cycle(string.ascii_uppercase[:5])
    return ta.Atomtypes(np.array(
        ['Type{}'.format(next(types))
         for _ in range(na)], dtype=object))

def make_names(size):
    """Atom names AAA -> ZZZ (all unique)"""
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    # produces, AAA, AAB, AAC, ABA etc
    names = itertools.product(*[string.ascii_uppercase] * 3)
    return ta.Atomnames(np.array(
        ['{}'.format(''.join(next(names)))
         for _ in range(na)], dtype=object))

def make_occupancies(size):
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    return ta.Occupancies(np.tile(np.array([1.0, 2, 3, 4, 5]), nr))

def make_radii(size):
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    return ta.Radii(np.tile(np.array([1.0, 2, 3, 4, 5]), nr))

def make_serials(size):
    """Serials go from 10 to size+10"""
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    return ta.Atomids(np.arange(na) + 10)

def make_masses(size):
    """Atom masses (5.1, 4.2, 3.3, 1.5, 0.5) repeated"""
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    masses = itertools.cycle([5.1, 4.2, 3.3, 1.5, 0.5])
    return ta.Masses(np.array([next(masses)
                               for _ in range(na)]))

def make_resnums(size):
    """Resnums 1 and upwards"""
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    return ta.Resnums(np.arange(nr, dtype=np.int64) + 1)

def make_resids(size):
    """Resids 1 and upwards"""
    import MDAnalysis.core.topologyattrs as ta
    na, nr, ns = size
    return ta.Resids(np.arange(nr, dtype=np.int64) + 1)

# Available extra TopologyAttrs to a dummy Universe
_MENU = {
    # Atoms
    'altLocs': make_altLocs,
    'bfactors': make_bfactors,
    'charges': make_charges,
    'names': make_names,
    'occupancies': make_occupancies,
    'radii': make_radii,
    'serials': make_serials,
    'tempfactors': make_tempfactors,
    'types': make_types,
    'masses': make_masses,
    # Residues
    'resnames': make_resnames,
    'resnums': make_resnums,
    'resids': make_resids,
    # Segments
    'segids': make_segids,
}

def make_FakeReader(n_atoms=None, velocities=False, forces=False):
    from MDAnalysis.coordinates.base import SingleFrameReaderBase

    class FakeReader(SingleFrameReaderBase):
        def __init__(self, n_atoms=None, velocities=False, forces=False):
            self.n_atoms = n_atoms if not n_atoms is None else _N_ATOMS
            self.filename = 'FakeReader'
            self.n_frames = 1
            self._read_first_frame(velocities, forces)

        def _read_first_frame(self, velocities=False, forces=False):
            ts = self.ts = self._Timestep(self.n_atoms, positions=True,
                                          velocities=velocities, forces=forces)
            ts.positions = np.arange(3 * self.n_atoms).reshape(self.n_atoms, 3)
            if velocities:
                ts.velocities = np.arange(3 * self.n_atoms).reshape(self.n_atoms, 3) + 100
            if forces:
                ts.forces = np.arange(3 * self.n_atoms).reshape(self.n_atoms, 3) + 10000
            ts.frame = 0

    return FakeReader(n_atoms, velocities, forces)
