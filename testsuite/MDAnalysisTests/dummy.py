# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
"""Mock Universe and Topology generated from scratch with default values

"""

import numpy as np
import string
import itertools


# Elements of the standard mock Universe
# 5 atoms per residue, 5 residues per segment
_N_ATOMS = 125
_N_RESIDUES = 25
_N_SEGMENTS = 5
_ATOMS_PER_RES = _N_ATOMS // _N_RESIDUES
_RESIDUES_PER_SEG = _N_RESIDUES // _N_SEGMENTS


def make_Universe(extras=None, size=None,
                  trajectory=False, velocities=False, forces=False):
    """Make a dummy reference Universe

    Allows the construction of arbitrary-sized Universes. Suitable for
    the generation of structures for output.

    Preferable for testing core components because:
      * minimises dependencies within the package
      * very fast compared to a "real" Universe

    Parameters
    ----------
    extras : list of strings, optional
      extra attributes to add to Universe:
      u = make_Universe(('masses', 'charges'))
      Creates a lightweight Universe with only masses and charges.
    size : tuple of int, optional
      number of elements of the Universe (n_atoms, n_residues, n_segments)
      defaults to (125, 25, 5)
    trajectory : bool, optional
      create a fake Reader object attached to Universe, default False
    velocities : bool, optional
      if the fake Reader provides velocities, default False
    force : bool, optional
      if the fake Reader provides forces, default False

    Returns
    -------
    MDAnalysis.core.universe.Universe object

    """
    import MDAnalysis as mda

    if size is None:
        size = _N_ATOMS, _N_RESIDUES, _N_SEGMENTS

    n_atoms, n_residues, n_segments = size
    u = mda.Universe.empty(
        # topology things
        n_atoms=n_atoms,
        n_residues=n_residues,
        n_segments=n_segments,
        atom_resindex=np.repeat(
            np.arange(n_residues), n_atoms // n_residues),
        residue_segindex=np.repeat(
            np.arange(n_segments), n_residues // n_segments),
        # trajectory things
        trajectory=trajectory,
        velocities=velocities,
        forces=forces,
    )
    if extras is None:
        extras = []
    for ex in extras:
        u.add_TopologyAttr(ex, values=_MENU[ex](size))

    if trajectory:
        ts = u.trajectory.ts
        ts.positions = np.arange(3 * n_atoms).reshape(n_atoms, 3)
        ts.frame = 0
        if velocities:
            ts.velocities = np.arange(3 * n_atoms).reshape(n_atoms, 3) + 100
        if forces:
            ts.forces = np.arange(3 * n_atoms).reshape(n_atoms, 3) + 10000

    return u

def make_altLocs(size):
    """AltLocs cycling through A B C D E"""
    na, nr, ns = size
    alts = itertools.cycle(('A', 'B', 'C', 'D', 'E'))
    return np.array(['{}'.format(next(alts)) for _ in range(na)],
                    dtype=object)


def make_bfactors(size):
    na, nr, ns = size
    return np.tile(np.array([1.0, 2, 3, 4, 5]), nr)


def make_tempfactors(size):
    na, nr, ns = size
    return np.tile(np.array([1.0, 2, 3, 4, 5]), nr)


def make_charges(size):
    """Atom charges (-1.5, -0.5, 0.0, 0.5, 1.5) repeated"""
    na, nr, ns = size
    charges = itertools.cycle([-1.5, -0.5, 0.0, 0.5, 1.5])
    return np.array([next(charges) for _ in range(na)])

def make_resnames(size):
    """Creates residues named RsA RsB ... """
    na, nr, ns = size
    return np.array(['Rs{}'.format(string.ascii_uppercase[i])
                     for i in range(nr)], dtype=object)

def make_segids(size):
    """Segids SegA -> SegY"""
    na, nr, ns = size
    return np.array(['Seg{}'.format(string.ascii_uppercase[i])
                     for i in range(ns)], dtype=object)

def make_types(size):
    """Atoms are given types TypeA -> TypeE on a loop"""
    na, nr, ns = size
    types = itertools.cycle(string.ascii_uppercase[:5])
    return np.array(['Type{}'.format(next(types)) for _ in range(na)],
                    dtype=object)

def make_names(size):
    """Atom names AAA -> ZZZ (all unique)"""
    na, nr, ns = size
    # produces, AAA, AAB, AAC, ABA etc
    names = itertools.product(*[string.ascii_uppercase] * 3)
    return np.array(['{}'.format(''.join(next(names))) for _ in range(na)],
                    dtype=object)

def make_occupancies(size):
    na, nr, ns = size
    return np.tile(np.array([1.0, 2, 3, 4, 5]), nr)

def make_radii(size):
    na, nr, ns = size
    return np.tile(np.array([1.0, 2, 3, 4, 5]), nr)

def make_serials(size):
    """Serials go from 10 to size+10"""
    na, nr, ns = size
    return np.arange(na) + 10

def make_masses(size):
    """Atom masses (5.1, 4.2, 3.3, 1.5, 0.5) repeated"""
    na, nr, ns = size
    masses = itertools.cycle([5.1, 4.2, 3.3, 1.5, 0.5])
    return np.array([next(masses) for _ in range(na)])

def make_resnums(size):
    """Resnums 1 and upwards"""
    na, nr, ns = size
    return np.arange(nr, dtype=np.int64) + 1

def make_resids(size):
    """Resids 1 and upwards"""
    na, nr, ns = size
    return np.arange(nr, dtype=np.int64) + 1


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
