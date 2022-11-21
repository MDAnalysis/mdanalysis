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
import itertools
import numpy as np
import pytest

import MDAnalysis as mda
from MDAnalysis.topology.MinimalParser import MinimalParser

from MDAnalysisTests.datafiles import (
    DCD,
    INPCRD,
    LAMMPSdcd2,
    NCDF,
    TRJ,
    TRJncdf,
    TRR,
    TRZ,
    XTC,
)


working_readers = pytest.mark.parametrize(
    'filename,expected_n_atoms', [
        (DCD, 3341),
        (INPCRD, 5),
        (LAMMPSdcd2, 12421),
        (NCDF, 2661),
        (TRR, 47681),
        (XTC, 47681),
        (np.zeros((1, 10, 3)), 10),  # memory reader default
    ])

@working_readers
def test_minimal_parser(filename, expected_n_atoms):
    with MinimalParser(filename) as p:
        top = p.parse()
    assert top.n_atoms == expected_n_atoms


@working_readers
def test_universe_with_minimal(filename, expected_n_atoms):
    u = mda.Universe(filename, to_guess=())

    assert len(u.atoms) == expected_n_atoms


nonworking_readers = pytest.mark.parametrize('filename,n_atoms', [
    (TRJ, 252),
    (TRJncdf, 2661),
    (TRZ, 8184),
])

@nonworking_readers
def test_minimal_parser_fail(filename,n_atoms):
    with MinimalParser(filename) as p:
        with pytest.raises(NotImplementedError):
            p.parse()


@nonworking_readers
def test_minimal_n_atoms_kwarg(filename, n_atoms):
    # test that these can load when we supply the number of atoms
    u = mda.Universe(filename, n_atoms=n_atoms, to_guess=())

    assert len(u.atoms) == n_atoms


def memory_possibilities():
    # iterate over all possible shapes for a MemoryReader array
    # number of frames, atoms and coordinates
    n = {'f': 1, 'a': 10, 'c': 3}
    for permutation in itertools.permutations('fac', 3):
        order = ''.join(permutation)
        array = np.zeros([n[val] for val in permutation])

        yield array, order


memory_reader = pytest.mark.parametrize('array,order', list(memory_possibilities()))

@memory_reader
def test_memory_minimal_parser(array, order):
    with MinimalParser(array) as p:
        top = p.parse(order=order)
    assert top.n_atoms == 10

@memory_reader
def test_memory_universe(array, order):
    u = mda.Universe(array, order=order, to_guess=())

    assert len(u.atoms) == 10
