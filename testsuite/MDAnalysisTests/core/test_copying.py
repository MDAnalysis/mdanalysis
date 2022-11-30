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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import numpy as np
from numpy.testing import assert_equal
import pytest

from MDAnalysisTests.datafiles import (
    PSF, DCD, PDB_small
)

import MDAnalysis as mda
from MDAnalysis.core import topology
from MDAnalysis.core import topologyattrs as ta


@pytest.fixture()
def refTT():
    ref = topology.TransTable(
        9, 6, 3,
        atom_resindex=np.array([0, 0, 1, 1, 2, 2, 3, 4, 5]),
        residue_segindex=np.array([0, 1, 2, 0, 1, 1])
    )
    return ref


class TestTransTableCopy(object):
    def test_size(self, refTT):
        new = refTT.copy()
        assert new.n_atoms == refTT.n_atoms
        assert new.n_residues == refTT.n_residues
        assert new.n_segments == refTT.n_segments

    def test_size_independent(self, refTT):
        # check changing 
        new = refTT.copy()
        old = refTT.n_atoms
        refTT.n_atoms = -10
        assert new.n_atoms == old

    @pytest.mark.parametrize('attr', ['_AR', 'RA', '_RS', 'SR'])
    def test_AR(self, refTT, attr):
        new = refTT.copy()
        ref = getattr(refTT, attr)
        other = getattr(new, attr)
        # arrays are dtype='object'
        for a, b in zip(ref, other):
            assert_equal(a, b)

    @pytest.mark.parametrize('attr', ['_AR', 'RA', '_RS', 'SR'])
    def test_AR_independent(self, refTT, attr):
        new = refTT.copy()
        ref = getattr(refTT, attr)
        other = getattr(new, attr)
        # check that arrays don't share memory instead maybe?
        assert ref is not other

    def test_move_atom(self, refTT):
        new = refTT.copy()
        # move atom 0 to residue 3
        refTT.move_atom(0, 3)
        assert refTT.atoms2residues(0) == 3
        assert new.atoms2residues(0) == 0

    def test_move_residue(self, refTT):
        new = refTT.copy()
        refTT.move_residue(1, 2)
        assert refTT.residues2segments(1) == 2
        assert new.residues2segments(1) == 1


TA_FILLER = {
    object: np.array(['dave', 'steve', 'hugo'], dtype=object),
    int: np.array([5, 4, 6]),
    float: np.array([15.4, 5.7, 22.2]),
    'record': np.array(['ATOM', 'ATOM', 'HETATM'], dtype='object'),
    'bond': [(0, 1), (1, 2), (5, 6)],
    'angles': [(0, 1, 2), (1, 2, 3), (4, 5, 6)],
    'dihe': [(0, 1, 2, 3), (1, 2, 3, 4), (5, 6, 7, 8)],
}

@pytest.fixture(params=[
    (ta.Atomids, int),
    (ta.Atomnames, object),
    (ta.Atomtypes, object),
    (ta.Elements, object),
    (ta.Radii, float),
    (ta.RecordTypes, 'record'),
    (ta.ChainIDs, object),
    (ta.Tempfactors, float),
    (ta.Masses, float),
    (ta.Charges, float),
    (ta.Occupancies, float),
    (ta.AltLocs, object),
    (ta.Resids, int),
    (ta.Resnames, object),
    (ta.Resnums, int),
    (ta.ICodes, object),
    (ta.Segids, object),
    (ta.Bonds, 'bond'),
    (ta.Angles, 'angles'),
    (ta.Dihedrals, 'dihe'),
    (ta.Impropers, 'dihe'),
])
def refTA(request):
    cls, dt = request.param
    return cls(TA_FILLER[dt])


def test_copy_attr(refTA):
    new = refTA.copy()

    assert_equal(new.values, refTA.values)
    assert new.values is not refTA.values


@pytest.fixture()
def refTop():
    return topology.Topology(
        3, 2, 2,
        attrs = [
            ta.Atomnames(TA_FILLER[object]),
            ta.Masses(TA_FILLER[float]),
            ta.Resids(TA_FILLER[int]),
            ta.Bonds(TA_FILLER['bond']),
        ],
        atom_resindex=np.array([0, 0, 1]),
        residue_segindex=np.array([0, 1])
    )

def test_topology_copy_n_attrs(refTop):
    new = refTop.copy()
    assert len(new.attrs) == 7  # 4 + 3 indices

@pytest.mark.parametrize('attr', [
    'names',
    'masses',
    'resids',
    'bonds',
    'tt',
])
def test_topology_copy_unique_attrs(refTop, attr):
    new = refTop.copy()
    assert getattr(refTop, attr) is not getattr(new, attr)



@pytest.fixture(scope='module')
def refUniverse():
    return mda.Universe(PSF, DCD)

class TestCopyUniverse(object):
    def test_universe_copy(self, refUniverse):
        new = refUniverse.copy()

        assert new is not refUniverse
        assert len(new.atoms) == len(refUniverse.atoms)

    def test_positions(self, refUniverse):
        new = refUniverse.copy()
        
        assert_equal(new.atoms.positions, refUniverse.atoms.positions)

    def test_change_positions(self, refUniverse):
        # check that coordinates act independently
        new = refUniverse.copy()

        previous = new.atoms[0].position.copy()
        refUniverse.atoms[0].position = 1, 2, 3

        assert_equal(new.atoms[0].position, previous)
        assert_equal(refUniverse.atoms[0].position, [1, 2, 3])
        
    def test_topology(self, refUniverse):
        new = refUniverse.copy()

        assert_equal(new.atoms.names, refUniverse.atoms.names)

    def test_change_topology(self, refUniverse):
        new = refUniverse.copy()

        previous = new.atoms[0].name
        refUniverse.atoms[0].name = 'newname'

        assert new.atoms[0].name == previous
        assert refUniverse.atoms[0].name == 'newname'

def test_pdb_copy():
    u = mda.Universe(PDB_small)

    u2 = u.copy()

    assert_equal(u.atoms.record_types, u2.atoms.record_types)
