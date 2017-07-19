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
from __future__ import division, absolute_import
from six.moves import zip

import numpy as np
from numpy.testing import assert_equal
import pytest

import MDAnalysis as mda
from MDAnalysis.core.topology import TransTable



@pytest.fixture()
def refTT():
    ref = TransTable(9, 6, 3,
                     atom_resindex=np.array([0, 0, 1, 1, 2, 2, 3, 4, 5]),
                     residue_segindex=np.array([0, 1, 2, 0, 1, 1]))
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

    @pytest.mark.parametrize('attr', ['_AR', '_RA', '_RS', '_SR'])
    def test_AR(self, refTT, attr):
        new = refTT.copy()
        ref = getattr(refTT, attr)
        other = getattr(new, attr)
        # arrays are dtype='object'
        for a, b in zip(ref, other):
            assert_equal(a, b)

    @pytest.mark.parametrize('attr', ['_AR', '_RA', '_RS', '_SR'])
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
