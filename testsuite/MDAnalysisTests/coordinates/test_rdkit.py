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
from __future__ import absolute_import

import pytest
import MDAnalysis as mda
import numpy as np
from numpy.testing import (assert_equal,
                           assert_almost_equal)

from MDAnalysisTests.datafiles import mol2_molecule

Chem = pytest.importorskip("rdkit.Chem")


class TestRDKitReader(object):
    @pytest.fixture
    def mol(self):
        return Chem.MolFromMol2File(mol2_molecule, removeHs=False)

    @pytest.fixture
    def universe(self, mol):
        return mda.Universe(mol)

    def test_coordinates(self, mol, universe):
        expected = np.array([
            conf.GetPositions() for conf in mol.GetConformers()], 
            dtype=np.float32)
        assert universe.trajectory.n_frames == 1
        assert_equal(expected, universe.trajectory.coordinate_array)

    def test_no_coordinates(self):
        pass

    def test_compare_mol2reader(self, universe):
        mol2 = mda.Universe(mol2_molecule)
        assert universe.trajectory.n_frames == mol2.trajectory.n_frames
        assert_equal(universe.trajectory.ts.positions, 
                     mol2.trajectory.ts.positions)
        
