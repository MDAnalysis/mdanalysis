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
import warnings

import pytest
import MDAnalysis as mda
import numpy as np
from numpy.testing import (assert_equal,
                           assert_almost_equal)

from MDAnalysisTests.datafiles import mol2_molecule

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")

def mol2_mol():
    return Chem.MolFromMol2File(mol2_molecule, removeHs=False)

def smiles_mol():
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=3)
    return mol

class TestRDKitReader(object):
    @pytest.mark.parametrize("rdmol, n_frames", [
        (mol2_mol(), 1),
        (smiles_mol(), 3),
    ])
    def test_coordinates(self, rdmol, n_frames):
        universe = mda.Universe(rdmol)
        assert universe.trajectory.n_frames == n_frames
        expected = np.array([
            conf.GetPositions() for conf in rdmol.GetConformers()], 
            dtype=np.float32)
        assert_equal(expected, universe.trajectory.coordinate_array)

    def test_no_coordinates(self):
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            u = mda.Universe.from_smiles("CCO", generate_coordinates=False)
            # Verify the warning
            assert len(w) == 1
            assert "No coordinates found" in str(
                w[-1].message)
        expected = np.empty((1,u.atoms.n_atoms,3), dtype=np.float32)
        expected[:] = np.nan
        assert_equal(u.trajectory.coordinate_array, expected)

    def test_compare_mol2reader(self):
        universe = mda.Universe(mol2_mol())
        mol2 = mda.Universe(mol2_molecule)
        assert universe.trajectory.n_frames == mol2.trajectory.n_frames
        assert_equal(universe.trajectory.ts.positions, 
                     mol2.trajectory.ts.positions)
        
