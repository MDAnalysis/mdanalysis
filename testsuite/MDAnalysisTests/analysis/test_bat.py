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

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal
import pytest

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (PSF, DCD, mol2_comments_header, XYZ_mini,
                                       BATArray)
from MDAnalysis.analysis.bat import BAT


class TestBAT(object):
    @pytest.fixture()
    def selected_residues(self):
        u = mda.Universe(mol2_comments_header)
        ag = u.select_atoms("all")
        return ag

    @pytest.fixture()
    def bat(self, selected_residues):
        R = BAT(selected_residues)
        R.run()
        return R.bat

    @pytest.fixture
    def bat_npz(self, tmpdir, selected_residues):
        filename = str(tmpdir / 'test_bat_IO.npy')
        R = BAT(selected_residues)
        R.run()
        R.save(filename)
        return filename

    def test_bat_root_selection(self, selected_residues):
        R = BAT(selected_residues)
        assert_equal(R._root.indices, [8, 2, 1],
                     err_msg="error: incorrect root atoms selected")

    def test_bat_number_of_frames(self, bat):
        assert_equal(len(bat),
                     2,
                     err_msg="error: list is not length of trajectory")

    def test_bat_coordinates(self, bat):
        test_bat = np.load(BATArray)
        assert_almost_equal(
            bat,
            test_bat,
            5,
            err_msg="error: BAT coordinates should match test values")

    def test_bat_coordinates_single_frame(self, selected_residues):
        bat = BAT(selected_residues).run(start=1, stop=2).bat
        test_bat = [np.load(BATArray)[1]]
        assert_almost_equal(
            bat,
            test_bat,
            5,
            err_msg="error: BAT coordinates should match test values")

    def test_bat_reconstruction(self, selected_residues, bat):
        R = BAT(selected_residues)
        XYZ = R.Cartesian(bat[0])
        assert_almost_equal(XYZ, selected_residues.positions, 5,
            err_msg="error: Reconstructed Cartesian coordinates " + \
                    "don't match original")

    def test_bat_IO(self, bat_npz, selected_residues, bat):
        R2 = BAT(selected_residues, filename=bat_npz)
        test_bat = R2.bat
        assert_almost_equal(
            bat,
            test_bat,
            5,
            err_msg="error: Loaded BAT coordinates should match test values")

    def test_bat_nobonds(self):
        u = mda.Universe(XYZ_mini)
        with pytest.raises(AttributeError):
            Z = BAT(u.atoms)

    def test_bat_bad_initial_atom(self, selected_residues):
        with pytest.raises(ValueError):
            R = BAT(selected_residues, initial_atom = selected_residues[0])

    def test_bat_disconnected_atom_group(self):
        u = mda.Universe(PSF, DCD)
        selected_residues = u.select_atoms("resid 1-3") + \
            u.select_atoms("resid 5-7")
        with pytest.raises(ValueError):
            R = BAT(selected_residues)

    def test_bat_incorrect_dims(self, bat_npz):
        u = mda.Universe(PSF, DCD)
        selected_residues = u.select_atoms("resid 1-3")
        with pytest.raises(ValueError):
            R = BAT(selected_residues, filename=bat_npz)
