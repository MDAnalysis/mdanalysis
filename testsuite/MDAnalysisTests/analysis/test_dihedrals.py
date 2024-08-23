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
import numpy as np
from numpy.testing import assert_equal, assert_allclose
import matplotlib
import pytest

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (GRO, XTC, TPR, DihedralArray,
                                       DihedralsArray, RamaArray, GLYRamaArray,
                                       JaninArray, LYSJaninArray, PDB_rama,
                                       PDB_janin)
import MDAnalysis.analysis.dihedrals
from MDAnalysis.analysis.dihedrals import Dihedral, Ramachandran, Janin


class TestDihedral(object):

    @pytest.fixture()
    def atomgroup(self):
        u = mda.Universe(GRO, XTC)
        ag = u.select_atoms("(resid 4 and name N CA C) or (resid 5 and name N)")
        return ag


    def test_dihedral(self, atomgroup, client_Dihedral):
        # client_Dihedral is defined in testsuite/analysis/conftest.py
        # among with other testing fixtures. During testing, it will
        # collect all possible backends and reasonable number of workers
        # for a given AnalysisBase subclass, and extend the tests
        # to run with all of them.

        dihedral = Dihedral([atomgroup]).run(**client_Dihedral)
        test_dihedral = np.load(DihedralArray)

        assert_allclose(dihedral.results.angles, test_dihedral, rtol=0, atol=1.5e-5,
                            err_msg="error: dihedral angles should "
                            "match test values")

    def test_dihedral_single_frame(self, atomgroup):
        dihedral = Dihedral([atomgroup]).run(start=5, stop=6)
        test_dihedral = [np.load(DihedralArray)[5]]

        assert_allclose(dihedral.results.angles, test_dihedral, rtol=0, atol=1.5e-5,
                            err_msg="error: dihedral angles should "
                            "match test vales")

    def test_atomgroup_list(self, atomgroup, client_Dihedral):
        dihedral = Dihedral([atomgroup, atomgroup]).run(**client_Dihedral)
        test_dihedral = np.load(DihedralsArray)

        assert_allclose(dihedral.results.angles, test_dihedral, rtol=0, atol=1.5e-5,
                            err_msg="error: dihedral angles should "
                            "match test values")

    def test_enough_atoms(self, atomgroup):
        with pytest.raises(ValueError):
            dihedral = Dihedral([atomgroup[:2]]).run()

    def test_dihedral_attr_warning(self, atomgroup):
        dihedral = Dihedral([atomgroup]).run(stop=2)

        wmsg = "The `angle` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(dihedral.angles, dihedral.results.angles)


class TestRamachandran(object):

    @pytest.fixture()
    def universe(self):
        return mda.Universe(GRO, XTC)

    @pytest.fixture()
    def rama_ref_array(self):
        return np.load(RamaArray)

    def test_ramachandran(self, universe, rama_ref_array, client_Ramachandran):
        rama = Ramachandran(universe.select_atoms("protein")).run(
            **client_Ramachandran)

        assert_allclose(rama.results.angles, rama_ref_array, rtol=0, atol=1.5e-5,
                            err_msg="error: dihedral angles should "
                            "match test values")

    def test_ramachandran_single_frame(self, universe, rama_ref_array):
        rama = Ramachandran(universe.select_atoms("protein")).run(
            start=5, stop=6)

        assert_allclose(rama.results.angles[0], rama_ref_array[5], rtol=0, atol=1.5e-5,
                            err_msg="error: dihedral angles should "
                            "match test values")

    def test_ramachandran_residue_selections(self, universe, client_Ramachandran):
        rama = Ramachandran(universe.select_atoms("resname GLY")).run(
            **client_Ramachandran)
        test_rama = np.load(GLYRamaArray)

        assert_allclose(rama.results.angles, test_rama, rtol=0, atol=1.5e-5,
                            err_msg="error: dihedral angles should "
                            "match test values")

    def test_outside_protein_length(self, universe):
        with pytest.raises(ValueError):
            rama = Ramachandran(universe.select_atoms("resid 220"),
                                check_protein=True).run()

    def test_outside_protein_unchecked(self, universe):
        rama = Ramachandran(universe.select_atoms("resid 220"),
                            check_protein=False).run()

    def test_protein_ends(self, universe):
        with pytest.warns(UserWarning) as record:
            rama = Ramachandran(universe.select_atoms("protein"),
                                check_protein=True).run()
        assert len(record) == 1

    def test_None_removal(self):
        with pytest.warns(UserWarning):
            u = mda.Universe(PDB_rama)
            rama = Ramachandran(u.select_atoms("protein").residues[1:-1])

    def test_plot(self, universe):
        ax = Ramachandran(universe.select_atoms("resid 5-10")).run().plot(ref=True)
        assert isinstance(ax, matplotlib.axes.Axes), \
            "Ramachandran.plot() did not return and Axes instance"

    def test_ramachandran_attr_warning(self, universe):
        rama = Ramachandran(universe.select_atoms("protein")).run(stop=2)

        wmsg = "The `angle` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(rama.angles, rama.results.angles)


class TestJanin(object):

    @pytest.fixture()
    def universe(self):
        return mda.Universe(GRO, XTC)

    @pytest.fixture()
    def universe_tpr(self):
        return mda.Universe(TPR, XTC)

    @pytest.fixture()
    def janin_ref_array(self):
        return np.load(JaninArray)

    def test_janin(self, universe, janin_ref_array, client_Janin):
        self._test_janin(universe, janin_ref_array, client_Janin)

    def test_janin_tpr(self, universe_tpr, janin_ref_array, client_Janin):
        """Test that CYSH are filtered (#2898)"""
        self._test_janin(universe_tpr, janin_ref_array, client_Janin)

    def _test_janin(self, u, ref_array, client_Janin):
        janin = Janin(u.select_atoms("protein")).run(**client_Janin)

        # Test precision lowered to account for platform differences with osx
        assert_allclose(janin.results.angles, ref_array, rtol=0, atol=1.5e-3,
                            err_msg="error: dihedral angles should "
                            "match test values")

    def test_janin_single_frame(self, universe, janin_ref_array):
        janin = Janin(universe.select_atoms("protein")).run(start=5, stop=6)

        assert_allclose(janin.results.angles[0], janin_ref_array[5], rtol=0, atol=1.5e-3,
                            err_msg="error: dihedral angles should "
                            "match test values")

    def test_janin_residue_selections(self, universe, client_Janin):
        janin = Janin(universe.select_atoms("resname LYS")).run(**client_Janin)
        test_janin = np.load(LYSJaninArray)

        assert_allclose(janin.results.angles, test_janin, rtol=0, atol=1.5e-3,
                            err_msg="error: dihedral angles should "
                            "match test values")

    def test_outside_protein_length(self, universe):
        with pytest.raises(ValueError):
            janin = Janin(universe.select_atoms("resid 220")).run()

    def test_remove_residues(self, universe):
        with pytest.warns(UserWarning):
            janin = Janin(universe.select_atoms("protein")).run()

    def test_atom_selection(self):
        with pytest.raises(ValueError):
            u = mda.Universe(PDB_janin)
            janin = Janin(u.select_atoms("protein and not resname ALA CYS GLY "
                                         "PRO SER THR VAL"))

    def test_plot(self, universe):
        ax = Janin(universe.select_atoms("resid 5-10")).run().plot(ref=True)
        assert isinstance(ax, matplotlib.axes.Axes), \
            "Ramachandran.plot() did not return and Axes instance"

    def test_janin_attr_warning(self, universe):
        janin = Janin(universe.select_atoms("protein")).run(stop=2)

        wmsg = "The `angle` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(janin.angles, janin.results.angles)


# tests for parallelization

@pytest.mark.parametrize(
    "classname,is_parallelizable",
    [
        (MDAnalysis.analysis.dihedrals.Dihedral, True),
        (MDAnalysis.analysis.dihedrals.Ramachandran, True),
        (MDAnalysis.analysis.dihedrals.Janin, True),
    ]
)
def test_class_is_parallelizable(classname, is_parallelizable):
    assert classname._analysis_algorithm_is_parallelizable == is_parallelizable


@pytest.mark.parametrize(
    "classname,backends",
    [
        (MDAnalysis.analysis.dihedrals.Dihedral,
         ('serial', 'multiprocessing', 'dask',)),
        (MDAnalysis.analysis.dihedrals.Ramachandran,
         ('serial', 'multiprocessing', 'dask',)),
        (MDAnalysis.analysis.dihedrals.Janin,
         ('serial', 'multiprocessing', 'dask',)),
    ]
)
def test_supported_backends(classname, backends):
    assert classname.get_supported_backends() == backends
