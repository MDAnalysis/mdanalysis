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
import pytest
from pathlib import Path
import MDAnalysis as mda
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    ITP,  # GROMACS itp
    ITP_nomass,  # from Automated Topology Builder
    ITP_atomtypes,
    ITP_charges,
    ITP_edited,
    ITP_tip5p,
    ITP_spce,
    GMX_TOP,
    GMX_DIR,
    GMX_TOP_BAD,
    ITP_no_endif,
)


class BaseITP(ParserBase):
    parser = mda.topology.ITPParser.ITPParser
    expected_attrs = ['ids', 'names', 'types', 'masses',
                      'charges', 'chargegroups',
                      'resids', 'resnames',
                      'segids', 'moltypes', 'molnums',
                      'bonds', 'angles', 'dihedrals', 'impropers']
    expected_n_atoms = 63
    expected_n_residues = 10
    expected_n_segments = 1

    expected_n_bonds = 0
    expected_n_angles = 0
    expected_n_dihedrals = 0
    expected_n_impropers = 0

    @pytest.fixture
    def universe(self, filename):
        return mda.Universe(filename)

    def test_bonds_total_counts(self, top):
        assert len(top.bonds.values) == self.expected_n_bonds
    
    def test_angles_total_counts(self, top):
        assert len(top.angles.values) == self.expected_n_angles

    def test_dihedrals_total_counts(self, top):
        assert len(top.dihedrals.values) == self.expected_n_dihedrals
    
    def test_impropers_total_counts(self, top):
        assert len(top.impropers.values) == self.expected_n_impropers


class TestITP(BaseITP):
    ref_filename = ITP

    expected_n_atoms = 63
    expected_n_residues = 10
    expected_n_segments = 1

    expected_n_bonds = 62
    expected_n_angles = 91
    expected_n_dihedrals = 30
    expected_n_impropers = 29
    
    def test_bonds_atom_counts(self, universe):
        assert len(universe.atoms[[0]].bonds) == 3
        assert len(universe.atoms[[42]].bonds) == 1

    def test_bonds_values(self, top):
        vals = top.bonds.values
        for b in ((0, 1), (0, 2), (0, 3), (3, 4)):
            assert b in vals
        
    def test_bonds_type(self, universe):
        assert universe.bonds[0].type == 2

    def test_angles_atom_counts(self, universe):
        assert len(universe.atoms[[0]].angles) == 5
        assert len(universe.atoms[[42]].angles) == 2

    def test_angles_values(self, top):
        vals = top.angles.values
        for b in ((1, 0, 2), (1, 0, 3), (2, 0, 3)):
            assert (b in vals) or (b[::-1] in vals)
    
    def test_angles_type(self, universe):
        assert universe.angles[0].type == 2

    def test_dihedrals_atom_counts(self, universe):
        assert len(universe.atoms[[0]].dihedrals) == 2

    def test_dihedrals_multiple_types(self, universe):
        ag = universe.atoms[[0, 3, 5, 7]]
        dih = universe.dihedrals.atomgroup_intersection(ag, strict=True)[0]
        assert len(dih.type) == 2

    def test_dihedrals_values(self, top):
        vals = top.dihedrals.values
        for b in ((1, 0, 3, 5), (0, 3, 5, 7)):
            assert (b in vals) or (b[::-1] in vals)
    
    def test_dihedrals_type(self, universe):
        assert universe.dihedrals[0].type == (1, 1)

    def test_impropers_atom_counts(self, universe):
        assert len(universe.atoms[[0]].impropers) == 1

    def test_impropers_values(self, top):
        vals = top.impropers.values
        for b in ((3, 0, 5, 4), (5, 3, 7, 6)):
            assert (b in vals) or (b[::-1] in vals)
    
    def test_impropers_type(self, universe):
        assert universe.impropers[0].type == 2


class TestITPNoMass(ParserBase):
    parser = mda.topology.ITPParser.ITPParser
    ref_filename = ITP_nomass
    expected_attrs = ['ids', 'names', 'types', 'masses',
                      'charges', 'chargegroups',
                      'resids', 'resnames',
                      'segids', 'moltypes', 'molnums',
                      'bonds', 'angles', 'dihedrals', 'impropers']
    guessed_attrs = ['masses']
    expected_n_atoms = 60
    expected_n_residues = 1
    expected_n_segments = 1

    @pytest.fixture
    def universe(self, filename):
        return mda.Universe(filename)

    def test_mass_guess(self, universe):
        assert universe.atoms[0].mass not in ('', None)


class TestITPAtomtypes(ParserBase):
    parser = mda.topology.ITPParser.ITPParser
    ref_filename = ITP_atomtypes
    expected_attrs = ['ids', 'names', 'types', 'masses',
                      'charges', 'chargegroups',
                      'resids', 'resnames',
                      'segids', 'moltypes', 'molnums',
                      'bonds', 'angles', 'dihedrals', 'impropers']
    guessed_attrs = ['masses']
    expected_n_atoms = 4
    expected_n_residues = 1
    expected_n_segments = 1

    @pytest.fixture
    def universe(self, filename):
        return mda.Universe(filename)

    def test_charge_parse(self, universe):
        assert_allclose(universe.atoms[0].charge, 4)
        assert_allclose(universe.atoms[1].charge, 1.1)
        assert_allclose(universe.atoms[2].charge, -3.000)
        assert_allclose(universe.atoms[3].charge, 1.)

    def test_mass_parse_or_guess(self, universe):
        # read from [ atoms ] section
        assert_allclose(universe.atoms[0].mass, 8.0)
        # read from [ atomtypes ] section
        assert_allclose(universe.atoms[1].mass, 20.989)
        # read from [ atomtypes ] section
        assert_allclose(universe.atoms[2].mass, 20.989)
        # guessed from element
        assert_allclose(universe.atoms[3].mass, 1.008)


class TestITPCharges(ParserBase):
    parser = mda.topology.ITPParser.ITPParser
    ref_filename = ITP_charges
    expected_attrs = ['ids', 'names', 'types', 'masses',
                      'charges', 'chargegroups',
                      'resids', 'resnames',
                      'segids', 'moltypes', 'molnums',
                      'bonds', 'angles', 'dihedrals', 'impropers']
    guessed_attrs = []
    expected_n_atoms = 9
    expected_n_residues = 3
    expected_n_segments = 1

    @pytest.fixture
    def universe(self, filename):
        return mda.Universe(filename)

    def test_charge_parse(self, universe):
        assert_allclose(universe.atoms[0].charge, -1.0)
        assert_allclose(universe.atoms[1].charge, 0)
        assert_allclose(universe.atoms[2].charge, 0)
        assert_allclose(universe.atoms[3].charge, -1.)

    def test_masses_are_read(self, universe):
        assert_allclose(universe.atoms.masses, [100] * 9)

class TestDifferentDirectivesITP(BaseITP):

    ref_filename = ITP_edited

    expected_n_bonds = 62 + 118
    expected_n_angles = 88
    expected_n_dihedrals = 28
    expected_n_impropers = 29

    def test_no_extra_angles(self, top):
        for a in ((57, 59, 61), (60, 59, 61), (59, 61, 62)):
            assert a not in top.angles.values

    def test_bonds_atom_counts(self, universe):
        assert len(universe.atoms[[0]].bonds) == 5
        assert len(universe.atoms[[42]].bonds) == 5

    def test_dihedrals_atom_counts(self, universe):
        assert len(universe.atoms[[0]].dihedrals) == 1

    def test_dihedrals_identity(self, universe):
        assert universe.dihedrals[0].type == (1, 1)


class TestITPNoKeywords(BaseITP):
    """
    Test reading ITP files *without* defined keywords.

    tip5p.itp has the lines:

        #ifndef HW1_CHARGE
            #define HW1_CHARGE 0.241
        #endif
        
        [ atoms ]
            1       opls_118     1       SOL              OW             1       0
            2       opls_119     1       SOL             HW1             1       HW1_CHARGE
    """
    ref_filename = ITP_tip5p
    expected_n_atoms = 5
    expected_n_residues = 1
    expected_n_segments = 1

    guessed_attrs = ['masses']

    expected_n_bonds = 2
    # FLEXIBLE not set -> SETTLE constraint -> water has no angle
    expected_n_angles = 0
    expected_n_dihedrals = 0
    expected_n_impropers = 0

    def test_whether_settles_types(self, universe):
        for param in list(universe.bonds) + list(universe.angles):
            assert param.type == 'settles'

    def test_bonds_values(self, top):
        vals = top.bonds.values
        for b in [(0, 1), (0, 2)]:
            assert b in vals

    def test_defines(self, top):
        assert_allclose(top.charges.values[1], 0.241)
        assert_allclose(top.charges.values[2], 0.241)

    
class TestITPKeywords(TestITPNoKeywords):
    """
    Test reading ITP files *with* defined keywords.
    """

    expected_n_atoms = 7
    # FLEXIBLE is set -> no SETTLE constraint -> water should have angle
    expected_n_angles = 1

    @pytest.fixture
    def universe(self, filename):
        return mda.Universe(filename, FLEXIBLE=True, EXTRA_ATOMS=True, 
                            HW1_CHARGE=1, HW2_CHARGE=3)

    @pytest.fixture()
    def top(self, filename):
        with self.parser(filename) as p:
            yield p.parse(FLEXIBLE=True, EXTRA_ATOMS=True, 
                          HW1_CHARGE=1, HW2_CHARGE=3)

    def test_whether_settles_types(self, universe):
        for param in list(universe.bonds) + list(universe.angles):
            assert param.type == 1

    def test_angles_values(self, top):
        assert (1, 0, 2) in top.angles.values

    def test_defines(self, top):
        assert_allclose(top.charges.values[1], 1)

    def test_kwargs_overrides_defines(self, top):
        assert_allclose(top.charges.values[2], 3)


class TestNestedIfs(BaseITP):
    """
    Test reading ITP files with nested ifdef/ifndef conditions.
    """
    ref_filename = ITP_spce
    expected_n_atoms = 7
    expected_n_residues = 1
    expected_n_segments = 1

    expected_n_bonds = 2
    expected_n_angles = 0
    expected_n_dihedrals = 0
    expected_n_impropers = 0

    @pytest.fixture
    def universe(self, filename):
        return mda.Universe(filename, HEAVY_H=True, EXTRA_ATOMS=True, HEAVY_SIX=True)

    @pytest.fixture()
    def top(self, filename):
        with self.parser(filename) as p:
            yield p.parse(HEAVY_H=True, EXTRA_ATOMS=True, HEAVY_SIX=True)
    
    def test_heavy_atom(self, universe):
        assert universe.atoms[5].mass > 40


class TestReadTop(BaseITP):
    """
    Test reading top files
    """

    ref_filename = GMX_TOP

    expected_n_atoms = 135
    expected_n_residues = 23
    expected_n_segments = 5

    expected_n_bonds = 130
    expected_n_angles = 182
    expected_n_dihedrals = 60
    expected_n_impropers = 58

    @pytest.fixture()
    def top(self, filename):
        with self.parser(filename) as p:
            yield p.parse(include_dir=GMX_DIR)

    @pytest.fixture()
    def universe(self, filename):
        return mda.Universe(filename, topology_format='ITP', include_dir=GMX_DIR)

    def test_output(self, filename):
        """Testing the call signature"""
        with self.parser(filename) as p:
            top = p.parse(include_dir=GMX_DIR)

    def test_creates_universe(self, filename):
        """Check that Universe works with this Parser"""
        u = mda.Universe(filename, topology_format='ITP', include_dir=GMX_DIR)

    def test_sequential(self, universe):
        resids = np.array(list(range(2, 12)) + list(range(13, 23)))
        assert_equal(universe.residues.resids[:20], resids)
        assert_equal(universe.residues.resindices, np.arange(self.expected_n_residues))
        assert_equal(universe.atoms.chargegroups[-1], 63)


class TestErrors:

    parser = mda.topology.ITPParser.ITPParser

    def test_no_include(self):
        with pytest.raises(IOError):
            with self.parser(GMX_TOP_BAD) as p:
                top = p.parse(include_dir=GMX_DIR)

    def test_missing_endif(self):
        with pytest.raises(IOError):
            with self.parser(ITP_no_endif) as p:
                top = p.parse(include_dir=GMX_DIR)


class TestRelativePath:
    def test_relstring(self, tmpdir):
        content = """ #include "../sub3/test2.itp"
        [ atoms ]
         1      H      1    SOL    HW1      1       0.41    1.00800
        """
        content2 = """[ atoms ]
         1      H      1    SOL    HW1      1       0.41    1.00800
        """
        p = tmpdir.mkdir("sub1").join("test.itp")
        p.write(content)
        p3 = tmpdir.mkdir("sub3").join("test2.itp")
        p3.write(content2)
        p2 = tmpdir.mkdir("sub2")
        p2.chdir()
        with p2.as_cwd() as pchange:
            u = mda.Universe(str("../sub1/test.itp"), format='ITP')

    def test_relpath(self, tmpdir):
        content = """
        [ atoms ]
         1      H      1    SOL    HW1      1       0.41    1.00800
        """
        p = tmpdir.mkdir("sub1").join("test.itp")
        p.write(content)
        p2 = tmpdir.mkdir("sub2")
        p2.chdir()
        with p2.as_cwd() as pchange:
            relpath = Path("../sub1/test.itp")
            u = mda.Universe(relpath, format='ITP')

    def test_relative_path(self, tmpdir):
        test_itp_content = '#include "../atoms.itp"'
        atoms_itp_content = """
        [ moleculetype ]
        UNK 3

        [ atoms ]
        1      H      1    SOL    HW1      1       0.41    1.00800
        """
        with tmpdir.as_cwd():
            with open("atoms.itp", "w") as f:
                f.write(atoms_itp_content)
            subdir = tmpdir.mkdir("subdir")
            with subdir.as_cwd():
                with open("test.itp", "w") as f:
                    f.write(test_itp_content)
                subsubdir = subdir.mkdir("subsubdir")
                with subsubdir.as_cwd():
                    u = mda.Universe("../test.itp")
                    assert len(u.atoms) == 1
