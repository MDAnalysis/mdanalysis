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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import division, absolute_import

from six.moves import range

import os
import itertools
import numpy as np
from numpy.testing import(
    assert_equal,
)

import MDAnalysis
import MDAnalysis as mda
import MDAnalysis.core.selection
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.core.selection import Parser
from MDAnalysis import SelectionError

from MDAnalysis.tests.datafiles import (
    PSF, DCD,
    PRMpbc, TRJpbc_bz2,
    PSF_NAMD, PDB_NAMD,
    GRO, NUCL, NUCLsel, TPR, XTC,
    TRZ_psf, TRZ,
    PDB_icodes,
)
from MDAnalysisTests import make_Universe

import pytest


class TestSelectionsCHARMM(object):
    @pytest.fixture(scope='class')
    def universe(self):
        """Set up the standard AdK system in implicit solvent.
        Geometry here is orthogonal
        """
        return MDAnalysis.Universe(PSF, DCD)

    @pytest.fixture()
    def universe_copy(self, universe):
        return MDAnalysis.Universe(PSF, DCD)

    def test_segid(self, universe):
        sel = universe.select_atoms('segid 4AKE')
        assert_equal(sel.n_atoms, 3341, "failed to select segment 4AKE")
        assert_equal(sorted(sel.indices),
                     sorted(universe.select_atoms('segid 4AKE').indices),
                     "selected segment 4AKE is not the same as auto-generated segment s4AKE")

    def test_protein(self, universe):
        sel = universe.select_atoms('protein')
        assert_equal(sel.n_atoms, 3341, "failed to select protein")
        assert_equal(sorted(sel.indices),
                     sorted(universe.select_atoms('segid 4AKE').indices),
                     "selected protein is not the same as auto-generated protein segment s4AKE")

    @pytest.mark.parametrize('resname', MDAnalysis.core.selection.ProteinSelection.prot_res)
    def test_protein_resnames(self, resname):
        u = make_Universe(('resnames',))
        # set half the residues' names to the resname we're testing
        myprot = u.residues[::2]
        # Windows note: the parametrized test input string objects
        # are actually of type 'numpy.str_' and coercion to str
        # proper is needed for unit test on Windows
        myprot.resnames = str(resname)
        # select protein
        sel = u.select_atoms('protein')
        # check that contents (atom indices) are identical afterwards
        assert_equal(myprot.atoms.ix, sel.ix)

    def test_backbone(self, universe):
        sel = universe.select_atoms('backbone')
        assert_equal(sel.n_atoms, 855)

    def test_resid_single(self, universe):
        sel = universe.select_atoms('resid 100')
        assert_equal(sel.n_atoms, 7)
        assert_equal(sel.residues.resnames, ['GLY'])

    def test_resid_range(self, universe):
        sel = universe.select_atoms('resid 100:105')
        assert_equal(sel.n_atoms, 89)
        assert_equal(sel.residues.resnames,
                     ['GLY', 'ILE', 'ASN', 'VAL', 'ASP', 'TYR'])

    def test_selgroup(self, universe):
        sel = universe.select_atoms('not resid 100')
        sel2 = universe.select_atoms('not group notr100', notr100=sel)
        assert_equal(sel2.n_atoms, 7)
        assert_equal(sel2.residues.resnames, ['GLY'])

    def test_fullselgroup(self, universe):
        sel1 = universe.select_atoms('resid 101')
        sel2 = universe.select_atoms('resid 100')
        sel3 = sel1.select_atoms('global group r100', r100=sel2)
        assert_equal(sel2.n_atoms, 7)
        assert_equal(sel2.residues.resnames, ['GLY'])

    # resnum selections are boring here because we haven't really a mechanism yet
    # to assign the canonical PDB resnums
    def test_resnum_single(self, universe):
        sel = universe.select_atoms('resnum 100')
        assert_equal(sel.n_atoms, 7)
        assert_equal(sel.residues.resids, [100])
        assert_equal(sel.residues.resnames, ['GLY'])

    def test_resnum_range(self, universe):
        sel = universe.select_atoms('resnum 100:105')
        assert_equal(sel.n_atoms, 89)
        assert_equal(sel.residues.resids, range(100, 106))
        assert_equal(sel.residues.resnames,
                     ['GLY', 'ILE', 'ASN', 'VAL', 'ASP', 'TYR'])

    def test_resname(self, universe):
        sel = universe.select_atoms('resname LEU')
        assert_equal(sel.n_atoms, 304,
                     "Failed to find all 'resname LEU' atoms.")
        assert_equal(sel.n_residues, 16,
                     "Failed to find all 'resname LEU' residues.")
        assert_equal(sorted(sel.indices),
                     sorted(universe.select_atoms('segid 4AKE and resname LEU').indices),
                     "selected 'resname LEU' atoms are not the same as auto-generated s4AKE.LEU")

    def test_name(self, universe):
        sel = universe.select_atoms('name CA')
        assert_equal(sel.n_atoms, 214)

    def test_atom(self, universe):
        sel = universe.select_atoms('atom 4AKE 100 CA')
        assert_equal(len(sel), 1)
        assert_equal(sel.resnames, ['GLY'])
        assert_equal(
            sel.positions,
            np.array([[20.38685226, -3.44224262, -5.92158318]],
                     dtype=np.float32))

    def test_atom_empty(self, universe):
        sel = universe.select_atoms('atom 4AKE 100 XX')  # Does not exist
        assert_equal(len(sel), 0)

    def test_type(self, universe):
        sel = universe.select_atoms("type 1")
        assert_equal(len(sel), 253)

    def test_and(self, universe):
        sel = universe.select_atoms('resname GLY and resid 100')
        assert_equal(len(sel), 7)

    def test_or(self, universe):
        sel = universe.select_atoms('resname LYS or resname ARG')
        assert_equal(sel.n_residues, 31)

    def test_not(self, universe):
        sel = universe.select_atoms('not backbone')
        assert_equal(len(sel), 2486)

    def test_around(self, universe):
        sel = universe.select_atoms('around 4.0 bynum 1943')
        assert_equal(len(sel), 32)

    def test_sphlayer(self, universe):
        sel = universe.select_atoms('sphlayer 4.0 6.0 bynum 1281')
        assert_equal(len(sel), 66)

    def test_sphzone(self, universe):
        sel = universe.select_atoms('sphzone 6.0 bynum 1281')
        assert_equal(len(sel), 86)

    def test_cylayer(self, universe):
        sel = universe.select_atoms('cylayer 4.0 6.0 10 -10 bynum 1281')
        assert_equal(len(sel), 88)

    def test_cyzone(self, universe):
        sel = universe.select_atoms('cyzone 6.0 10 -10 bynum 1281')
        assert_equal(len(sel), 166)

    def test_point(self, universe):
        ag = universe.select_atoms('point 5.0 5.0 5.0 3.5')

        d = distance_array(np.array([[5.0, 5.0, 5.0]], dtype=np.float32),
                           universe.atoms.positions,
                           box=universe.dimensions)

        idx = np.where(d < 3.5)[1]

        assert_equal(set(ag.indices), set(idx))

    def test_prop(self, universe):
        sel = universe.select_atoms('prop y <= 16')
        sel2 = universe.select_atoms('prop abs z < 8')
        assert_equal(len(sel), 3194)
        assert_equal(len(sel2), 2001)

    def test_bynum(self, universe):
        "Tests the bynum selection, also from AtomGroup instances (Issue 275)"
        sel = universe.select_atoms('bynum 5')
        assert_equal(sel[0].index, 4)
        sel = universe.select_atoms('bynum 1:10')
        assert_equal(len(sel), 10)
        assert_equal(sel[0].index, 0)
        assert_equal(sel[-1].index, 9)
        subsel = sel.select_atoms('bynum 5')
        assert_equal(subsel[0].index, 4)
        subsel = sel.select_atoms('bynum 2:5')
        assert_equal(len(subsel), 4)
        assert_equal(subsel[0].index, 1)
        assert_equal(subsel[-1].index, 4)

    def test_byres(self, universe):
        sel = universe.select_atoms('byres bynum 0:5')

        assert_equal(len(sel), len(universe.residues[0].atoms))

    def test_same_resname(self, universe):
        """Test the 'same ... as' construct (Issue 217)"""
        sel = universe.select_atoms("same resname as resid 10 or resid 11")
        assert_equal(len(sel), 331,
                     ("Found a wrong number of atoms with same resname as "
                      "resids 10 or 11"))
        target_resids = np.array([7, 8, 10, 11, 12, 14, 17, 25, 32, 37, 38,
                                  42, 46, 49, 55, 56, 66, 73, 80, 85, 93, 95,
                                  99, 100, 122, 127, 130, 144, 150, 176, 180,
                                  186, 188, 189, 194, 198, 203, 207, 214])
        assert_equal(sel.residues.resids, target_resids,
                     ("Found wrong residues with same resname as "
                      "resids 10 or 11"))

    def test_same_segment(self, universe_copy):
        """Test the 'same ... as' construct (Issue 217)"""

        SNew_A = universe_copy.add_Segment(segid='A')
        SNew_B = universe_copy.add_Segment(segid='B')
        SNew_C = universe_copy.add_Segment(segid='C')

        universe_copy.residues[:100].segments = SNew_A
        universe_copy.residues[100:150].segments = SNew_B
        universe_copy.residues[150:].segments = SNew_C

        target_resids = np.arange(100) + 1
        sel = universe_copy.select_atoms("same segment as resid 10")
        assert_equal(len(sel), 1520,
                     "Found a wrong number of atoms in the same segment of resid 10")
        assert_equal(sel.residues.resids,
                     target_resids,
                     "Found wrong residues in the same segment of resid 10")

        target_resids = np.arange(100, 150) + 1
        sel = universe_copy.select_atoms("same segment as resid 110")
        assert_equal(len(sel), 797,
                     "Found a wrong number of atoms in the same segment of resid 110")
        assert_equal(sel.residues.resids, target_resids,
                     "Found wrong residues in the same segment of resid 110")

        target_resids = np.arange(150, universe_copy.atoms.n_residues) + 1
        sel = universe_copy.select_atoms("same segment as resid 160")
        assert_equal(len(sel), 1024,
                     "Found a wrong number of atoms in the same segment of resid 160")
        assert_equal(sel.residues.resids, target_resids,
                     "Found wrong residues in the same segment of resid 160")

    def test_empty_same(self, universe):
        ag = universe.select_atoms('resname MET')

        # No GLY, so 'as resname GLY' is empty
        ag2 = ag.select_atoms('same mass as resname GLY')

        assert len(ag2) == 0

    def test_empty_selection(self, universe):
        """Test that empty selection can be processed (see Issue 12)"""
        # no Trp in AdK
        assert_equal(len(universe.select_atoms('resname TRP')), 0)

    def test_parenthesized_expression(self, universe):
        sel = universe.select_atoms(
            '( name CA or name CB ) and resname LEU')
        assert_equal(len(sel), 32)

    def test_no_space_around_parentheses(self, universe):
        """Test that no space is needed around parentheses (Issue 43)."""
        # note: will currently be ERROR because it throws a ParseError
        sel = universe.select_atoms('(name CA or name CB) and resname LEU')
        assert_equal(len(sel), 32)

    # TODO:
    # test for checking ordering and multiple comma-separated selections
    def test_concatenated_selection(self, universe):
        E151 = universe.select_atoms('segid 4AKE').select_atoms('resid 151')
        # note that this is not quite phi... HN should be C of prec. residue
        phi151 = E151.atoms.select_atoms('name HN', 'name N', 'name CA',
                                         'name CB')
        assert_equal(len(phi151), 4)
        assert_equal(phi151[0].name, 'HN',
                     "wrong ordering in selection, should be HN-N-CA-CB")

    def test_global(self, universe):
        """Test the `global` modifier keyword (Issue 268)"""
        ag = universe.select_atoms("resname LYS and name NZ")
        # Lys amines within 4 angstrom of the backbone.
        ag1 = universe.select_atoms(
            "resname LYS and name NZ and around 4 backbone")
        ag2 = ag.select_atoms("around 4 global backbone")
        assert_equal(ag2.indices, ag1.indices)


class TestSelectionsAMBER(object):
    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(PRMpbc, TRJpbc_bz2)

    def test_protein(self, universe):
        sel = universe.select_atoms('protein')
        assert_equal(sel.n_atoms, 22, "failed to select protein")

    def test_backbone(self, universe):
        sel = universe.select_atoms('backbone')
        assert_equal(sel.n_atoms, 7)

    def test_type(self, universe):
        sel = universe.select_atoms('type HC')
        assert_equal(len(sel), 6)
        assert_equal(sel.names, ['HH31', 'HH32', 'HH33', 'HB1', 'HB2', 'HB3'])


@pytest.mark.xfail(os.name == 'nt',
                   strict=True,
                   reason="Not supported on Windows yet.")
class TestSelectionsNAMD(object):
    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(PSF_NAMD, PDB_NAMD)

    def test_protein(self, universe):
        # must include non-standard residues
        sel = universe.select_atoms(
            'protein or resname HAO or resname ORT')
        assert_equal(sel.n_atoms, universe.atoms.n_atoms,
                     "failed to select peptide")
        assert_equal(sel.n_residues, 6,
                     "failed to select all peptide residues")

    def test_resid_single(self, universe):
        sel = universe.select_atoms('resid 12')
        assert_equal(sel.n_atoms, 26)
        assert_equal(sel.residues.resnames, ['HAO'])

    def test_type(self, universe):
        sel = universe.select_atoms('type H')
        assert_equal(len(sel), 5)
        # note 4th HH
        assert_equal(sel.names, ['HN', 'HN', 'HN', 'HH', 'HN'])


class TestSelectionsGRO(object):
    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(GRO)

    def test_protein(self, universe):
        sel = universe.select_atoms('protein')
        assert_equal(sel.n_atoms, 3341, "failed to select protein")

    def test_backbone(self, universe):
        sel = universe.select_atoms('backbone')
        assert_equal(sel.n_atoms, 855)

    def test_resid_single(self, universe):
        sel = universe.select_atoms('resid 100')
        assert_equal(sel.n_atoms, 7)
        assert_equal(sel.residues.resnames, ['GLY'])

    def test_same_coordinate(self, universe):
        """Test the 'same ... as' construct (Issue 217)"""
        sel = universe.select_atoms("same x as bynum 1 or bynum 10")
        assert_equal(len(sel), 12,
                     "Found a wrong number of atoms with same x as ids 1 or 10")
        target_ids = np.array([0, 8, 9, 224, 643, 3515,
                               11210, 14121, 18430, 25418, 35811, 43618])
        assert_equal(sel.indices, target_ids,
                     "Found wrong atoms with same x as ids 1 or 10")

    def test_cylayer(self, universe):
        """Cylinder layer selections with tricilinic periodicity (Issue 274)"""
        atgp = universe.select_atoms('name OW')
        sel = atgp.select_atoms('cylayer 10 20 20 -20 bynum 3554')
        assert_equal(len(sel), 1155)

    def test_cyzone(self, universe):
        """Cylinder zone selections with tricilinic periodicity (Issue 274)"""
        atgp = universe.select_atoms('name OW')
        sel = atgp.select_atoms('cyzone 20 20 -20 bynum 3554')
        assert_equal(len(sel), 1556)


class TestSelectionsTPR(object):
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(TPR,XTC)

    def test_same_fragment(self, universe):
        """Test the 'same ... as' construct (Issue 217)"""
        # This test comes here because it's a system with solvent,
        # and thus multiple fragments.
        sel = universe.select_atoms("same fragment as bynum 1")
        assert_equal(
            len(sel), 3341,
            "Found a wrong number of atoms on the same fragment as id 1")
        assert_equal(
            sel.indices, universe.atoms[0].fragment.indices,
            "Found a different set of atoms when using the 'same fragment as' construct vs. the .fragment prperty")

    def test_moltype(self, universe):
        sel = universe.select_atoms("moltype NA+")
        ref = np.array([47677, 47678, 47679, 47680], dtype=np.int32)
        assert_equal(sel.ids, ref)

    @pytest.mark.parametrize(
        'selection_string,reference',
        (('molnum 1', [3341, 3342, 3343, 3344]),
         ('molnum 2:4', [3345, 3346, 3347, 3348, 3349, 3350,
                         3351, 3352, 3353, 3354, 3355, 3356]),
         )
    )
    def test_molnum(self, universe, selection_string, reference):
        sel = universe.select_atoms(selection_string)
        assert_equal(sel.ids, np.array(reference, dtype=np.int32))


class TestSelectionsNucleicAcids(object):
    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(NUCL)

    def test_nucleic(self, universe):
        rna = universe.select_atoms("nucleic")
        assert_equal(rna.n_atoms, 739)
        assert_equal(rna.n_residues, 23)

    def test_nucleic_all(self, universe):
        u = mda.Universe(NUCLsel)

        sel = u.select_atoms('nucleic')

        assert len(sel) == 34

    def test_nucleicbackbone(self, universe):
        rna = universe.select_atoms("nucleicbackbone")
        assert_equal(rna.n_residues, 23)
        assert_equal(rna.n_atoms, rna.n_residues * 5 - 1)
        # -1 because this is a single strand of RNA and on P is missing at the 5' end

    # todo: need checks for other selection resnames such as DT DA DG DC DU

    def test_nucleicbase(self, universe):
        rna = universe.select_atoms("nucleicbase")
        assert_equal(rna.n_residues, 23)
        assert_equal(rna.n_atoms, 214)

    def test_nucleicsugar(self, universe):
        rna = universe.select_atoms("nucleicsugar")
        assert_equal(rna.n_residues, 23)
        assert_equal(rna.n_atoms, rna.n_residues * 5)


class BaseDistanceSelection(object):
    """Both KDTree and distmat selections on orthogonal system

    Selections to check:
     - Around
     - SphericalLayer
     - SphericalZone
     - Point

    Cylindrical methods don't use KDTree
    """

    methods = [('kdtree', False),
               ('kdtree', True),
               ('distmat', True),
               ('distmat', False)]

    @staticmethod
    def choosemeth(sel, meth, periodic):
        """hack in the desired apply method"""
        if meth == 'kdtree':
            sel.apply = sel._apply_KDTree
        elif meth == 'distmat':
            sel.apply = sel._apply_distmat

        if periodic:
            sel.periodic = True
        else:
            sel.periodic = False

        return sel

    @pytest.mark.parametrize('meth, periodic', methods)
    def test_around(self, u, meth, periodic):
        sel = Parser.parse('around 5.0 resid 1', u.atoms)
        sel = self.choosemeth(sel, meth, periodic)
        result = sel.apply(u.atoms)

        r1 = u.select_atoms('resid 1')
        cog = r1.center_of_geometry().reshape(1, 3)

        box = u.dimensions if periodic else None
        d = distance_array(u.atoms.positions, r1.positions,
                           box=box)
        ref = set(np.where(d < 5.0)[0])

        # Around doesn't include atoms from the reference group
        ref.difference_update(set(r1.indices))
        assert ref == set(result.indices)

    @pytest.mark.parametrize('meth, periodic', methods)
    def test_spherical_layer(self, u, meth, periodic):
        sel = Parser.parse('sphlayer 2.4 6.0 resid 1', u.atoms)
        sel = self.choosemeth(sel, meth, periodic)
        result = sel.apply(u.atoms)

        r1 = u.select_atoms('resid 1')
        box = u.dimensions if periodic else None
        cog = r1.center_of_geometry(pbc=periodic).reshape(1, 3)
        d = distance_array(u.atoms.positions, cog, box=box)
        ref = set(np.where((d > 2.4) & (d < 6.0))[0])

        assert ref == set(result.indices)

    @pytest.mark.parametrize('meth, periodic', methods)
    def test_spherical_zone(self, u, meth, periodic):
        sel = Parser.parse('sphzone 5.0 resid 1', u.atoms)
        sel = self.choosemeth(sel, meth, periodic)
        result = sel.apply(u.atoms)

        r1 = u.select_atoms('resid 1')
        box = u.dimensions if periodic else None
        cog = r1.center_of_geometry(pbc=periodic).reshape(1, 3)
        d = distance_array(u.atoms.positions, cog, box=box)
        ref = set(np.where(d < 5.0)[0])

        assert ref == set(result.indices)

    @pytest.mark.parametrize('meth, periodic', methods)
    def test_point(self, u, meth, periodic):
        sel = Parser.parse('point 5.0 5.0 5.0  3.0', u.atoms)
        sel = self.choosemeth(sel, meth, periodic)
        result = sel.apply(u.atoms)

        box = u.dimensions if periodic else None
        d = distance_array(np.array([[5.0, 5.0, 5.0]], dtype=np.float32),
                           u.atoms.positions,
                           box=box)
        ref = set(np.where(d < 3.0)[1])

        assert ref == set(result.indices)


class TestOrthogonalDistanceSelections(BaseDistanceSelection):
    @pytest.fixture()
    def u(self):
        return mda.Universe(TRZ_psf, TRZ)

    @pytest.mark.parametrize('meth, periodic', [
        ('distmat', True),
        ('distmat', False)
    ])
    def test_cyzone(self, u, meth, periodic):
        sel = Parser.parse('cyzone 5 4 -4 resid 2', u.atoms)
        sel.periodic = periodic
        result = sel.apply(u.atoms)

        other = u.select_atoms('resid 2')
        pos = other.center_of_geometry()

        vecs = u.atoms.positions - pos
        if periodic:
            box = u.dimensions[:3]
            vecs -= box * np.rint(vecs / box)

        mask = (vecs[:, 2] > -4) & (vecs[:, 2] < 4)

        radii = vecs[:, 0] ** 2 + vecs[:, 1] ** 2
        mask &= radii < 5 ** 2

        ref = set(u.atoms[mask].indices)

        assert ref == set(result.indices)


class TestTriclinicDistanceSelections(BaseDistanceSelection):
    @pytest.fixture()
    def u(self):
        return mda.Universe(GRO)


class TestTriclinicSelections(object):
    """Non-KDTree based selections

    This system has triclinic geometry so won't use KDTree based selections
    """

    @pytest.fixture()
    def u(self):
        return mda.Universe(GRO)

    def test_around(self, u):
        r1 = u.select_atoms('resid 1')

        ag = u.select_atoms('around 5.0 resid 1')

        d = distance_array(u.atoms.positions, r1.positions, box=u.dimensions)
        idx = set(np.where(d < 5.0)[0])

        # Around doesn't include atoms from the reference group
        idx.difference_update(set(r1.indices))
        assert idx == set(ag.indices)

    def test_sphlayer(self, u):
        r1 = u.select_atoms('resid 1')
        cog = r1.center_of_geometry().reshape(1, 3)

        ag = u.select_atoms('sphlayer 2.4 6.0 resid 1')

        d = distance_array(u.atoms.positions, cog, box=u.dimensions)
        idx = set(np.where((d > 2.4) & (d < 6.0))[0])

        assert idx == set(ag.indices)

    def test_sphzone(self, u):
        r1 = u.select_atoms('resid 1')
        cog = r1.center_of_geometry().reshape(1, 3)

        ag = u.select_atoms('sphzone 5.0 resid 1')

        d = distance_array(u.atoms.positions, cog, box=u.dimensions)
        idx = set(np.where(d < 5.0)[0])

        assert idx == set(ag.indices)

    def test_point_1(self, u):
        # The example selection
        ag = u.select_atoms('point 5.0 5.0 5.0 3.5')

        d = distance_array(np.array([[5.0, 5.0, 5.0]], dtype=np.float32),
                           u.atoms.positions,
                           box=u.dimensions)

        idx = np.where(d < 3.5)[1]

        assert_equal(set(ag.indices), set(idx))

    def test_point_2(self, u):
        ag1 = u.atoms[:10000]

        ag2 = ag1.select_atoms('point 5.0 5.0 5.0 3.5')

        d = distance_array(np.array([[5.0, 5.0, 5.0]], dtype=np.float32),
                           ag1.positions,
                           box=u.dimensions)

        idx = np.where(d < 3.5)[1]

        assert_equal(set(ag2.indices), set(idx))


def gen_sel_strings(prop, oper):
    """Generate all possible combinations of spaces in selection strings

    ie:
      'prop x < 1.5'
      'prop x< 1.5'
      'prop x <1.5'
      'prop x<1.5'

    """
    for x, y in itertools.product([' ', ''], [' ', '']):
        yield 'prop {prop}{spc1}{oper}{spc2}1.5'.format(
            prop=prop, spc1=x, oper=oper, spc2=y)


class TestPropSelection(object):
    plurals = {'mass': 'masses',
               'charge': 'charges'}
    op_funcs = {
        '<': np.less,
        '<=': np.less_equal,
        '>': np.greater,
        '>=': np.greater_equal,
        '==': np.equal,
        '!=': np.not_equal
    }
    opposites = {
        '==': '==', '!=': '!=',
        '>': '<=', '<=': '>',
        '<': '>=', '>=': '<',
    }

    @pytest.fixture(params=[slice(None, None), slice(None, 100)])
    def ag(self, request):
        u = make_Universe(('masses', 'charges'))
        u.atoms[::2].masses = 1.5
        u.atoms[::2].charges = 1.5
        return u.atoms[request.param]

    @pytest.mark.parametrize('prop, selstr', [
        (prop, sel)
        for prop in ['mass', 'charge']
        for sel in gen_sel_strings(prop, '<')
    ])
    def test_lt(self, prop, selstr, ag):
        sel = ag.select_atoms(selstr)
        assert_equal(set(sel.indices),
                     set(ag[getattr(ag, self.plurals[prop]) < 1.5].indices))

    @pytest.mark.parametrize('prop, selstr', [
        (prop, sel)
        for prop in ['mass', 'charge']
        for sel in gen_sel_strings(prop, '<=')
    ])
    def test_le(self, prop, selstr, ag):
        sel = ag.select_atoms(selstr)
        assert_equal(set(sel.indices),
                     set(ag[getattr(ag,
                                    self.plurals[prop]) <= 1.5].indices))

    @pytest.mark.parametrize('prop, selstr', [
        (prop, sel)
        for prop in ['mass', 'charge']
        for sel in gen_sel_strings(prop, '>')
    ])
    def test_gt(self, prop, selstr, ag):
        sel = ag.select_atoms(selstr)
        assert_equal(set(sel.indices),
                     set(ag[getattr(ag, self.plurals[prop]) > 1.5].indices))

    @pytest.mark.parametrize('prop, selstr', [
        (prop, sel)
        for prop in ['mass', 'charge']
        for sel in gen_sel_strings(prop, '>=')
    ])
    def test_ge(self, prop, selstr, ag):
        sel = ag.select_atoms(selstr)
        assert_equal(set(sel.indices),
                     set(ag[getattr(ag,
                                    self.plurals[prop]) >= 1.5].indices))

    @pytest.mark.parametrize('prop, selstr', [
        (prop, sel)
        for prop in ['mass', 'charge']
        for sel in gen_sel_strings(prop, '==')
    ])
    def test_eq(self, prop, selstr, ag):
        sel = ag.select_atoms(selstr)
        assert_equal(set(sel.indices),
                     set(ag[getattr(ag,
                                    self.plurals[prop]) == 1.5].indices))

    @pytest.mark.parametrize('prop, selstr', [
        (prop, sel)
        for prop in ['mass', 'charge']
        for sel in gen_sel_strings(prop, '!=')
    ])
    def test_ne(self, prop, selstr, ag):
        sel = ag.select_atoms(selstr)
        assert_equal(set(sel.indices),
                     set(ag[getattr(ag,
                                    self.plurals[prop]) != 1.5].indices))

    @pytest.mark.parametrize('prop, op', [
        (prop, op)
        for prop in ['mass', 'charge']
        for op in ('<', '>', '<=', '>=', '==', '!=')
    ])
    def test_flip(self, prop, ag, op):
        func = self.op_funcs[op]

        # reference group, doing things forwards
        ref = ag[func(getattr(ag, self.plurals[prop]), 1.5)]

        selstr = 'prop 1.5 {op} {prop}'.format(
            op=self.opposites[op], prop=prop)
        sel = ag.select_atoms(selstr)

        assert_equal(set(ref.indices), set(sel.indices))


class TestBondedSelection(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(PSF, DCD)

    def test_bonded_1(self, u):
        ag = u.select_atoms('type 2 and bonded name N')

        assert len(ag) == 3

    def test_nobonds_warns(self, u):
        u = make_Universe(('names',))

        # empty bond topology attr
        batt = mda.core.topologyattrs.Bonds([])
        u.add_TopologyAttr(batt)

        with pytest.warns(UserWarning):
            u.select_atoms('bonded name AAA')


class TestSelectionErrors(object):
    @staticmethod
    @pytest.fixture()
    def universe():
        return make_Universe(
            ('names', 'masses', 'resids', 'resnames', 'resnums'))

    @pytest.mark.parametrize('selstr', [
        'name and H',  # string selection
        'name )',
        'resid abcd',  # resid arg parsing selection
        'resnum 7a7',  # rangeselection arg parsing
        'resid 1-',
        'prop chicken == tasty',
        'prop chicken <= 7.4',
        'prop mass ^^ 12.0',
        'same this as resid 1',  # same selection
        'same resid resname mass 5.0',  # same / expect
        'name H and',  # check all tokens used
        'naem H',  # unkonwn (misplet) opertaor
        'resid and name C',  # rangesel not finding vals
        'resnum ',
        'bynum or protein',
        'prop mass < 4.0 hello',  # unused token
        'prop mass > 10. and group this',  # missing group
        'prop mass > 10. and fullgroup this',  # missing fullgroup
    ])
    def test_selection_fail(self, selstr, universe):
        with pytest.raises(SelectionError):
            universe.select_atoms(selstr)


def test_segid_and_resid():
    u = make_Universe(('segids', 'resids'))

    ag = u.select_atoms('segid SegB and resid 1-100')

    ref = ag.select_atoms('segid SegB').select_atoms('resid 1-100')

    assert_equal(ag.indices, ref.indices)


class TestImplicitOr(object):
    @staticmethod
    @pytest.fixture()
    def universe():
        return make_Universe(
            ('names', 'types', 'resids', 'resnums', 'resnames', 'segids'))

    def _check_sels(self, ref, sel, universe):
        ref = universe.select_atoms(ref)
        sel = universe.select_atoms(sel)

        assert_equal(ref.indices, sel.indices)

    @pytest.mark.parametrize('ref, sel', [
        ('name NameABA or name NameACA or name NameADA',
         'name NameABA NameACA NameADA'),
        ('type TypeE or type TypeD or type TypeB', 'type TypeE TypeD TypeB'),
        ('resname RsC or resname RsY', 'resname RsC RsY'),
        ('name NameAB* or name NameACC', 'name NameAB* NameACC'),
        ('segid SegA or segid SegC', 'segid SegA SegC'),
        ('(name NameABC or name NameABB) and (resname RsD or resname RsF)',
         'name NameABC NameABB and resname RsD RsF'),
    ])
    def test_string_selections(self, ref, sel, universe):
        self._check_sels(ref, sel, universe)

    @pytest.mark.parametrize("seltype", ['resid', 'resnum', 'bynum'])
    @pytest.mark.parametrize('ref, sel', [
        ('{typ} 1 or {typ} 2', '{typ} 1 2'),
        ('{typ} 1:10 or {typ} 22', '{typ} 1:10 22'),
        ('{typ} 1:10 or {typ} 20:30', '{typ} 1:10 20:30'),
        ('{typ} 1-5 or {typ} 7', '{typ} 1-5 7'),
        ('{typ} 1-5 or {typ} 7:10 or {typ} 12', '{typ} 1-5 7:10 12'),
        ('{typ} 1 or {typ} 3 or {typ} 5:10', '{typ} 1 3 5:10'),
    ])
    def test_range_selections(self, seltype, ref, sel, universe):
        self._check_sels(ref.format(typ=seltype), sel.format(typ=seltype),
                         universe)


class TestICodeSelection(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(PDB_icodes)

    def test_select_icode(self, u):
        ag = u.select_atoms('resid 163A')

        assert len(ag) == 7
        assert_equal(ag.ids, np.arange(7) + 1230)

    def test_select_resid_implicit_icode(self, u):
        ag = u.select_atoms('resid 163')

        assert len(ag) == 6
        assert_equal(ag.ids, np.arange(6) + 1224)

    def test_select_icode_range_1(self, u):
        # testing range within a single resid integer value
        u = u
        ag = u.select_atoms('resid 163B-163D')

        # do it manually without selection language...
        ref = u.residues[u.residues.resids == 163]
        ref = ref[(ref.icodes >= 'B') & (ref.icodes <= 'D')]
        ref = ref.atoms

        assert_equal(ag.ids, ref.ids)

        assert len(ag) == 19
        assert_equal(ag.ids, np.arange(19) + 1237)

    def test_select_icode_range_2(self, u):
        u = u

        ag = u.select_atoms('resid 163B-165')

        resids = u.residues.resids

        start = u.residues[resids == 163]
        start = start[start.icodes >= 'B']

        mid = u.residues[resids == 164]

        end = u.residues[resids == 165]
        end = end[end.icodes == '']

        ref = start.atoms + mid.atoms + end.atoms

        assert_equal(ag.ids, ref.ids)

    def test_select_icode_range_3(self, u):
        # same as #2 but with no "middle" icodes
        u = u

        ag = u.select_atoms('resid 163B-164')

        resids = u.residues.resids

        start = u.residues[resids == 163]
        start = start[start.icodes >= 'B']

        end = u.residues[resids == 164]
        end = end[end.icodes == '']

        ref = start.atoms + end.atoms

        assert_equal(ag.ids, ref.ids)

    def test_select_icode_range_4(self, u):
        u = u

        ag = u.select_atoms('resid 160-163G')

        resids = u.residues.resids

        start = u.residues[resids == 160]
        start = start[start.icodes >= '']

        mid = u.residues[(resids == 161) | (resids == 162)]

        end = u.residues[resids == 163]
        end = end[end.icodes <= 'G']

        ref = start.atoms + mid.atoms + end.atoms

        assert_equal(ag.ids, ref.ids)

    def test_select_icode_range_5(self, u):
        # same as #4 but with no "middle" icodes in range
        u = u

        ag = u.select_atoms('resid 162-163G')

        resids = u.residues.resids

        start = u.residues[resids == 162]
        start = start[start.icodes >= '']

        end = u.residues[resids == 163]
        end = end[end.icodes <= 'G']

        ref = start.atoms + end.atoms

        assert_equal(ag.ids, ref.ids)

    def test_missing_icodes_VE(self, u):
        # trying a selection with icodes in a Universe without raises VA
        u = make_Universe(('resids',))
        with pytest.raises(ValueError):
            u.select_atoms('resid 10A')

    def test_missing_icodes_range_VE(self, u):
        u = make_Universe(('resids',))
        with pytest.raises(ValueError):
            u.select_atoms('resid 10A-12')


def test_arbitrary_atom_group_raises_error():
    u = make_Universe(trajectory=True)
    with pytest.raises(TypeError):
        u.select_atoms('around 2.0 group this', this=u.atoms[0])
