# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from six.moves import range

import numpy as np
from numpy.testing import(
    TestCase,
    dec,
    assert_equal,
    assert_array_almost_equal,
    assert_,
    assert_array_equal,
    assert_warns,
    assert_raises,
)
from nose.plugins.attrib import attr
import warnings

import MDAnalysis
import MDAnalysis as mda
import MDAnalysis.core.Selection
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.core.topologyobjects import TopologyGroup
from MDAnalysis.core.Selection import Parser
from MDAnalysis import SelectionError

from MDAnalysis.tests.datafiles import (
    PSF, DCD,
    PRMpbc, TRJpbc_bz2,
    PSF_NAMD, PDB_NAMD,
    GRO, NUCL, NUCLsel, TPR, XTC,
    TRZ_psf, TRZ,
    PDB_full,
)
from MDAnalysisTests.plugins.knownfailure import knownfailure
from MDAnalysisTests import parser_not_found


class TestSelectionsCHARMM(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        """Set up the standard AdK system in implicit solvent.

        Geometry here is orthogonal
        """
        self.universe = MDAnalysis.Universe(PSF, DCD)

    def tearDown(self):
        self.universe.trajectory.close()
        del self.universe

    def test_segid(self):
        sel = self.universe.select_atoms('segid 4AKE')
        assert_equal(sel.n_atoms, 3341, "failed to select segment 4AKE")
        assert_array_equal(sel.indices, self.universe.s4AKE.indices,
                     ("selected segment 4AKE is not the same as "
                      "auto-generated segment s4AKE"))

    def test_protein(self):
        sel = self.universe.select_atoms('protein')
        assert_equal(sel.n_atoms, 3341, "failed to select protein")
        assert_array_equal(sel.indices, self.universe.s4AKE.indices,
                     ("selected protein is not the same as auto-generated"
                      " protein segment s4AKE"))

    def test_backbone(self):
        sel = self.universe.select_atoms('backbone')
        assert_equal(sel.n_atoms, 855)

    def test_resid_single(self):
        sel = self.universe.select_atoms('resid 100')
        assert_equal(sel.n_atoms, 7)
        assert_equal(sel.residues.resnames, ['GLY'])

    def test_resid_range(self):
        sel = self.universe.select_atoms('resid 100:105')
        assert_equal(sel.n_atoms, 89)
        assert_equal(sel.residues.resnames, ['GLY', 'ILE', 'ASN', 'VAL', 'ASP', 'TYR'])

    def test_selgroup(self):
        sel = self.universe.select_atoms('not resid 100')
        sel2 = self.universe.select_atoms('not group notr100', notr100=sel)
        assert_equal(sel2.n_atoms, 7)
        assert_equal(sel2.residues.resnames, ['GLY'])

    def test_fullselgroup(self):
        sel1 = self.universe.select_atoms('resid 101')
        sel2 = self.universe.select_atoms('resid 100')
        sel3 = sel1.select_atoms('fullgroup r100', r100=sel2)
        assert_equal(sel2.n_atoms, 7)
        assert_equal(sel2.residues.resnames, ['GLY'])

    # resnum selections are boring here because we haven't really a mechanism
    # yet to assign the canonical PDB resnums
    def test_resnum_single(self):
        sel = self.universe.select_atoms('resnum 100')
        assert_equal(sel.n_atoms, 7)
        assert_equal(sel.residues.resids, [100])
        assert_equal(sel.residues.resnames, ['GLY'])

    def test_resnum_range(self):
        sel = self.universe.select_atoms('resnum 100:105')
        assert_equal(sel.n_atoms, 89)
        assert_equal(sel.residues.resids, range(100, 106))
        assert_equal(sel.residues.resnames,
                     ['GLY', 'ILE', 'ASN', 'VAL', 'ASP', 'TYR'])

    def test_resname(self):
        sel = self.universe.select_atoms('resname LEU')
        assert_equal(sel.n_atoms, 304,
                     "Failed to find all 'resname LEU' atoms.")
        assert_equal(sel.n_residues, 16,
                     "Failed to find all 'resname LEU' residues.")
        assert_array_equal(sel.indices, self.universe.s4AKE.LEU.indices,
                           ("selected 'resname LEU' atoms are not the same as "
                            "auto-generated s4AKE.LEU"))

    def test_name(self):
        sel = self.universe.select_atoms('name CA')
        assert_equal(sel.n_atoms, 214)

    def test_atom(self):
        sel = self.universe.select_atoms('atom 4AKE 100 CA')
        assert_equal(len(sel), 1)
        assert_equal(sel.resnames, ['GLY'])
        assert_array_almost_equal(
            sel.positions,
            np.array([[20.38685226, -3.44224262, -5.92158318]],
                     dtype=np.float32))

    def test_atom_empty(self):
        sel = self.universe.select_atoms('atom 4AKE 100 XX')  # Does not exist
        assert_equal(len(sel), 0)

    def test_type(self):
        sel = self.universe.select_atoms("type 1")
        assert_equal(len(sel), 253)

    def test_and(self):
        sel = self.universe.select_atoms('resname GLY and resid 100')
        assert_equal(len(sel), 7)

    def test_or(self):
        sel = self.universe.select_atoms('resname LYS or resname ARG')
        assert_equal(sel.n_residues, 31)

    def test_not(self):
        sel = self.universe.select_atoms('not backbone')
        assert_equal(len(sel), 2486)

    def test_around(self):
        sel = self.universe.select_atoms('around 4.0 bynum 1943')
        assert_equal(len(sel), 32)

    def test_sphlayer(self):
        sel = self.universe.select_atoms('sphlayer 4.0 6.0 bynum 1281')
        assert_equal(len(sel), 66)

    def test_sphzone(self):
        sel = self.universe.select_atoms('sphzone 6.0 bynum 1281')
        assert_equal(len(sel), 86)

    def test_cylayer(self):
        sel = self.universe.select_atoms('cylayer 4.0 6.0 10 -10 bynum 1281')
        assert_equal(len(sel), 88)

    def test_cyzone(self):
        sel = self.universe.select_atoms('cyzone 6.0 10 -10 bynum 1281')
        assert_equal(len(sel), 166)

    def test_point(self):
        ag = self.universe.select_atoms('point 5.0 5.0 5.0 3.5')

        d = distance_array(np.array([[5.0, 5.0, 5.0]], dtype=np.float32),
                           self.universe.atoms.positions,
                           box=self.universe.dimensions)

        idx = np.where(d < 3.5)[1]

        assert_equal(set(ag.indices), set(idx))

    def test_prop(self):
        sel = self.universe.select_atoms('prop y <= 16')
        sel2 = self.universe.select_atoms('prop abs z < 8')
        assert_equal(len(sel), 3194)
        assert_equal(len(sel2), 2001)

    def test_bynum(self):
        "Tests the bynum selection, also from AtomGroup instances (Issue 275)"
        sel = self.universe.select_atoms('bynum 5')
        assert_equal(sel[0].index, 4)
        sel = self.universe.select_atoms('bynum 1:10')
        assert_equal(len(sel), 10)
        assert_equal(sel[0].index, 0)
        assert_equal(sel[-1].index, 9)
        subsel = sel.select_atoms('bynum 5')
        assert_equal(subsel[0].index, 4)
        subsel = sel.select_atoms('bynum 2:5')
        assert_equal(len(subsel), 4)
        assert_equal(subsel[0].index, 1)
        assert_equal(subsel[-1].index, 4)

    # TODO:
    # and also for selection keywords such as 'nucleic'
    def test_byres(self):
        sel = self.universe.select_atoms('byres bynum 0:5')

        assert_equal(len(sel), len(self.universe.residues[0]))

    def test_same_resname(self):
        """Test the 'same ... as' construct (Issue 217)"""
        sel = self.universe.select_atoms("same resname as resid 10 or resid 11")
        assert_equal(len(sel), 331,
                     ("Found a wrong number of atoms with same resname as "
                      "resids 10 or 11"))
        target_resids = np.array([ 7, 8, 10, 11, 12, 14, 17, 25, 32, 37, 38,
                                   42, 46, 49, 55, 56, 66, 73, 80, 85, 93, 95,
                                   99, 100, 122, 127, 130, 144, 150, 176, 180,
                                   186, 188, 189, 194, 198, 203, 207, 214])
        assert_array_equal(sel.residues.resids, target_resids,
                           ("Found wrong residues with same resname as "
                            "resids 10 or 11"))

    def test_same_segment(self):
        """Test the 'same ... as' construct (Issue 217)"""
        self.universe.residues[:100].set_segids("A")  # make up some segments
        self.universe.residues[100:150].set_segids("B")
        self.universe.residues[150:].set_segids("C")

        target_resids = np.arange(100)+1
        sel = self.universe.select_atoms("same segment as resid 10")
        assert_equal(len(sel), 1520,
                     "Found a wrong number of atoms in the same "
                     "segment of resid 10")
        assert_array_equal(sel.residues.resids, target_resids,
                           ("Found wrong residues in the same segment"
                            " of resid 10"))

        target_resids = np.arange(100,150)+1
        sel = self.universe.select_atoms("same segment as resid 110")
        assert_equal(len(sel), 797,
                     ("Found a wrong number of atoms in the same"
                      " segment of resid 110"))
        assert_array_equal(sel.residues.resids, target_resids,
                           ("Found wrong residues in the same segment"
                            " of resid 110"))

        target_resids = np.arange(150,self.universe.atoms.n_residues)+1
        sel = self.universe.select_atoms("same segment as resid 160")
        assert_equal(len(sel), 1024,
                     ("Found a wrong number of atoms in the same"
                      " segment of resid 160"))
        assert_array_equal(sel.residues.resids, target_resids,
                           ("Found wrong residues in the same segment of"
                            " resid 160"))

        #cleanup
        self.universe.residues.set_segids("4AKE")

    def test_empty_same(self):
        ag = self.universe.select_atoms('resname MET')

        # No GLY, so 'as resname GLY' is empty
        ag2 = ag.select_atoms('same mass as resname GLY')

        assert_(len(ag2) == 0)

    def test_empty_selection(self):
        """Test that empty selection can be processed (see Issue 12)"""
        # no Trp in AdK
        assert_equal(len(self.universe.select_atoms('resname TRP')), 0)

    def test_parenthesized_expression(self):
        sel = self.universe.select_atoms(
            '( name CA or name CB ) and resname LEU')
        assert_equal(len(sel), 32)

    def test_no_space_around_parentheses(self):
        """Test that no space is needed around parentheses (Issue 43)."""
        # note: will currently be ERROR because it throws a ParseError
        sel = self.universe.select_atoms('(name CA or name CB) and resname LEU')
        assert_equal(len(sel), 32)

    # TODO:
    # test for checking ordering and multiple comma-separated selections
    def test_concatenated_selection(self):
        E151 = self.universe.s4AKE.r151
        # note that this is not quite phi... HN should be C of prec. residue
        phi151 = E151.select_atoms('name HN', 'name N', 'name CA', 'name CB')
        assert_equal(len(phi151), 4)
        assert_equal(phi151[0].name, 'HN',
                     "wrong ordering in selection, should be HN-N-CA-CB")

    def test_global(self):
        """Test the `global` modifier keyword (Issue 268)"""
        ag = self.universe.select_atoms("resname LYS and name NZ")
        # Lys amines within 4 angstrom of the backbone.
        ag1 = self.universe.select_atoms(
            "resname LYS and name NZ and around 4 backbone")
        ag2 = ag.select_atoms("around 4 global backbone")
        assert_array_equal(ag2.indices, ag1.indices)


class TestSelectionsAMBER(TestCase):
    def setUp(self):
        """Set up AMBER system"""
        self.universe = MDAnalysis.Universe(PRMpbc, TRJpbc_bz2)

    def tearDown(self):
        self.universe.trajectory.close()
        del self.universe

    def test_protein(self):
        sel = self.universe.select_atoms('protein')
        assert_equal(sel.n_atoms, 22, "failed to select protein")

    def test_backbone(self):
        sel = self.universe.select_atoms('backbone')
        assert_equal(sel.n_atoms, 7)

    def test_resid_single(self):
        sel = self.universe.select_atoms('resid 3')
        assert_equal(sel.n_atoms, 6)
        assert_equal(sel.residues.resnames, ['NME'])

    def test_type(self):
        sel = self.universe.select_atoms('type 1')
        assert_equal(len(sel), 6)
        assert_equal(sel.names, ['HH31', 'HH32', 'HH33', 'HB1', 'HB2', 'HB3'])


class TestSelectionsNAMD(TestCase):
    def setUp(self):
        """Set up NAMD system"""
        self.universe = MDAnalysis.Universe(PSF_NAMD, PDB_NAMD)

    def tearDown(self):
        self.universe.trajectory.close()
        del self.universe

    def test_protein(self):
        # must include non-standard residues
        sel = self.universe.select_atoms(
            'protein or resname HAO or resname ORT')
        assert_equal(sel.n_atoms, self.universe.atoms.n_atoms,
                     "failed to select peptide")
        assert_equal(sel.n_residues, 6,
                     "failed to select all peptide residues")

    def test_resid_single(self):
        sel = self.universe.select_atoms('resid 12')
        assert_equal(sel.n_atoms, 26)
        assert_equal(sel.residues.resnames, ['HAO'])

    def test_type(self):
        sel = self.universe.select_atoms('type H')
        assert_equal(len(sel), 5)
        # note 4th HH
        assert_array_equal(sel.names, ['HN', 'HN', 'HN', 'HH', 'HN'])


class TestSelectionsGRO(TestCase):
    def setUp(self):
        """Set up GRO system (implicit types, charges, masses, ...)"""
        self.universe = MDAnalysis.Universe(GRO)

    @dec.slow
    def test_protein(self):
        sel = self.universe.select_atoms('protein')
        assert_equal(sel.n_atoms, 3341, "failed to select protein")

    @dec.slow
    def test_backbone(self):
        sel = self.universe.select_atoms('backbone')
        assert_equal(sel.n_atoms, 855)

    @dec.slow
    def test_type(self):
        sel = self.universe.select_atoms('type H')
        assert_equal(len(sel), 23853)
        sel = self.universe.select_atoms('type S')
        assert_equal(len(sel), 7)
        assert_equal(
                sel.residues.resnames,
                self.universe.select_atoms(
                    "resname CYS or resname MET").residues.resnames)

    @dec.slow
    def test_resid_single(self):
        sel = self.universe.select_atoms('resid 100')
        assert_equal(sel.n_atoms, 7)
        assert_equal(sel.residues.resnames, ['GLY'])

    @dec.slow
    def test_atom(self):
        sel = self.universe.select_atoms('atom SYSTEM 100 CA')
        assert_equal(len(sel), 1)
        assert_equal(sel.resnames, ['GLY'])

    @dec.slow
    def test_same_coordinate(self):
        """Test the 'same ... as' construct (Issue 217)"""
        sel = self.universe.select_atoms("same x as bynum 1 or bynum 10")
        assert_equal(len(sel), 12,
                     "Found a wrong number of atoms with same x as ids 1 or 10")
        target_ids = np.array([ 0, 8, 9, 224, 643, 3515,
                                11210, 14121, 18430, 25418, 35811, 43618])
        assert_array_equal(sel.indices, target_ids,
                           "Found wrong atoms with same x as ids 1 or 10")

    def test_cylayer(self):
        """Cylinder layer selections with tricilinic periodicity (Issue 274)"""
        atgp = self.universe.select_atoms('name OW')
        sel = atgp.select_atoms('cylayer 10 20 20 -20 bynum 3554')
        assert_equal(len(sel), 1155)

    def test_cyzone(self):
        """Cylinder zone selections with tricilinic periodicity (Issue 274)"""
        atgp = self.universe.select_atoms('name OW')
        sel = atgp.select_atoms('cyzone 20 20 -20 bynum 3554')
        assert_equal(len(sel), 1556)


class TestSelectionsXTC(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(TPR,XTC)

    def test_same_fragment(self):
        """Test the 'same ... as' construct (Issue 217)"""
        # This test comes here because it's a system with solvent,
        # and thus multiple fragments.
        try:
            sel = self.universe.select_atoms("same fragment as bynum 1")
            assert_equal(len(sel), 3341,
                         ("Found a wrong number of atoms on the same"
                          " fragment as id 1"))
            assert_array_equal(
                sel.indices, self.universe.atoms[0].fragment.indices,
                ("Found a different set of atoms when using the"
                 " 'same fragment as' construct vs. the .fragment "
                 "property"))
        except MDAnalysis.NoDataError:
            assert_equal(True, False)


class TestSelectionsNucleicAcids(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(NUCL)

    def test_nucleic(self):
        rna = self.universe.select_atoms("nucleic")
        assert_equal(rna.n_atoms, 739)
        assert_equal(rna.n_residues, 23)

    def test_nucleic_all(self):
        u = mda.Universe(NUCLsel)

        sel = u.select_atoms('nucleic')

        assert_(len(sel) == 34)

    def test_nucleicbackbone(self):
        rna = self.universe.select_atoms("nucleicbackbone")
        assert_equal(rna.n_residues, 23)
        assert_equal(rna.n_atoms, rna.n_residues * 5 - 1)
        # -1 because this is a single strand of RNA and on P is missing at the 5' end

    # todo: need checks for other selection resnames such as DT DA DG DC DU

    def test_nucleicbase(self):
        rna = self.universe.select_atoms("nucleicbase")
        assert_equal(rna.n_residues, 23)
        assert_equal(rna.n_atoms, 214)

    def test_nucleicsugar(self):
        rna = self.universe.select_atoms("nucleicsugar")
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
               ('distmat', True),
               ('distmat', False)]

    def choosemeth(self, sel, meth, periodic):
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

    def _check_around(self, meth, periodic):
        sel = Parser.parse('around 5.0 resid 1', self.u.atoms)
        sel = self.choosemeth(sel, meth, periodic)
        result = sel.apply(self.u.atoms)

        r1 = self.u.select_atoms('resid 1')
        cog = r1.center_of_geometry().reshape(1, 3)

        box = self.u.dimensions if periodic else None
        d = distance_array(self.u.atoms.positions, r1.positions,
                           box=box)
        ref = set(np.where(d < 5.0)[0])

        # Around doesn't include atoms from the reference group
        ref.difference_update(set(r1.indices))
        assert_(ref == set(result.indices))

    def test_around(self):
        for meth, periodic in self.methods:
            yield self._check_around, meth, periodic

    def _check_spherical_layer(self, meth, periodic):
        sel = Parser.parse('sphlayer 2.4 6.0 resid 1' , self.u.atoms)
        sel = self.choosemeth(sel, meth, periodic)
        result = sel.apply(self.u.atoms)

        r1 = self.u.select_atoms('resid 1')
        cog = r1.center_of_geometry().reshape(1, 3)

        box = self.u.dimensions if periodic else None
        d = distance_array(self.u.atoms.positions, cog, box=box)
        ref = set(np.where((d > 2.4) & (d < 6.0))[0])

        assert_(ref == set(result.indices))

    def test_spherical_layer(self):
        for meth, periodic in self.methods:
            yield self._check_spherical_layer, meth, periodic

    def _check_spherical_zone(self, meth, periodic):
        sel = Parser.parse('sphzone 5.0 resid 1', self.u.atoms)
        sel = self.choosemeth(sel, meth, periodic)
        result = sel.apply(self.u.atoms)

        r1 = self.u.select_atoms('resid 1')
        cog = r1.center_of_geometry().reshape(1, 3)

        box = self.u.dimensions if periodic else None
        d = distance_array(self.u.atoms.positions, cog, box=box)
        ref = set(np.where(d < 5.0)[0])

        assert_(ref == set(result.indices))

    def test_spherical_zone(self):
        for meth, periodic in self.methods:
            yield self._check_spherical_zone, meth, periodic

    def _check_point(self, meth, periodic):
        sel = Parser.parse('point 5.0 5.0 5.0  3.0', self.u.atoms)
        sel = self.choosemeth(sel, meth, periodic)
        result = sel.apply(self.u.atoms)

        box = self.u.dimensions if periodic else None
        d = distance_array(np.array([[5.0, 5.0, 5.0]], dtype=np.float32),
                           self.u.atoms.positions,
                           box=box)
        ref = set(np.where(d < 3.0)[1])

        assert_(ref == set(result.indices))

    def test_point(self):
        for meth, periodic in self.methods:
            yield self._check_point, meth, periodic


class TestOrthogonalDistanceSelections(BaseDistanceSelection):
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def setUp(self):
        self.u = mda.Universe(TRZ_psf, TRZ)

    def tearDown(self):
        del self.u

    def _check_cyzone(self, meth, periodic):
        sel = Parser.parse('cyzone 5 4 -4 resid 2', self.u.atoms)
        sel.periodic = periodic
        result = sel.apply(self.u.atoms)

        other = self.u.select_atoms('resid 2')
        pos = other.center_of_geometry()

        vecs = self.u.atoms.positions - pos
        if periodic:
            box = self.u.dimensions[:3]
            vecs -= box * np.rint(vecs / box)

        mask = (vecs[:,2] > -4) & (vecs[:,2] < 4)

        radii = vecs[:,0] ** 2 + vecs[:, 1] ** 2
        mask &= radii < 5**2

        ref = set(self.u.atoms[mask].indices)

        assert_(ref == set(result.indices))

    def test_cyzone(self):
        for meth, periodic in self.methods[1:]:
            yield self._check_cyzone, meth, periodic


class TestTriclinicDistanceSelections(BaseDistanceSelection):
    def setUp(self):
        self.u = mda.Universe(GRO)

    def tearDown(self):
        del self.u

class TestTriclinicSelections(object):
    """Non-KDTree based selections

    This system has triclinic geometry so won't use KDTree based selections
    """
    def setUp(self):
        self.u = mda.Universe(GRO)
        self.box = self.u.dimensions

    def tearDown(self):
        del self.u

    def test_around(self):
        r1 = self.u.select_atoms('resid 1')

        ag = self.u.select_atoms('around 5.0 resid 1')

        d = distance_array(self.u.atoms.positions, r1.positions, box=self.box)
        idx = set(np.where(d < 5.0)[0])

        # Around doesn't include atoms from the reference group
        idx.difference_update(set(r1.indices))
        assert_(idx == set(ag.indices))

    def test_sphlayer(self):
        r1 = self.u.select_atoms('resid 1')
        cog = r1.center_of_geometry().reshape(1, 3)

        ag = self.u.select_atoms('sphlayer 2.4 6.0 resid 1')

        d = distance_array(self.u.atoms.positions, cog, box=self.box)
        idx = set(np.where((d > 2.4) & (d < 6.0))[0])

        assert_(idx == set(ag.indices))

    def test_sphzone(self):
        r1 = self.u.select_atoms('resid 1')
        cog = r1.center_of_geometry().reshape(1, 3)

        ag = self.u.select_atoms('sphzone 5.0 resid 1')

        d = distance_array(self.u.atoms.positions, cog, box=self.box)
        idx = set(np.where(d < 5.0)[0])

        assert_(idx == set(ag.indices))

    def test_point_1(self):
        # The example selection
        ag = self.u.select_atoms('point 5.0 5.0 5.0 3.5')

        d = distance_array(np.array([[5.0, 5.0, 5.0]], dtype=np.float32),
                           self.u.atoms.positions,
                           box=self.box)

        idx = np.where(d < 3.5)[1]

        assert_equal(set(ag.indices), set(idx))

    def test_point_2(self):
        ag1 = self.u.atoms[:10000]

        ag2 = ag1.select_atoms('point 5.0 5.0 5.0 3.5')

        d = distance_array(np.array([[5.0, 5.0, 5.0]], dtype=np.float32),
                           ag1.positions,
                           box=self.box)

        idx = np.where(d < 3.5)[1]

        assert_equal(set(ag2.indices), set(idx))


class TestPropSelection(object):
    plurals = {'mass': 'masses',
               'charge': 'charges'}

    def _check_lt(self, prop, ag):
        setattr(ag[::2], self.plurals[prop], 500.0)

        sel = ag.select_atoms('prop {0} < 500.0'.format(prop))

        assert_equal(set(sel.indices),
                     set(ag[getattr(ag, self.plurals[prop]) < 500.0].indices))

    def _check_le(self, prop, ag):
        setattr(ag[::2], self.plurals[prop], 500.0)

        sel = ag.select_atoms('prop {0} <= 500.0'.format(prop))

        assert_equal(set(sel.indices),
                     set(ag[getattr(ag, self.plurals[prop]) <= 500.0].indices))

    def _check_gt(self, prop, ag):
        setattr(ag[::2], self.plurals[prop], 500.0)

        sel = ag.select_atoms('prop {0} > 500.0'.format(prop))

        assert_equal(set(sel.indices),
                     set(ag[getattr(ag, self.plurals[prop]) > 500.0].indices))

    def _check_ge(self, prop, ag):
        setattr(ag[::2], self.plurals[prop], 500.0)

        sel = ag.select_atoms('prop {0} >= 500.0'.format(prop))

        assert_equal(set(sel.indices),
                     set(ag[getattr(ag, self.plurals[prop]) >= 500.0].indices))

    def _check_eq(self, prop, ag):
        setattr(ag[::2], self.plurals[prop], 500.0)

        sel = ag.select_atoms('prop {0} == 500.0'.format(prop))

        assert_equal(set(sel.indices),
                     set(ag[getattr(ag, self.plurals[prop]) == 500.0].indices))

    def _check_ne(self, prop, ag):
        setattr(ag[::2], self.plurals[prop], 500.0)

        sel = ag.select_atoms('prop {0} != 500.0'.format(prop))

        assert_equal(set(sel.indices),
                     set(ag[getattr(ag, self.plurals[prop]) != 500.0].indices))

    def test_props(self):
        u = mda.Universe(GRO)

        for prop in ['mass', 'charge']:
            for ag in [u.atoms, u.atoms[:100]]:
                yield self._check_lt, prop, ag
                yield self._check_le, prop, ag
                yield self._check_gt, prop, ag
                yield self._check_ge, prop, ag
                yield self._check_eq, prop, ag
                yield self._check_ne, prop, ag


class TestBondedSelection(object):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    def test_bonded_1(self):
        ag = self.u.select_atoms('type 2 and bonded name N')

        assert_(len(ag) == 3)

    def test_nobonds_warns(self):
        self.u.bonds = TopologyGroup([])

        assert_warns(UserWarning,
                     self.u.select_atoms, 'type 2 and bonded name N')


class TestSelectionErrors(object):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    def selection_fail(self, selstr):
        assert_raises(SelectionError, self.u.select_atoms,
                      selstr)

    def test_expected_errors(self):
        for selstr in [
                'name and H',  # string selection
                'name )',
                'resid abcd',  # range selection
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
        ]:
            yield self.selection_fail, selstr


def test_segid_and_resid():
    u = mda.Universe(PDB_full)

    ag = u.select_atoms('segid B and resid 1-100')

    ref = ag.select_atoms('segid B').select_atoms('resid 1-100')

    assert_array_equal(ag.indices, ref.indices)


class TestImplicitOr(object):
    def setUp(self):
        self.u = mda.Universe(PSF)

    def tearDown(self):
        del self.u

    def _check_sels(self, ref, sel):
        ref = self.u.select_atoms(ref)
        sel = self.u.select_atoms(sel)

        assert_array_equal(ref.indices, sel.indices)

    def test_string_selections(self):
        for ref, sel in (
                ('name HT1 or name HT2 or name HT3', 'name HT1 HT2 HT3'),
                ('type 2 or type 3 or type 4', 'type 2 3 4'),
                ('resname MET or resname GLY', 'resname MET GLY'),
                ('name H* or name N', 'name H* N'),
                ('(name N or name HT1) and (resname MET or resname GLY)',
                 'name N HT1 and resname MET GLY'),
        ):
            yield self._check_sels, ref, sel

    def test_segids(self):
        # Don't put me into generator, I use a different Universe
        # (With multiple segids)
        u = mda.Universe(PDB_full)
        ref = u.select_atoms('segid A or segid B')
        ag = u.select_atoms('segid A B')

        assert_array_equal(ref.indices, ag.indices)

    def test_range_selections(self):
        # All these selections just use numeric types,
        # So loop over what type of selections,
        # And apply the same numeric constraints to all
        for seltype in ['resid', 'resnum', 'bynum']:
            for ref, sel in (
                    ('{typ} 1 or {typ} 2', '{typ} 1 2'),
                    ('{typ} 1:10 or {typ} 22', '{typ} 1:10 22'),
                    ('{typ} 1:10 or {typ} 20:30', '{typ} 1:10 20:30'),
                    ('{typ} 1-5 or {typ} 7', '{typ} 1-5 7'),
                    ('{typ} 1-5 or {typ} 7:10 or {typ} 12',
                     '{typ} 1-5 7:10 12'),
                    ('{typ} 1 or {typ} 3 or {typ} 5:10', '{typ} 1 3 5:10'),
            ):
                yield (self._check_sels,
                       ref.format(typ=seltype),
                       sel.format(typ=seltype))
