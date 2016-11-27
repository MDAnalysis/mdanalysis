# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
from six.moves import zip

import MDAnalysis
from MDAnalysis.tests.datafiles import (PSF, DCD, PDB_small, GRO, TRR,
                                        TRZ, TRZ_psf, PSF_notop,
                                        XYZ_mini, two_water_gro,
                                        two_water_gro_nonames,
                                        COORDINATES_XYZ, COORDINATES_TRR)
import MDAnalysis.core.groups
from MDAnalysis.core.groups import Atom, AtomGroup
from MDAnalysis import NoDataError

import numpy as np
from numpy.testing import (TestCase, dec, raises, assert_equal,
                           assert_almost_equal, assert_raises, assert_,
                           assert_array_almost_equal, assert_array_equal,
                           assert_allclose)
from nose.plugins.attrib import attr

from unittest import skip

import os
import itertools

from MDAnalysisTests import parser_not_found, tempdir
from MDAnalysisTests.core.groupbase import make_Universe


class TestAtom(TestCase):
    """Tests of Atom."""

    # VALID
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.atom = self.universe.atoms[1000]  # Leu67:CG
        self.known_pos = np.array([3.94543672, -12.4060812, -7.26820087],
                                  dtype=np.float32)

    # VALID
    def tearDown(self):
        del self.universe
        del self.atom
        del self.known_pos

    # VALID
    def test_attributes_names(self):
        a = self.atom
        assert_equal(a.name, 'CG')
        assert_equal(a.resname, 'LEU')

    def test_setting_attribute_name(self):
        self.atom.name = 'AA'
        assert_equal(self.atom.name, 'AA')

    def test_setting_attribute_type(self):
        self.atom.type = 'Z'
        assert_equal(self.atom.type, 'Z')

    #INVALID
    @skip
    def test_setting_attribute_mass(self):
        assert_equal(self.atom.mass, 13)

    def test_setting_attributes_charge(self):
        self.atom.charge = 6
        assert_equal(self.atom.charge, 6)

    # VALID
    def test_attributes_positions(self):
        a = self.atom
        # new position property (mutable)
        assert_almost_equal(a.position, self.known_pos)
        pos = a.position + 3.14
        a.position = pos
        assert_almost_equal(a.position, pos)

    # VALID: does this selection style generalize under new topology system"
    def test_atom_selection(self):
        asel = self.universe.select_atoms('atom 4AKE 67 CG').atoms[0]
        assert_equal(self.atom, asel)

    # VALID
    def test_hierarchy(self):
        u = self.universe
        a = self.atom
        assert_equal(a.segment, u.s4AKE)
        assert_equal(a.residue, u.residues[66])

    # VALID: shows we need to raise TypError here
    @skip
    def test_bad_add(self):
        def bad_add():
            return self.atom + 1

        assert_raises(TypeError, bad_add)

    # VALID
    def test_add_AG(self):
        ag = self.universe.atoms[:2]

        ag2 = self.atom + ag

        for at in [self.atom, ag[0], ag[1]]:
            assert_equal(at in ag2, True)

    # VALID
    def test_no_velo(self):
        def lookup_velo():
            return self.atom.velocity

        assert_raises(NoDataError, lookup_velo)

    def test_bonded_atoms(self):
        at = self.universe.atoms[0]
        ref = [b.partner(at) for b in at.bonds]
        assert_equal(ref, list(at.bonded_atoms))

    # VALID
    @raises(AttributeError)
    def test_undefined_occupancy(self):
        self.universe.atoms[0].occupancy

    # INVALID: Atom objects are created when doing, e.g. self.universe.atoms[0].
    # attributes are not propagated to tables of attributes in this way
    @skip
    def test_set_undefined_occupancy(self):
        self.universe.atoms[0].occupancy = .5
        assert_equal(self.universe.atoms[0].occupancy, .5)
        assert_equal(self.universe.atoms[1].occupancy, 1)


# VALID: these attributes exist always and are pulled from trajectory;
# they give NoDataError if the reader doesn't give this information
class TestAtomNoForceNoVel(TestCase):
    def setUp(self):
        self.u = MDAnalysis.Universe(XYZ_mini)
        self.a = self.u.atoms[0]

    def tearDown(self):
        del self.u

    def test_velocity_fail(self):
        assert_raises(NoDataError, getattr, self.a, 'velocity')

    def test_force_fail(self):
        assert_raises(NoDataError, getattr, self.a, 'force')

    def test_velocity_set_fail(self):
        assert_raises(NoDataError, setattr, self.a, 'velocity',
                      [1.0, 1.0, 1.0])

    def test_force_set_fail(self):
        assert_raises(NoDataError, setattr, self.a, 'force', [1.0, 1.0, 1.0])


class TestAtomGroup(TestCase):
    """Tests of AtomGroup; selections are tested separately."""
    # all tests are done with the AdK system (PSF and DCD) sequence:
    # http://www.uniprot.org/uniprot/P69441.fasta
    # >sp|P69441|KAD_ECOLI Adenylate kinase OS=Escherichia coli (strain K12) GN=adk PE=1 SV=1
    ref_adk_sequence = (
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT"
        "DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI"
        "VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG"
        "YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"
    )

    # VALID
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.ag = self.universe.atoms  # prototypical AtomGroup
        self.dih_prec = 2

    # INVALID: cannot use class in core to build AtomGroups for any Universe
    @skip
    def test_newAtomGroup(self):
        newag = MDAnalysis.core.groups.AtomGroup(self.ag[1000:2000:200])
        assert_equal(type(newag), type(self.ag),
                     "Failed to make a new AtomGroup: type mismatch")
        assert_equal(newag.n_atoms, len(self.ag[1000:2000:200]))
        assert_equal(newag.n_residues, 5)
        # check any special method
        assert_almost_equal(newag.total_mass(), 40.044999999999995)


    def test_getitem_int(self):
        assert_equal(self.universe.atoms[0].ix, self.universe.atoms.ix[0])

    def test_getitem_slice(self):
        assert_array_equal(self.universe.atoms[0:4].ix,
                           self.universe.atoms.ix[:4])

    def test_getitem_slice2(self):
        assert_equal(self.universe.atoms[0:8:2].ix,
                     self.universe.atoms.ix[0:8:2])

    # INVALID: we don't support getitem with names anymore.
    # it could be supported by making the Atomnames topologyattr transplant
    # a wrapped __getitem__ onto AtomGroup
    @skip
    def test_getitem_str(self):
        ag1 = self.universe.atoms['HT1']
        # select_atoms always returns an AtomGroup even if single result
        ag2 = self.universe.select_atoms('name HT1')[0]
        assert_equal(ag1, ag2)

    # INVALID: can't use `AtomGroup` class directly;
    # getitem with list works fine though
    @skip
    def test_getitem_list(self):
        sel = [0, 1, 4]
        ag1 = self.universe.atoms[sel]
        ag2 = AtomGroup([self.universe.atoms[i] for i in sel])
        assert_equal(ag1._atoms, ag2._atoms)

    # INVALID: can't use `AtomGroup` class directly;
    # getitem with array works fine though
    @skip
    def test_getitem_nparray(self):
        sel = np.arange(5)
        ag1 = self.universe.atoms[sel]
        ag2 = AtomGroup([self.universe.atoms[i] for i in sel])
        assert_equal(ag1._atoms, ag2._atoms)

    # VALID: we should raise TypeError for this; currently get IndexError
    @skip
    def test_getitem_TE(self):
        d = {'A': 1}
        raises(TypeError, self.universe.atoms.__getitem__, d)

    # INVALID: can't build AtomGroup directly like this
    @skip
    def test_bad_make(self):
        raises(TypeError, AtomGroup, ['these', 'are', 'not', 'atoms'])

    # VALID
    def test_n_atoms(self):
        assert_equal(self.ag.n_atoms, 3341)

    # VALID
    def test_n_residues(self):
        assert_equal(self.ag.n_residues, 214)

    # VALID
    def test_n_segments(self):
        assert_equal(self.ag.n_segments, 1)

    # VALID
    def test_resids_dim(self):
        assert_equal(len(self.ag.resids), len(self.ag))

    # INVALID: this topology doesn't include resnums, so we don't have any
    @skip
    def test_resnums_dim(self):
        assert_equal(len(self.ag.resnums), len(self.ag))

    # VALID
    def test_segids_dim(self):
        assert_equal(len(self.ag.segids), len(self.ag))

    # VALID
    def test_len(self):
        """testing that len(atomgroup) == atomgroup.n_atoms"""
        assert_equal(len(self.ag), self.ag.n_atoms, "len and n_atoms disagree")

    def test_center(self):
        weights = np.zeros(self.ag.n_atoms)
        weights[0] = 1
        assert_array_almost_equal(self.ag.center(weights),
                                  self.ag.positions[0])

        weights[:4] = 1. / 4.
        assert_array_almost_equal(self.ag.center(weights),
                                  self.ag.positions[:4].mean(axis=0))

    # VALID
    def test_center_of_geometry(self):
        assert_array_almost_equal(self.ag.center_of_geometry(),
                                  np.array([-0.04223963, 0.0141824,
                                            -0.03505163], dtype=np.float32))

    # INVALID: AtomGroup only has center_of_mass method if topology has masses
    @skip
    def test_center_of_mass(self):
        assert_array_almost_equal(self.ag.center_of_mass(),
                                  np.array([-0.01094035, 0.05727601,
                                            -0.12885778]))

    # VALID
    def test_coordinates(self):
        assert_array_almost_equal(
            self.ag.positions[1000:2000:200],
            np.array([[3.94543672, -12.4060812, -7.26820087],
                      [13.21632767, 5.879035, -14.67914867],
                      [12.07735443, -9.00604534, 4.09301519],
                      [11.35541916, 7.0690732, -0.32511973],
                      [-13.26763439, 4.90658951, 10.6880455]],
                     dtype=np.float32))

    # VALID
    def test_principal_axes(self):
        assert_array_almost_equal(
            self.ag.principal_axes(),
            np.array([[-9.99925632e-01, 1.21546132e-02, 9.98264877e-04],
                      [1.20986911e-02, 9.98951474e-01, -4.41539838e-02],
                      [1.53389276e-03, 4.41386224e-02, 9.99024239e-01]]))

    # VALID: need to set precision to correct within 5 decimal points
    @skip
    def test_total_charge(self):
        assert_almost_equal(self.ag.total_charge(), -4.0)

    # VALID: need to set precision to correct within 2 decimal points
    # perhaps we should use a higher precision float for masses in
    # PSF parser
    @skip
    def test_total_mass(self):
        assert_almost_equal(self.ag.total_mass(), 23582.043)

    # VALID
    def test_indices_ndarray(self):
        assert_equal(isinstance(self.ag.indices, np.ndarray), True)

    # VALID
    def test_indices(self):
        assert_array_equal(self.ag.indices[:5], np.array([0, 1, 2, 3, 4]))

    # VALID
    def test_resids_ndarray(self):
        assert_equal(isinstance(self.ag.resids, np.ndarray), True)

    # VALID
    def test_resids(self):
        assert_array_equal(self.ag.residues.resids, np.arange(1, 215))

    # INVALID: this topology doesn't give resnums, so no such property exists
    @skip
    def test_resnums_ndarray(self):
        assert_equal(isinstance(self.ag.residues.resnums, np.ndarray), True)

    # INVALID: see above
    @skip
    def test_resnums(self):
        assert_array_equal(self.ag.residues.resnums, np.arange(1, 215))

    # VALID
    def test_resnames_ndarray(self):
        assert_equal(isinstance(self.ag.residues.resnames, np.ndarray), True)

    # VALID
    def test_resnames(self):
        resnames = self.ag.residues.resnames
        assert_array_equal(resnames[0:3], np.array(["MET", "ARG", "ILE"]))

    # VALID
    def test_names_ndarray(self):
        assert_equal(isinstance(self.ag.names, np.ndarray), True)

    # VALID
    def test_names(self):
        names = self.ag.names
        assert_array_equal(names[0:3], np.array(["N", "HT1", "HT2"]))

    # VALID
    def test_segids_ndarray(self):
        assert_equal(isinstance(self.ag.segids, np.ndarray), True)

    # VALID
    def test_segids(self):
        segids = self.ag.segids
        assert_array_equal(segids[0], np.array(["4AKE"]))

    # INVALID: this topology doesn't give masses
    @skip
    def test_masses_ndarray(self):
        assert_equal(isinstance(self.ag.masses, np.ndarray), True)

    # INVALID: this topology doesn't give masses
    @skip
    def test_masses(self):
        masses = self.ag.masses
        assert_array_equal(masses[0:3], np.array([14.007, 1.008, 1.008]))

    # VALID
    def test_charges_ndarray(self):
        assert_equal(isinstance(self.ag.charges, np.ndarray), True)

    # VALID
    def test_charges(self):
        assert_array_almost_equal(self.ag.charges[1000:2000:200],
                                  np.array([-0.09, 0.09, -0.47, 0.51, 0.09]))

    # INVALID: this topology doesn't give atom radii
    @skip
    def test_radii_ndarray(self):
        assert_equal(isinstance(self.ag.radii, np.ndarray), True)

    # INVALID: this topology doesn't give atom radii
    @skip
    def test_radii(self):
        radii = self.ag.radii
        assert_array_equal(radii[0:3], np.array([None, None, None]))

    # INVALID: this topology doesn't give bfactors
    @skip
    def test_bfactors_ndarray(self):
        assert_equal(isinstance(self.ag.bfactors, np.ndarray), True)

    # INVALID: this topology doesn't give bfactors
    @skip
    def test_bfactors(self):
        bfactors = self.ag.bfactors  # property, not method!
        assert_array_equal(bfactors[0:3], np.array([None, None, None]))

    # INVALID: we give AttributeError since occupancies attribute not present
    # if topology doesn't provide them
    @skip
    def test_occupancies(self):
        assert_raises(NoDataError, getattr, self.ag, 'occupancies')
        self.ag.occupancies = 0.25
        assert_array_almost_equal(self.ag.occupancies,
                                  np.ones(len(self.ag)) * 0.25)

    # INVALID: AtomGroup no longer gets a sequence method; only in ResidueGroup
    # if resnames in topology
    @skip
    def test_sequence_from_atoms(self):
        p = self.universe.select_atoms("protein")
        assert_equal(p.sequence(format="string"),
                     p.residues.sequence(format="string"),
                     err_msg="sequence() yields different results for "
                     "residues and atoms")

    # VALID
    def test_sequence_string(self):
        p = self.universe.select_atoms("protein")
        assert_equal(p.residues.sequence(format="string"),
                     self.ref_adk_sequence)

    # VALID
    def test_sequence_SeqRecord(self):
        p = self.universe.select_atoms("protein")
        s = p.residues.sequence(format="SeqRecord",
                                id="P69441", name="KAD_ECOLI Adenylate kinase",
                                description="EcAdK from pdb 4AKE")
        assert_equal(s.id, "P69441")
        assert_equal(s.seq.tostring(), self.ref_adk_sequence)

    # VALID
    def test_sequence_SeqRecord_default(self):
        p = self.universe.select_atoms("protein")
        s = p.residues.sequence(id="P69441", name="KAD_ECOLI Adenylate kinase",
                                description="EcAdK from pdb 4AKE")
        assert_equal(s.id, "P69441")
        assert_equal(s.seq.tostring(), self.ref_adk_sequence)

    # VALID
    def test_sequence_Seq(self):
        p = self.universe.select_atoms("protein")
        s = p.residues.sequence(format="Seq")
        assert_equal(s.tostring(), self.ref_adk_sequence)

    # INVALID: currently raises AttributeError
    @skip
    def test_sequence_nonIUPACresname(self):
        """test_sequence_nonIUPACresname: non recognized amino acids raise
        ValueError"""
        # fake non-IUPAC residue name for this test
        residues = self.universe.select_atoms("resname MET").residues
        residues.resnames = "MSE"

        def wrong_res():
            self.universe.atoms.sequence()

        assert_raises(ValueError, wrong_res)

    # INVALID: cannot create Atoms without a Universe
    @skip
    def test_no_uni_1(self):
        at1 = Atom(1, 'dave', 'C', 'a', 1, 1, 0.1, 0.0)
        at2 = Atom(2, 'dave', 'C', 'a', 1, 1, 0.1, 0.0)
        ag = AtomGroup([at1, at2])
        assert_raises(NoDataError, getattr, ag, 'universe')

    # INVALID: cannot create AtomGroups without a Universe
    @skip
    def test_no_uni_2(self):
        ag = AtomGroup([])
        assert_raises(NoDataError, getattr, ag, 'universe')

    # VALID: AtomGroup should raise TypeError when adding with non-AtomGroup
    @skip
    def test_bad_add_AG(self):
        def bad_add():
            return self.ag + [1, 2, 3]
        assert_raises(TypeError, bad_add)

    # INVALID: cannot create AtomGroup without a Universe
    @skip
    def test_bool_false(self):
        # Issue #304
        ag = AtomGroup([])
        assert_equal(bool(ag), False)

    # VALID
    def test_bool_true(self):
        # Issue #304
        ag = self.universe.atoms[:3]
        assert_equal(bool(ag), True)

    # VALID
    def test_repr(self):
        # Should make sure that the user facing info stays as expected
        assert_equal(repr(self.ag), "<AtomGroup with 3341 atoms>")

    # INVALID: set resnums with `ag.residues.resnums = new`
    # also resnums won't be present with this topology, perhaps
    # Issue 202 following 4 tests
    @skip
    def test_set_resnum_single(self):
        ag = self.universe.atoms[:3]
        new = 5
        ag.set_resnums(new)
        for at in ag:
            assert_equal(at.resnum, new)
        assert_equal(all(ag.resnums == new), True)

    # INVALID: set resnums with `ag.residues.resnums = new`
    # also resnums won't be present with this topology, perhaps
    @skip
    def test_set_resnum_many(self):
        ag = self.universe.atoms[:3]
        new = [22, 23, 24]
        ag.set_resnums(new)
        for at, v in zip(ag, new):
            assert_equal(at.resnum, v)
        assert_equal(all(ag.resnums == new), True)

    def test_set_resname_single(self):
        ag = self.universe.atoms[:3]
        new = 'abc'
        ag.residues.resnames = new
        for at in ag:
            assert_equal(at.resname, new)
        assert_equal(all(ag.resnames == new), True)

    # INVALID: set resnames with `ag.residues.resnames = new`
    @skip
    def test_set_resname_many(self):
        ag = self.universe.atoms[:3]
        new = ['aa', 'bb', 'cc']
        ag.set_resnames(new)
        for at, v in zip(ag, new):
            assert_equal(at.resname, v)
        assert_equal(all(ag.resnames == new), True)

    # TODO: add all other methods except select_atoms(), see
    # test_atomselections.py
    # INVALID: cannot create Atoms like this; must come from a Universe
    @skip
    def test_set_charges(self):
        # Charges are initially 0
        at1 = Atom(1, 'dave', 'C', 'a', 1, 1, 0.1, 0.0)
        at2 = Atom(2, 'dave', 'C', 'a', 1, 1, 0.1, 0.0)
        ag = AtomGroup([at1, at2])

        charges = [1.0, 2.0]
        ag.set_charges(charges)
        for at, val in zip(ag, charges):
            assert_equal(at.charge, val)

    # INVALID: cannot create Atoms like this; must come from a Universe
    @skip
    def test_set_radii(self):
        at1 = Atom(1, 'dave', 'C', 'a', 1, 1, 0.1, 0.0)
        at2 = Atom(2, 'dave', 'C', 'a', 1, 1, 0.1, 0.0)
        ag = AtomGroup([at1, at2])

        radii = [1.0, 2.0]
        ag.set_radii(radii)
        for at, val in zip(ag, radii):
            assert_equal(at.radius, val)

    # INVALID: cannot create Atoms like this; must come from a Universe
    @skip
    def test_set_bfactors(self):
        at1 = Atom(1, 'dave', 'C', 'a', 1, 1, 0.1, 0.0)
        at2 = Atom(2, 'dave', 'C', 'a', 1, 1, 0.1, 0.0)
        ag = AtomGroup([at1, at2])

        bfacs = [1.0, 2.0]
        ag.set_bfactors(bfacs)
        for at, val in zip(ag, bfacs):
            assert_equal(at.bfactor, val)

    # add new methods here...
    # VALID
    def test_packintobox_badshape(self):
        ag = self.universe.atoms[:10]
        box = np.zeros(9, dtype=np.float32).reshape(3, 3)

        def badpack(a):
            return a.pack_into_box(box=box)

        assert_raises(ValueError, badpack, ag)

    # VALID
    def test_packintobox_noshape(self):
        ag = self.universe.atoms[:10]

        def badpack(a):
            return a.pack_into_box()

        assert_raises(ValueError, badpack, ag)

    # VALID
    def test_packintobox(self):
        """test AtomGroup.pack_into_box(): Tests application of periodic boundary
        conditions on coordinates

        Reference system doesn't have dimensions, so an arbitrary box is
        imposed on the system

        """
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        ag = u.atoms[1000:2000:200]
        # Provide arbitrary box
        ag.pack_into_box(box=np.array([5., 5., 5.], dtype=np.float32))
        assert_array_almost_equal(
            ag.positions,
            np.array([[3.94543672, 2.5939188, 2.73179913],
                      [3.21632767, 0.879035, 0.32085133],
                      [2.07735443, 0.99395466, 4.09301519],
                      [1.35541916, 2.0690732, 4.67488003],
                      [1.73236561, 4.90658951, 0.6880455]], dtype=np.float32))

    # INVALID: no `_atoms` property for ResidueGroup or AtomGroup; use `.atoms`
    # which won't work because they are compared by reference in the default.
    # But is also doesn't work to use `ix` as they will differ.
    @skip
    def test_residues(self):
        u = self.universe
        assert_equal(u.residues[100]._atoms,
                     u.select_atoms('resname ILE and resid 101')._atoms,
                     "Direct selection from residue group does not match "
                     "expected I101.")

    def test_segments(self):
        u = self.universe
        assert_equal(len(u.segments.s4AKE.atoms),
                     len(u.select_atoms('segid 4AKE').atoms),
                     "Direct selection of segment 4AKE from segments failed.")

    # VALID
    def test_index_integer(self):
        u = self.universe
        a = u.atoms[100]
        assert_(isinstance(a, Atom), "integer index did not return Atom")

    # VALID
    def test_index_slice(self):
        u = self.universe
        a = u.atoms[100:200:10]
        assert_(isinstance(a, AtomGroup),
                "slice index did not return AtomGroup")

    # VALID
    def test_index_slice_empty(self):
        u = self.universe
        assert_array_equal(u.atoms[0:0], [],
                           "making an empty AtomGroup failed")

    # VALID:
    def test_index_advancedslice(self):
        u = self.universe
        aslice = [0, 10, 20, -1, 10]
        ag = u.atoms[aslice]
        assert_(isinstance(ag, AtomGroup),
                "advanced slicing does not produce a AtomGroup")
        assert_equal(ag[1], ag[-1], "advanced slicing does not preserve order")

    # VALID
    def test_boolean_indexing(self):
        # index an array with a sequence of bools
        # issue #282
        sel = np.array([True, False, True])
        ag = self.universe.atoms[10:13]
        ag2 = ag[sel]
        assert_equal(len(ag2), 2)
        for at in [ag[0], ag[2]]:
            assert_equal(at in ag2, True)

    # VALID
    def test_boolean_indexing_2(self):
        # index an array with a sequence of bools
        # issue #282
        sel = [True, False, True]
        ag = self.universe.atoms[10:13]
        ag2 = ag[sel]
        assert_equal(len(ag2), 2)
        for at in [ag[0], ag[2]]:
            assert_equal(at in ag2, True)

    # VALID
    def test_bool_IE(self):
        # indexing with empty list doesn't run foul of bool check
        sel = []
        ag = self.universe.atoms[10:30]
        ag2 = ag[sel]
        assert_equal(len(ag2), 0)

    # VALID
    def test_dihedral_ValueError(self):
        """test that AtomGroup.dihedral() raises ValueError if not exactly
        4 atoms given"""
        nodih = self.universe.select_atoms("resid 3:10")
        assert_raises(ValueError, getattr, nodih, 'dihedral')
        nodih = self.universe.select_atoms("resid 3:5")
        assert_raises(ValueError, getattr, nodih, 'dihedral')

    # VALID
    def test_improper(self):
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        peptbond = u.select_atoms("atom 4AKE 20 C", "atom 4AKE 21 CA",
                                  "atom 4AKE 21 N", "atom 4AKE 21 HN")
        assert_almost_equal(peptbond.improper.value(), 168.52952575683594,
                            self.dih_prec,
                            "Peptide bond improper dihedral for M21 "
                            "calculated wrongly.")

    # VALID
    def test_dihedral_equals_improper(self):
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        peptbond = u.select_atoms("atom 4AKE 20 C", "atom 4AKE 21 CA",
                                  "atom 4AKE 21 N", "atom 4AKE 21 HN")
        assert_equal(peptbond.improper.value(), peptbond.dihedral.value(),
                     "improper() and proper dihedral() give different results")

    # INVALID: we don't support getting residues by resid as attributes of segments anymore
    @skip
    def test_bond(self):
        self.universe.trajectory.rewind()  # just to make sure...
        sel2 = self.universe.s4AKE.r98.select_atoms("name OE1", "name OE2")
        assert_almost_equal(sel2.bond.value(), 2.1210737228393555, 3,
                            "distance of Glu98 OE1--OE2 wrong")

    # INVALID: we don't support getting residues by resid as attributes of segments anymore
    @skip
    def test_bond_pbc(self):
        self.universe.trajectory.rewind()
        sel2 = self.universe.s4AKE.r98.select_atoms("name OE1", "name OE2")
        assert_almost_equal(sel2.bond.value(pbc=True), 2.1210737228393555, 3,
                            "distance of Glu98 OE1--OE2 wrong")

    # VALID
    def test_bond_ValueError(self):
        ag = self.universe.atoms[:4]
        assert_raises(ValueError, getattr, ag, 'bond')

    # INVALID: we don't support getting residues by resid as attributes of segments anymore
    @skip
    def test_angle(self):
        self.universe.trajectory.rewind()  # just to make sure...
        sel3 = self.universe.s4AKE.r98.select_atoms("name OE1", "name CD",
                                                    "name OE2")
        assert_almost_equal(sel3.angle.value(), 117.46187591552734, 3,
                            "angle of Glu98 OE1-CD-OE2 wrong")

    # VALID
    def test_angle_ValueError(self):
        ag = self.universe.atoms[:2]
        assert_raises(ValueError, getattr, ag, 'angle')

    # INVALID: this topology has no masses, so no `shape_parameter` method
    @skip
    def test_shape_parameter(self):
        s = self.universe.s4AKE.shape_parameter()
        assert_almost_equal(s, 0.00240753939086033, 6)

    # INVALID: this topology has no masses, so no `asphericity` method
    @skip
    def test_asphericity(self):
        a = self.universe.s4AKE.asphericity()
        assert_almost_equal(a, 0.020227504542775828, 6)

    # TODO: tests for the coordinate manipulation methods
    # - transform
    # - translate
    # - rotate
    # - rotateby

    # VALID
    def test_positions(self):
        ag = self.universe.select_atoms("bynum 12:42")
        pos = ag.positions + 3.14
        ag.positions = pos
        # should work
        assert_almost_equal(ag.positions, pos,
                            err_msg="failed to update atoms 12:42 position "
                            "to new position")

        def set_badarr(pos=pos):
            # create wrong size array
            badarr = np.random.random((pos.shape[0] - 1, pos.shape[1] - 1))
            ag.positions = badarr
        assert_raises(ValueError, set_badarr)

    # INVALID: no `get_positions` method for AtomGroup; use `positions` property instead
    @skip
    def test_set_positions(self):
        ag = self.universe.select_atoms("bynum 12:42")
        pos = ag.get_positions() + 3.14
        ag.set_positions(pos)
        assert_almost_equal(ag.positions, pos,
                            err_msg="failed to update atoms 12:42 position "
                            "to new position")

    # INVALID: no `get_velocities` method for AtomGroup; use `velocities` property instead
    @skip
    def test_no_velocities_raises_NoDataError(self):
        def get_vel(ag=self.universe.select_atoms("bynum 12:42")):
            ag.get_velocities()
        # trj has no velocities
        assert_raises(NoDataError, get_vel)

    # INVALID: no `set_velocities` method for AtomGroup; use `velocities` property instead
    @skip
    def test_set_velocities_NoData(self):
        def set_vel():
            return self.universe.atoms[:2].set_velocities([0.2])
        assert_raises(NoDataError, set_vel)

    # INVALID: no `get_forces` method for AtomGroup; use `forces` property instead
    @skip
    def test_get_forces_NoData(self):
        def get_for():
            return self.universe.atoms[:2].get_forces()
        assert_raises(NoDataError, get_for)

    # INVALID: no `set_forces` method for AtomGroup; use `forces` property instead
    @skip
    def test_set_forces_NoData(self):
        def set_for():
            return self.universe.atoms[:2].set_forces([0.2])
        assert_raises(NoDataError, set_for)

    # INVALID: no `set_resids` method for AtomGroup; use
    # `AtomGroup.residues.resids` property instead
    @skip
    def test_set_resids(self):
        ag = self.universe.select_atoms("bynum 12:42")
        resid = 999
        ag.set_resids(resid)
        # check individual atoms
        assert_equal([a.resid for a in ag],
                     resid * np.ones(ag.n_atoms),
                     err_msg="failed to set_resid atoms 12:42 to same resid")
        # check residues
        assert_equal(ag.resids, 999 * np.ones(ag.n_residues),
                     err_msg="failed to set_resid of residues belonging to "
                     "atoms 12:42 to same resid")

    # INVALID: no `set_names` method for AtomGroup; use `names` property instead
    @skip
    def test_set_names(self):
        ag = self.universe.atoms[:2]
        names = ['One', 'Two']
        ag.set_names(names)
        for a, b in zip(ag, names):
            assert_equal(a.name, b)

    # INVALID: change name to not match name above
    # also, no `set_resids` method of AtomGroups; use
    # `AtomGroup.residues.resids` property instead
    @skip
    def test_set_resids(self):
        """test_set_resid: set AtomGroup resids on a per-atom basis"""
        ag = self.universe.select_atoms("bynum 12:42")
        resids = np.array([a.resid for a in ag]) + 1000
        ag.set_resids(resids)
        # check individual atoms
        assert_equal([a.resid for a in ag], resids,
                     err_msg="failed to set_resid atoms 12:42 to "
                     "resids {0}".format(resids))
        # check residues
        assert_equal(ag.residues.resids, np.unique(resids),
                     err_msg="failed to set_resid of residues belonging to "
                     "atoms 12:42 to same resid")

    # INVALID: setting resids does not set residue membership for atoms.
    # must change an atom's resindex, which is guaranteed to be unique for each
    # residue
    # also no `set_resids` method for AtomGroups
    @skip
    def test_merge_residues(self):
        ag = self.universe.select_atoms("resid 12:14")
        nres_old = self.universe.atoms.n_residues
        natoms_old = ag.n_atoms
        ag.set_resids(12)  # merge all into one with resid 12
        nres_new = self.universe.atoms.n_residues
        r_merged = self.universe.select_atoms("resid 12:14").residues
        natoms_new = self.universe.select_atoms("resid 12").n_atoms
        assert_equal(len(r_merged), 1,
                     err_msg="set_resid failed to merge residues: "
                     "merged = {0}".format(r_merged))
        assert_equal(nres_new, nres_old - 2,
                     err_msg="set_resid failed to merge residues: "
                     "merged = {0}".format(r_merged))
        assert_equal(natoms_new, natoms_old,
                     err_msg="set_resid lost atoms on merge".format(r_merged))

    # INVALID: no `set_masses` method for AtomGroup. `AtomGroup.masses` property only exists
    # if topology already has it, or if the topologyattr for masses is added
    @skip
    def test_set_masses(self):
        ag = self.universe.select_atoms("bynum 12:42 and name H*")
        mass = 2.0
        ag.set_masses(mass)
        # check individual atoms
        assert_equal([a.mass for a in ag],
                     mass * np.ones(ag.n_atoms),
                     err_msg="failed to set_mass H* atoms in resid "
                     "12:42 to {0}".format(mass))

    # INVALID: no `set_segids` method of AtomGroup. Residues are members
    # of segments, so one must do `AtomGroup.residues.segindices = u.LID.segindex`
    @skip
    def test_set_segids(self):
        u = self.universe
        u.select_atoms("(resid 1-29 or resid 60-121 or resid 160-214)").set_segids("CORE")
        u.select_atoms("resid 122-159").set_segids("LID")
        u.select_atoms("resid 30-59").set_segids("NMP")
        assert_equal(u.segments.segids, ["CORE", "NMP", "CORE", "LID", "CORE"],
                     err_msg="failed to change segids = {0}".format(u.atoms.segids))

    # INVALID: no `set_masses` method for AtomGroup. `AtomGroup.masses` property only exists
    # if topology already has it, or if the topologyattr for masses is added
    @skip
    def test_wronglen_set(self):
        """Give the setter function a list of wrong length"""
        assert_raises(ValueError, self.ag.set_masses, [0.1, 0.2])

    # VALID
    def test_split_atoms(self):
        ag = self.universe.select_atoms("resid 1:50 and not resname LYS and "
                                        "(name CA or name CB)")
        sg = ag.split("atom")
        assert_equal(len(sg), len(ag))
        for g, ref_atom in zip(sg, ag):
            atom = g[0]
            assert_equal(len(g), 1)
            assert_equal(atom.name, ref_atom.name)
            assert_equal(atom.resid, ref_atom.resid)

    # VALID
    def test_split_residues(self):
        ag = self.universe.select_atoms("resid 1:50 and not resname LYS and "
                                        "(name CA or name CB)")
        sg = ag.split("residue")
        assert_equal(len(sg), len(ag.residues.resids))
        for g, ref_resname in zip(sg, ag.residues.resnames):
            if ref_resname == "GLY":
                assert_equal(len(g), 1)
            else:
                assert_equal(len(g), 2)
            for atom in g:
                assert_equal(atom.resname, ref_resname)

    # VALID
    def test_split_segments(self):
        ag = self.universe.select_atoms("resid 1:50 and not resname LYS and "
                                        "(name CA or name CB)")
        sg = ag.split("segment")
        assert_equal(len(sg), len(ag.segments.segids))
        for g, ref_segname in zip(sg, ag.segments.segids):
            for atom in g:
                assert_equal(atom.segid, ref_segname)

    # VALID
    # instant selectors
    @attr("issue")
    def test_nonexistent_instantselector_raises_AttributeError(self):
        def access_nonexistent_instantselector():
            self.universe.atoms.NO_SUCH_ATOM
        assert_raises(AttributeError, access_nonexistent_instantselector)

    # VALID
    def test_atom_order(self):
        assert_equal(self.universe.atoms.indices,
                     sorted(self.universe.atoms.indices))


class TestAtomGroupNoTop(TestCase):

    # VALID
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.u = MDAnalysis.Universe(PSF_notop, DCD)
        self.ag = self.u.atoms[:10]

    # VALID
    def tearDown(self):
        del self.u
        del self.ag

    # INVALID: no bonds attribute if not bonds in topology
    @skip
    def test_nobonds(self):
        assert_equal(self.ag.bonds, [])

    # INVALID: no angles attribute if not bonds in topology
    @skip
    def test_noangles(self):
        assert_equal(self.ag.angles, [])

    # INVALID: no dihedrals attribute if not bonds in topology
    @skip
    def test_nodihedrals(self):
        assert_equal(self.ag.dihedrals, [])

    # INVALID: no impropers attribute if not bonds in topology
    @skip
    def test_noimps(self):
        assert_equal(self.ag.impropers, [])

    # Because I'm messing with atom info, I've put these here separated from
    # other tests
    # INVALID: AtomGroups keep no caches now
    @skip
    def test_clear_cache(self):
        self.ag._clear_caches()

        assert_equal(self.ag._cache, dict())

    # INVALID: no caches in AtomGroups anymore
    @skip
    def test_rebuild_cache_residues(self):
        assert_equal(len(self.ag.residues), 1)

        # Mess with stuff, add a different residues and segment for the first
        # atom
        self.ag[0].residue = self.u.atoms[100].residue

        # There's actually 2 residues now, but because of cache this isn't
        # detected
        assert_equal(len(self.ag.residues), 1)

        # After cache rebuild second residue is finally seen
        self.ag._rebuild_caches()
        assert_equal(len(self.ag.residues), 2)

    # INVALID: no caches in AtomGroups anymore
    @skip
    def test_rebuild_cache_segments(self):
        # This test is similar to above, but a second segment has to be taken
        # from a new universe
        assert_equal(len(self.ag.segments), 1)

        u2 = MDAnalysis.Universe(PSF_notop, DCD)
        self.ag[0].segment = u2.atoms[0].segment

        assert_equal(len(self.ag.segments), 1)
        self.ag._rebuild_caches()
        assert_equal(len(self.ag.segments), 2)

    # INVALID: no caches in AtomGroups anymore
    @skip
    def test_atom_cachesize_change(self):
        # By default 10,000 atoms are required to necessitate cache lookup, we
        # can change this though
        ag = self.u.atoms[:100]
        # run a __contains__ query
        self.u.atoms[10] in ag
        # Check that cache wasn't used
        assert_equal('atoms' in ag._cache, False)
        ag._atomcache_size = 50  # now will make cache if size > 50
        # Run another query
        self.u.atoms[10] in ag
        # Check if cache was built this time
        assert_equal('atoms' in ag._cache, True)

    # INVALID: no caches in AtomGroups anymore
    @skip
    def test_atomcache_use(self):
        # Tests that lookup with 'atoms' cache works
        ag = self.u.atoms[:100]
        ag._atomcache_size = 50
        assert_equal(self.u.atoms[50] in ag, True)

    # INVALID: no caches in AtomGroups anymore
    @skip
    def test_rebuild_atomcache_no(self):
        # Don't always add atoms into cache
        ag = self.u.atoms[:100]
        ag._rebuild_caches()
        assert_equal('atoms' in ag._cache, False)

    # INVALID: no caches in AtomGroups anymore
    @skip
    def test_rebuild_atomcache(self):
        # Tests that 'atoms' is built into cache if size is enough
        ag = self.u.atoms[:100]
        ag._atomcache_size = 50
        ag._rebuild_caches()
        assert_equal('atoms' in ag._cache, True)

    # INVALID: can only set dimensions from the Universe
    @skip
    def test_set_dimensions(self):
        u = MDAnalysis.Universe(PSF, DCD)
        box = np.array([10, 11, 12, 90, 90, 90])
        u.atoms.dimensions = np.array([10, 11, 12, 90, 90, 90])
        assert_allclose(u.dimensions, box)
        assert_allclose(u.atoms.dimensions, box)


class TestUniverseSetTopology(TestCase):
    """Tests setting of bonds/angles/dihedrals/impropers from Universe."""

    # VALID
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.u = MDAnalysis.Universe(PSF, DCD)

    # VALID
    def tearDown(self):
        del self.u

    # INVALID: universe doesn't have bonds; AtomGroups do, though
    @skip
    def test_set_bonds(self):
        assert_equal(len(self.u.bonds), 3365)
        assert_equal(len(self.u.atoms[0].bonds), 4)

        self.u.bonds = []

        assert_equal(len(self.u.bonds), 0)
        assert_equal(len(self.u.atoms[0].bonds), 0)

    # INVALID: universe doesn't have angles; AtomGroups do, though
    @skip
    def test_set_angles(self):
        assert_equal(len(self.u.angles), 6123)
        assert_equal(len(self.u.atoms[0].angles), 9)

        self.u.angles = []

        assert_equal(len(self.u.angles), 0)
        assert_equal(len(self.u.atoms[0].angles), 0)

    # INVALID: universe doesn't have dihedrals; AtomGroups do, though
    @skip
    def test_set_dihedrals(self):
        assert_equal(len(self.u.dihedrals), 8921)
        assert_equal(len(self.u.atoms[0].dihedrals), 14)

        self.u.dihedrals = []

        assert_equal(len(self.u.dihedrals), 0)
        assert_equal(len(self.u.atoms[0].dihedrals), 0)

    # INVALID: universe doesn't have impropers; AtomGroups do, though
    @skip
    def test_set_impropers(self):
        assert_equal(len(self.u.impropers), 541)
        assert_equal(len(self.u.atoms[4].impropers), 1)

        self.u.impropers = []

        assert_equal(len(self.u.impropers), 0)
        assert_equal(len(self.u.atoms[4].impropers), 0)

    # Test deleting topology information
    # In general, access it to make sure it's built
    # Assert it's in cache
    # Delete
    # Assert it's not in cache

    # INVALID: universe has no direct bonds access
    @skip
    def test_bonds_delete(self):
        self.u.bonds
        self.u.atoms[0].bonds

        assert_equal('bonds' in self.u._cache, True)
        assert_equal('bondDict' in self.u._cache, True)

        del self.u.bonds

        assert_equal('bonds' in self.u._cache, False)
        assert_equal('bondDict' in self.u._cache, False)

    # INVALID: universe has no direct angles access
    @skip
    def test_angles_delete(self):
        self.u.angles
        self.u.atoms[0].angles

        assert_equal('angles' in self.u._cache, True)
        assert_equal('angleDict' in self.u._cache, True)

        del self.u.angles

        assert_equal('angles' in self.u._cache, False)
        assert_equal('angleDict' in self.u._cache, False)

    # INVALID: universe has no direct dihedrals access
    @skip
    def test_dihedrals_delete(self):
        self.u.dihedrals
        self.u.atoms[0].dihedrals

        assert_equal('dihedrals' in self.u._cache, True)
        assert_equal('dihedralDict' in self.u._cache, True)

        del self.u.dihedrals

        assert_equal('dihedrals' in self.u._cache, False)
        assert_equal('dihedralDict' in self.u._cache, False)

    # INVALID: universe has no direct impropers access
    @skip
    def test_impropers_delete(self):
        self.u.impropers
        self.u.atoms[0].impropers

        assert_equal('impropers' in self.u._cache, True)
        assert_equal('improperDict' in self.u._cache, True)

        del self.u.impropers

        assert_equal('impropers' in self.u._cache, False)
        assert_equal('improperDict' in self.u._cache, False)


class TestResidue(TestCase):
    # VALID
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.res = self.universe.residues[100]

    # INVALID: a Universe has its own Residue class, so comparing
    # to base class definition isn't meaningful
    @skip
    def test_type(self):
        assert_equal(type(self.res), MDAnalysis.core.groups.Residue)
        assert_equal(self.res.name, "ILE")
        assert_equal(self.res.id, 101)

    # INVALID: residues do not behave like AtomGroups anymore,
    # so cannot index in this way; should use `res.atoms[2]`
    @skip
    def test_index(self):
        atom = self.res[2]
        assert_equal(type(atom), MDAnalysis.core.groups.Atom)
        assert_equal(atom.name, "CA")
        assert_equal(atom.index, 1522)
        assert_equal(atom.resid, 101)

    # INVALID: residues do not behave like AtomGroups anymore,
    # so cannot slice in this way
    @skip
    def test_slicing(self):
        atoms = self.res[2:10:2]
        assert_equal(len(atoms), 4)
        assert_equal(type(atoms), MDAnalysis.core.groups.AtomGroup)

    # INVALID: residues do not behave like AtomGroups anymore,
    # so cannot slice in this way
    @skip
    def test_advanced_slicing(self):
        atoms = self.res[[0, 2, -2, -1]]
        assert_equal(len(atoms), 4)
        assert_equal(type(atoms), MDAnalysis.core.groups.AtomGroup)
        assert_equal(atoms.names, ["N", "CA", "C", "O"])

    # VALID
    def test_atom_order(self):
        assert_equal(self.res.atoms.indices,
                     sorted(self.res.atoms.indices))


class TestResidueGroup(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    # VALID
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.rg = self.universe.atoms.residues

    # VALID
    def test_newResidueGroup(self):
        """test that slicing a ResidueGroup returns a new ResidueGroup
        (Issue 135)"""
        rg = self.universe.atoms.residues
        newrg = rg[10:20:2]
        assert_equal(type(newrg), type(rg),
                     "Failed to make a new ResidueGroup: type mismatch")
        assert_equal(len(newrg), len(rg[10:20:2]))

    # VALID
    def test_n_atoms(self):
        assert_equal(self.rg.n_atoms, 3341)

    # VALID
    def test_n_residues(self):
        assert_equal(self.rg.n_residues, 214)

    # VALID
    def test_resids_dim(self):
        assert_equal(len(self.rg.resids), len(self.rg))

    # INVALID: this topology has no resnums, so no resnums property
    @skip
    def test_resnums_dim(self):
        assert_equal(len(self.rg.resnums), len(self.rg))

    # VALID
    def test_segids_dim(self):
        assert_equal(len(self.rg.segids), len(self.rg))

    # VALID
    def test_len(self):
        """testing that len(residuegroup) == residuegroup.n_residues"""
        assert_equal(len(self.rg), self.rg.n_residues,
                     "len and n_residues disagree")

    # INVALID: set resids with `ResidueGroup.resids` property; no `set_resids` method
    @skip
    def test_set_resids(self):
        rg = self.universe.select_atoms("bynum 12:42").residues
        resid = 999
        rg.set_resids(resid)
        # check individual atoms
        assert_equal([a.resid for a in rg.atoms],
                     resid * np.ones(rg.n_atoms),
                     err_msg="failed to set_resid atoms 12:42 to same resid")
        # check residues
        assert_equal(rg.resids, resid * np.ones(rg.n_residues),
                     err_msg="failed to set_resid of residues belonging to "
                     "atoms 12:42 to same resid")

    # INVALID: set resids with `ResidueGroup.resids` property; no `set_resids` method
    @skip
    def test_set_resids(self):
        """test_set_resid: set ResidueGroup resids on a per-residue basis"""
        rg = self.universe.select_atoms("resid 10:18").residues
        resids = np.array(rg.resids) + 1000
        rg.set_resids(resids)
        # check individual atoms
        for r, resid in zip(rg, resids):
            assert_equal([a.resid for a in r.atoms],
                         resid * np.ones(r.n_atoms),
                         err_msg="failed to set_resid residues 10:18 to same "
                         "resid in residue {0}\n"
                         "(resids = {1}\nresidues = {2})".format(r, resids, rg))
        # check residues
        # NOTE: need to create a new selection because underlying Residue
        #       objects are not changed; only Atoms are changed, and Residues
        #       are rebuilt from Atoms.
        rgnew = self.universe.select_atoms("resid 1010:1018").residues
        assert_equal(rgnew.resids, np.unique(resids),
                     err_msg="failed to set_resid of residues belonging to "
                     "residues 10:18 to new resids")

    # INVALID: set resids with `ResidueGroup.resids` property; no `set_resids` method
    @skip
    def test_set_resids_updates_self(self):
        rg = self.universe.select_atoms("resid 10:18").residues
        resids = np.array(rg.resids) + 1000
        rg.set_resids(resids)
        assert_equal(rg.resids, np.unique(resids),
                     err_msg="old selection was not changed in place "
                     "after set_resid")

    # INVALID: no resnums in this topology, so no resnums property
    @skip
    def test_set_resnum_single(self):
        rg = self.universe.residues[:3]
        new = 22
        rg.set_resnums(new)

        assert_equal(all(rg.resnums == new), True)
        for r in rg:
            assert_equal(r.resnum, new)

    # INVALID: no resnums in this topology, so no resnums property
    @skip
    def test_set_resnum_many(self):
        rg = self.universe.residues[:3]
        new = [22, 23, 24]
        rg.set_resnums(new)

        assert_equal(all(rg.resnums == new), True)
        for r, v in zip(rg, new):
            assert_equal(r.resnum, v)

    # INVALID: no resnums in this topology, so no resnums property
    @skip
    def test_set_resnum_ValueError(self):
        rg = self.universe.residues[:3]
        new = [22, 23, 24, 25]

        assert_raises(ValueError, rg.set_resnums, new)

    # INVALID: no `set_resnames` method; use `resnames` property directly
    @skip
    def test_set_resname_single(self):
        rg = self.universe.residues[:3]
        new = 'newname'

        rg.set_resnames(new)
        assert_equal(all(rg.resnames == new), True)
        for r in rg:
            assert_equal(r.name, new)

    # INVALID: no `set_resnames` method; use `resnames` property directly
    @skip
    def test_set_resname_many(self):
        rg = self.universe.residues[:3]
        new = ['a', 'b', 'c']
        rg.set_resnames(new)

        assert_equal(all(rg.resnames == new), True)
        for r, v in zip(rg, new):
            assert_equal(r.name, v)

    # INVALID: no `set_resnames` method; use `resnames` property directly
    @skip
    def test_set_resname_ValueError(self):
        rg = self.universe.residues[:3]
        new = ['a', 'b', 'c', 'd']

        assert_raises(ValueError, rg.set_resnames, new)

    # INVALID: no `set_resids` method; also, residues are not mergeable
    # by setting resids; resids are not necessarily unique; atoms must
    # have their resindex set to change residue membership
    @skip
    def test_merge_residues(self):
        rg = self.universe.select_atoms("resid 12:14").residues
        nres_old = self.universe.atoms.n_residues
        natoms_old = rg.n_atoms
        rg.set_resids(12)  # merge all into one with resid 12
        nres_new = self.universe.atoms.n_residues
        r_merged = self.universe.select_atoms("resid 12:14").residues
        natoms_new = self.universe.select_atoms("resid 12").n_atoms
        assert_equal(len(r_merged), 1, err_msg="set_resid failed to merge "
                     "residues: merged = {0}".format(r_merged))
        assert_equal(nres_new, nres_old - 2,
                     err_msg="set_resid failed to merge residues: "
                     "merged = {0}".format(r_merged))
        assert_equal(natoms_new, natoms_old, err_msg="set_resid lost atoms "
                     "on merge".format(r_merged))

        assert_equal(self.universe.residues.n_residues,
                     self.universe.atoms.n_residues,
                     err_msg="Universe.residues and Universe.atoms.n_residues "
                     "do not agree after residue "
                             "merge.")

    # INVALID: no `set_masses` method; use `masses` property directly
    @skip
    def test_set_masses(self):
        rg = self.universe.select_atoms("bynum 12:42 and name H*").residues
        mass = 2.0
        rg.set_masses(mass)
        # check individual atoms
        assert_equal([a.mass for a in rg.atoms],
                     mass * np.ones(rg.n_atoms),
                     err_msg="failed to set_mass H* atoms in resid 12:42 to {0}".format(mass))

    # VALID
    def test_atom_order(self):
        assert_equal(self.universe.residues.atoms.indices,
                     sorted(self.universe.residues.atoms.indices))


class TestSegment(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')

    # INVALID: no `set_segids` method for ResidueGroup; c
    @skip
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.universe.residues[:100].set_segids("A")  # make up some segments
        self.universe.residues[100:150].set_segids("B")
        self.universe.residues[150:].set_segids("C")
        self.sB = self.universe.segments[1]

    # VALID but temporary
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)

    # INVALID: `Segment` from a particular Universe will not be the same class
    # in core; each Universe generates its own set of classes based on topology
    @skip
    def test_type(self):
        assert_equal(type(self.sB), MDAnalysis.core.groups.Segment)
        assert_equal(self.sB.name, "B")

    # INVALID: `Segment` doesn't behave as a `ResidueGroup`;
    # no index behavior
    @skip
    def test_index(self):
        s = self.sB
        res = s[5]
        assert_equal(type(res), MDAnalysis.core.groups.Residue)

    # INVALID: `Segment` doesn't behave as a `ResidueGroup`;
    # no slicing behavior
    @skip
    def test_slicing(self):
        res = self.sB[5:10]
        assert_equal(len(res), 5)
        assert_equal(type(res), MDAnalysis.core.groups.ResidueGroup)

    # INVALID: `Segment` doesn't behave as a `ResidueGroup`;
    # no slicing behavior
    @skip
    def test_advanced_slicing(self):
        res = self.sB[[3, 7, 2, 4]]
        assert_equal(len(res), 4)
        assert_equal(type(res), MDAnalysis.core.groups.ResidueGroup)

    # INVALID: no `name` or `id` attribute for `Segment`; only `segid`
    @skip
    def test_id(self):
        assert_equal(self.sB.name, self.sB.id)

    # INVALID: no `name` or `id` attribute for `Segment`; only `segid`
    @skip
    def test_set_id(self):
        # Test setting the name via the id attribute
        new = 'something'
        self.sB.id = new
        for val in [self.sB.id, self.sB.name]:
            assert_equal(val, new)

    # VALID
    def test_atom_order(self):
        assert_equal(self.universe.segments[0].atoms.indices,
                     sorted(self.universe.segments[0].atoms.indices))


class TestSegmentGroup(TestCase):
    # VALID
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.g = self.universe.atoms.segments

    # VALID
    def test_newSegmentGroup(self):
        """test that slicing a SegmentGroup returns a new SegmentGroup (Issue 135)"""
        g = self.universe.atoms.segments
        newg = g[:]
        assert_equal(type(newg), type(g), "Failed to make a new SegmentGroup: type mismatch")
        assert_equal(len(newg), len(g))

    # VALID
    def test_n_atoms(self):
        assert_equal(self.g.n_atoms, 3341)

    # VALID
    def test_n_residues(self):
        assert_equal(self.g.n_residues, 214)

    # INVALID: `SegmentGroup.resids` gives list of arrays, one array for each segment
    @skip
    def test_resids_dim(self):
        assert_equal(len(self.g.resids), len(self.g.residues))

    # INVALID: topology has no resnums
    @skip
    def test_resnums_dim(self):
        assert_equal(len(self.g.resnums), len(self.g.residues))

    # VALID
    def test_segids_dim(self):
        assert_equal(len(self.g.segids), len(self.g))

    # INVALID: cannot set resids from `SegmentGroup`; no `set_resids` method
    @skip
    def test_set_resids(self):
        g = self.universe.select_atoms("bynum 12:42").segments
        resid = 999
        g.set_resids(resid)
        # check individual atoms
        assert_equal([a.resid for a in g.atoms],
                     resid * np.ones(g.n_atoms),
                     err_msg="failed to set_resid for segment to same resid")
        # check residues
        assert_equal(g.residues.resids, resid * np.ones(g.n_residues),
                     err_msg="failed to set_resid of segments belonging to atoms 12:42 to same resid")

    # INVALID: cannot set resids from `SegmentGroup`; no `set_resids` method
    @skip
    def test_set_resids(self):
        g = self.universe.select_atoms("resid 10:18").segments
        resid = 999
        g.set_resids(resid * np.ones(len(g)))
        # note: all is now one residue... not meaningful but it is the correct behaviour
        assert_equal(g.resids, [resid],
                     err_msg="failed to set_resid  in Segment {0}".format(g))

    # INVALID: no `set_segids` method; use `segids` property directly
    @skip
    def test_set_segids(self):
        s = self.universe.select_atoms('all').segments
        s.set_segids(['ADK'])
        assert_equal(self.universe.segments.segids, ['ADK'],
                     err_msg="failed to set_segid on segments")

    # INVALID: no `set_segids` method; use `segids` property directly
    @skip
    def test_set_segid_updates_self(self):
        g = self.universe.select_atoms("resid 10:18").segments
        g.set_segids('ADK')
        assert_equal(g.segids, ['ADK'],
                     err_msg="old selection was not changed in place after set_segid")

    # INVALID: no `set_segids` method; use `segids` property directly
    @skip
    def test_set_masses(self):
        g = self.universe.select_atoms("bynum 12:42 and name H*").segments
        mass = 2.0
        g.set_masses(mass)
        # check individual atoms
        assert_equal([a.mass for a in g.atoms],
                     mass * np.ones(g.n_atoms),
                     err_msg="failed to set_mass in segment of  H* atoms in resid 12:42 to {0}".format(mass))

    # INVALID: no `set_segids` method; use `segids` property directly
    @skip
    def test_set_segid_ValueError(self):
        assert_raises(ValueError, self.g.set_resids, [1, 2, 3, 4])

    # VALID
    def test_atom_order(self):
        assert_equal(self.universe.segments.atoms.indices,
                     sorted(self.universe.segments.atoms.indices))


class TestAtomGroupVelocities(TestCase):
    """Tests of velocity-related functions in AtomGroup"""

    # VALID
    def setUp(self):
        self.universe = MDAnalysis.Universe(GRO, TRR)
        self.ag = self.universe.select_atoms("bynum 12:42")

    # INVALID: no `get_velocities` method; use `velocities` property directly
    @skip
    @dec.slow
    def test_get_velocities(self):
        v = self.ag.get_velocities()
        assert_(np.any(np.abs(v) > 1e-6), "velocities should be non-zero")

    # VALID
    @dec.slow
    def test_velocities(self):
        ag = self.universe.atoms[42:45]
        ref_v = np.array([
            [-3.61757946, -4.9867239, 2.46281552],
            [2.57792854, 3.25411797, -0.75065529],
            [13.91627216, 30.17778587, -12.16669178]])
        v = ag.velocities
        assert_almost_equal(v, ref_v, err_msg="velocities were not read correctly")

    # INVALID: no `get_velocities` method; use `velocities` property directly
    @skip
    @dec.slow
    def test_set_velocities(self):
        ag = self.ag
        v = ag.get_velocities() - 2.7271
        ag.set_velocities(v)
        assert_almost_equal(ag.get_velocities(), v,
                            err_msg="messages were not set to new value")


class TestAtomGroupForces(TestCase):
    """Tests of velocity-related functions in AtomGroup"""

    # VALID
    def setUp(self):
        self.universe = MDAnalysis.Universe(COORDINATES_XYZ, COORDINATES_TRR)
        self.ag = self.universe.select_atoms("bynum 12:42")

    # INVALID: no `get_forces` method; use `forces` property directly
    @skip
    @dec.slow
    def test_get_forces(self):
        v = self.ag.get_forces()
        assert_(np.any(np.abs(v) > 1e-6), "forces should be non-zero")

    # VALID
    @dec.slow
    def test_forces(self):
        ag = self.universe.atoms[1:4]
        ref_v = np.arange(9).reshape(3, 3) * .01 + .03
        v = ag.forces
        assert_almost_equal(v, ref_v, err_msg="forces were not read correctly")

    # INVALID: no `get_forces` method; use `forces` property directly
    @skip
    @dec.slow
    def test_set_forces(self):
        ag = self.ag
        v = ag.get_forces() - 2.7271
        ag.set_forces(v)
        assert_almost_equal(ag.get_forces(), v,
                            err_msg="messages were not set to new value")


class TestAtomGroupTimestep(TestCase):
    """Tests the AtomGroup.ts attribute (partial timestep)"""

    # VALID
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = MDAnalysis.Universe(TRZ_psf, TRZ)
        self.prec = 6

    # VALID
    def tearDown(self):
        del self.universe
        del self.prec


    # VALID: but should be testing non-hidden attributes
    def test_partial_timestep(self):
        ag = self.universe.select_atoms('name Ca')
        idx = ag.indices

        assert_equal(len(ag.ts._pos), len(ag))

        for ts in self.universe.trajectory[0:20:5]:
            assert_array_almost_equal(ts._pos[idx], ag.ts._pos, self.prec,
                                      err_msg="Partial timestep coordinates wrong")
            assert_array_almost_equal(ts._velocities[idx], ag.ts._velocities, self.prec,
                                      err_msg="Partial timestep coordinates wrong")


# INVALID: AtomGroups can't exist without a Universe
@skip
def test_empty_AtomGroup():
    """Test that an empty AtomGroup can be constructed (Issue 12)"""
    ag = MDAnalysis.core.groups.AtomGroup([])
    assert_equal(len(ag), 0)


class _WriteAtoms(TestCase):
    """Set up the standard AdK system in implicit solvent."""
    ext = None  # override to test various output writers
    precision = 3

    # VALID
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        suffix = '.' + self.ext
        self.tempdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tempdir.name, 'writeatoms' + suffix)

    # VALID
    def tearDown(self):
        del self.universe
        del self.tempdir

    # VALID
    def universe_from_tmp(self):
        return MDAnalysis.Universe(self.outfile, convert_units=True)

    # VALID
    def test_write_atoms(self):
        self.universe.atoms.write(self.outfile)
        u2 = self.universe_from_tmp()
        assert_array_almost_equal(self.universe.atoms.positions, u2.atoms.positions, self.precision,
                                  err_msg="atom coordinate mismatch between original and {0!s} file".format(self.ext))

    # VALID
    def test_write_empty_atomgroup(self):
        sel = self.universe.select_atoms('name doesntexist')
        assert_raises(IndexError, sel.write, self.outfile)

    # VALID
    def test_write_selection(self):
        CA = self.universe.select_atoms('name CA')
        CA.write(self.outfile)
        u2 = self.universe_from_tmp()
        CA2 = u2.select_atoms('all')  # check EVERYTHING, otherwise we might get false positives!
        assert_equal(len(u2.atoms), len(CA.atoms), "written CA selection does not match original selection")
        assert_almost_equal(CA2.positions, CA.positions, self.precision,
                            err_msg="CA coordinates do not agree with original")

    # INVALID: Only `AtomGroup`s have `write` method. Must do `G.atoms.write`
    @skip
    def test_write_Residue(self):
        G = self.universe.s4AKE.ARG[-2]  # 2nd but last Arg
        G.write(self.outfile)
        u2 = self.universe_from_tmp()
        G2 = u2.select_atoms('all')  # check EVERYTHING, otherwise we might get false positives!
        assert_equal(len(u2.atoms), len(G.atoms), "written R206 Residue does not match original ResidueGroup")
        assert_almost_equal(G2.positions, G.positions, self.precision,
                            err_msg="Residue R206 coordinates do not agree with original")

    # INVALID: Only `AtomGroup`s have `write` method. Must do `G.atoms.write`
    @skip
    def test_write_ResidueGroup(self):
        G = self.universe.s4AKE.LEU
        G.write(self.outfile)
        u2 = self.universe_from_tmp()
        G2 = u2.select_atoms('all')  # check EVERYTHING, otherwise we might get false positives!
        assert_equal(len(u2.atoms), len(G.atoms), "written LEU ResidueGroup does not match original ResidueGroup")
        assert_almost_equal(G2.positions, G.positions, self.precision,
                            err_msg="ResidueGroup LEU coordinates do not agree with original")

    # INVALID: Only `AtomGroup`s have `write` method. Must do `G.atoms.write`
    @skip
    def test_write_Segment(self):
        G = self.universe.s4AKE
        G.write(self.outfile)
        u2 = self.universe_from_tmp()
        G2 = u2.select_atoms('all')  # check EVERYTHING, otherwise we might get false positives!
        assert_equal(len(u2.atoms), len(G.atoms), "written s4AKE segment does not match original segment")
        assert_almost_equal(G2.positions, G.positions, self.precision,
                            err_msg="segment s4AKE coordinates do not agree with original")

    # VALID
    def test_write_Universe(self):
        U = self.universe
        W = MDAnalysis.Writer(self.outfile)
        W.write(U)
        W.close()
        u2 = self.universe_from_tmp()
        assert_equal(len(u2.atoms), len(U.atoms), "written 4AKE universe does not match original universe in size")
        assert_almost_equal(u2.atoms.positions, U.atoms.positions, self.precision,
                            err_msg="written universe 4AKE coordinates do not agree with original")


class TestWritePDB(_WriteAtoms):
    ext = "pdb"
    precision = 3


import MDAnalysis.coordinates


class TestWriteGRO(_WriteAtoms):
    ext = "gro"
    precision = 2

    # INVALID: flags should be retired(?)
    @skip
    def test_flag_convert_length(self):
        assert_equal(MDAnalysis.core.flags['convert_lengths'], True,
                     "The flag convert_lengths SHOULD be True by default! "
                     "(If it is not then this might indicate a race condition in the "
                     "testing suite.)")


# VALID
@attr("issue")
@dec.skipif(parser_not_found('DCD'),
            'DCD parser not available. Are you using python 3?')
def test_generated_residueselection():
    """Test that a generated residue group always returns a ResidueGroup (Issue 47)
    unless there is a single residue (Issue 363 change)"""
    universe = MDAnalysis.Universe(PSF, DCD)
    # only a single Cys in AdK
    cys = universe.s4AKE.CYS
    assert_(isinstance(cys, MDAnalysis.core.groups.Residue),
            "Single Cys77 is NOT returned as a single Residue (Issue 47)")

    # multiple Met
    met = universe.s4AKE.MET
    assert_(isinstance(met, MDAnalysis.core.groups.ResidueGroup),
            "Met selection does not return a ResidueGroup")

    del universe


# VALID
@attr('issue')
@dec.skipif(parser_not_found('DCD'),
            'DCD parser not available. Are you using python 3?')
def test_instantselection_termini():
    """Test that instant selections work, even for residues that are also termini (Issue 70)"""
    universe = MDAnalysis.Universe(PSF, DCD)
    assert_equal(universe.residues[20].CA.name, 'CA', "CA of MET21 is not selected correctly")
    del universe

class TestPBCFlag(object):
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def setUp(self):
        self.prec = 3
        self.universe = MDAnalysis.Universe(TRZ_psf, TRZ)
        self.ref_noPBC = {
            'COG': np.array([4.23789883, 0.62429816, 2.43123484], dtype=np.float32),
            'COM': np.array([4.1673783, 0.70507009, 2.21175832]),
            'ROG': 119.30368949900134, 'Shape': 0.6690026954813445,
            'Asph': 0.5305456387833748,
            'MOI': np.array([
                [152117.06620921, 55149.54042136, -26630.46034023],
                [55149.54042136, 72869.64061494, 21998.1778074],
                [-26630.46034023, 21998.1778074, 162388.70002471]]),
            'BBox': np.array([[-75.74159241, -144.86634827, -94.47974396], [95.83090973, 115.11561584, 88.09812927]],
                                dtype=np.float32),
            'BSph': (173.40482, np.array([4.23789883, 0.62429816, 2.43123484], dtype=np.float32)),
            'PAxes': np.array([
                [0.46294889, -0.85135849, 0.24671249],
                [0.40611024, 0.45112859, 0.7947059],
                [-0.78787867, -0.26771575, 0.55459488]])
        }
        self.ref_PBC = {
            'COG': np.array([26.82960892, 31.5592289, 30.98238945], dtype=np.float32),
            'COM': np.array([26.67781143, 31.2104336, 31.19796289]),
            'ROG': 27.713008969174918, 'Shape': 0.0017390512580463542,
            'Asph': 0.020601215358731016,
            'MOI': np.array([
                [7333.79167791, -211.8997285, -721.50785456],
                [-211.8997285, 7059.07470427, -91.32156884],
                [-721.50785456, -91.32156884, 6509.31735029]]),
            'BBox': np.array(
                [[1.45964116e-01, 1.85623169e-02, 4.31785583e-02], [5.53314018e+01, 5.54227829e+01, 5.54158211e+01]],
                dtype=np.float32),
            'BSph': (47.923367, np.array([26.82960892, 31.5592289, 30.98238945], dtype=np.float32)),
            'PAxes': np.array([
                [-0.50622389, -0.18364489, -0.84262206],
                [-0.07520116, -0.96394227, 0.25526473],
                [-0.85911708, 0.19258726, 0.4741603]])
        }
        self.ag = self.universe.residues[0:3]

    def tearDown(self):
        MDAnalysis.core.flags['use_pbc'] = False
        del self.universe
        del self.ref_noPBC
        del self.ref_PBC
        del self.ag

    def test_flag(self):
        # Test default setting of flag
        assert_equal(MDAnalysis.core.flags['use_pbc'], False)

    def test_default(self):
        # Test regular behaviour
        assert_almost_equal(self.ag.center_of_geometry(), self.ref_noPBC['COG'], self.prec)
        assert_almost_equal(self.ag.center_of_mass(), self.ref_noPBC['COM'], self.prec)
        assert_almost_equal(self.ag.radius_of_gyration(), self.ref_noPBC['ROG'], self.prec)
        assert_almost_equal(self.ag.shape_parameter(), self.ref_noPBC['Shape'], self.prec)
        assert_almost_equal(self.ag.asphericity(), self.ref_noPBC['Asph'], self.prec)
        assert_almost_equal(self.ag.moment_of_inertia(), self.ref_noPBC['MOI'], self.prec)
        assert_almost_equal(self.ag.bbox(), self.ref_noPBC['BBox'], self.prec)
        assert_almost_equal(self.ag.bsphere()[0], self.ref_noPBC['BSph'][0], self.prec)
        assert_almost_equal(self.ag.bsphere()[1], self.ref_noPBC['BSph'][1], self.prec)
        assert_almost_equal(self.ag.principal_axes(), self.ref_noPBC['PAxes'], self.prec)

    def test_pbcflag(self):
        # Test using ag method flag
        assert_almost_equal(self.ag.center_of_geometry(pbc=True), self.ref_PBC['COG'], self.prec)
        assert_almost_equal(self.ag.center_of_mass(pbc=True), self.ref_PBC['COM'], self.prec)
        assert_almost_equal(self.ag.radius_of_gyration(pbc=True), self.ref_PBC['ROG'], self.prec)
        assert_almost_equal(self.ag.shape_parameter(pbc=True), self.ref_PBC['Shape'], self.prec)
        assert_almost_equal(self.ag.asphericity(pbc=True), self.ref_PBC['Asph'], self.prec)
        assert_almost_equal(self.ag.moment_of_inertia(pbc=True), self.ref_PBC['MOI'], self.prec)
        assert_almost_equal(self.ag.bbox(pbc=True), self.ref_PBC['BBox'], self.prec)
        assert_almost_equal(self.ag.bsphere(pbc=True)[0], self.ref_PBC['BSph'][0], self.prec)
        assert_almost_equal(self.ag.bsphere(pbc=True)[1], self.ref_PBC['BSph'][1], self.prec)
        assert_almost_equal(self.ag.principal_axes(pbc=True), self.ref_PBC['PAxes'], self.prec)

    def test_usepbc_flag(self):
        # Test using the core.flags flag
        MDAnalysis.core.flags['use_pbc'] = True
        assert_almost_equal(self.ag.center_of_geometry(), self.ref_PBC['COG'], self.prec)
        assert_almost_equal(self.ag.center_of_mass(), self.ref_PBC['COM'], self.prec)
        assert_almost_equal(self.ag.radius_of_gyration(), self.ref_PBC['ROG'], self.prec)
        assert_almost_equal(self.ag.shape_parameter(), self.ref_PBC['Shape'], self.prec)
        assert_almost_equal(self.ag.asphericity(), self.ref_PBC['Asph'], self.prec)
        assert_almost_equal(self.ag.moment_of_inertia(), self.ref_PBC['MOI'], self.prec)
        assert_almost_equal(self.ag.bbox(), self.ref_PBC['BBox'], self.prec)
        assert_almost_equal(self.ag.bsphere()[0], self.ref_PBC['BSph'][0], self.prec)
        assert_almost_equal(self.ag.bsphere()[1], self.ref_PBC['BSph'][1], self.prec)
        assert_almost_equal(self.ag.principal_axes(), self.ref_PBC['PAxes'], self.prec)
        MDAnalysis.core.flags['use_pbc'] = False

    def test_override_flag(self):
        # Test using the core.flags flag, then overriding
        MDAnalysis.core.flags['use_pbc'] = True
        assert_almost_equal(self.ag.center_of_geometry(pbc=False), self.ref_noPBC['COG'], self.prec)
        assert_almost_equal(self.ag.center_of_mass(pbc=False), self.ref_noPBC['COM'], self.prec)
        assert_almost_equal(self.ag.radius_of_gyration(pbc=False), self.ref_noPBC['ROG'], self.prec)
        assert_almost_equal(self.ag.shape_parameter(pbc=False), self.ref_noPBC['Shape'], self.prec)
        assert_almost_equal(self.ag.asphericity(pbc=False), self.ref_noPBC['Asph'], self.prec)
        assert_almost_equal(self.ag.moment_of_inertia(pbc=False), self.ref_noPBC['MOI'], self.prec)
        assert_almost_equal(self.ag.bbox(pbc=False), self.ref_noPBC['BBox'], self.prec)
        assert_almost_equal(self.ag.bsphere(pbc=False)[0], self.ref_noPBC['BSph'][0], self.prec)
        assert_almost_equal(self.ag.bsphere(pbc=False)[1], self.ref_noPBC['BSph'][1], self.prec)
        assert_almost_equal(self.ag.principal_axes(pbc=False), self.ref_noPBC['PAxes'], self.prec)
        MDAnalysis.core.flags['use_pbc'] = False


# INVALID: not including as_Universe, since not clear what it does that's different from `Universe()`
@skip
class TestAsUniverse(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.u = MDAnalysis.Universe(PSF_notop, DCD)

    def tearDown(self):
        del self.u

    def test_empty_TypeError(self):
        assert_raises(TypeError, as_Universe)

    def test_passback(self):
        returnval = as_Universe(self.u)

        assert_equal(returnval is self.u, True)

    def test_makeuni(self):
        returnval = as_Universe(PSF_notop, DCD)

        ## __eq__ method for Universe doesn't exist, make one up here
        assert_equal(set(returnval.atoms), set(self.u.atoms))


class TestUniverseCache(TestCase):
    # VALID: we should support building Universes with empty Topologies somehow
    def setUp(self):
        self.u = MDAnalysis.Universe()  # not using atoms so just blank universe
        self.fill = [1, 2, 3]

    # VALID
    def tearDown(self):
        del self.u
        del self.fill

    # INVALID: no caches at all in Universe
    @skip
    def test_add_to_cache(self):
        # add an item to cache and see if it sticks
        cache = 'aa'
        self.u._fill_cache(cache, self.fill)

        assert_equal('aa' in self.u._cache, True)
        assert_equal(self.u._cache[cache], self.fill)

    # INVALID: no caches at all in Universe
    @skip
    def test_remove_single(self):
        # remove a single item from cache
        cache = 'bb'

        self.u._fill_cache(cache, self.fill)

        assert_equal(cache in self.u._cache, True)

        self.u._clear_caches(cache)

        assert_equal(cache in self.u._cache, False)

    # INVALID: no caches at all in Universe
    @skip
    def test_remove_list(self):
        # remove a few things from cache
        caches = ['cc', 'dd']
        for c in caches:
            self.u._fill_cache(c, self.fill)

        for c in caches:
            assert_equal(c in self.u._cache, True)

        self.u._clear_caches(*caches)

        for c in caches:
            assert_equal(c in self.u._cache, False)

    # INVALID: no caches at all in Universe
    @skip
    def test_clear_all(self):
        # remove everything from cache
        caches = ['ee', 'ff', 'gg']
        for c in caches:
            self.u._fill_cache(c, self.fill)

        self.u._clear_caches()

        assert_equal(self.u._cache, dict())


class TestCustomReaders(TestCase):
    """
    Can pass a reader as kwarg on Universe creation
    """
    # VALID
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def test_custom_reader(self):
        # check that reader passing works
        u = MDAnalysis.Universe(TRZ_psf, TRZ, format=MDAnalysis.coordinates.TRZ.TRZReader)
        assert_equal(len(u.atoms), 8184)

    # VALID
    def test_custom_reader_singleframe(self):
        T = MDAnalysis.topology.GROParser.GROParser
        R = MDAnalysis.coordinates.GRO.GROReader
        u = MDAnalysis.Universe(two_water_gro, two_water_gro,
                                topology_format=T, format=R)
        assert_equal(len(u.atoms), 6)

    # VALID
    def test_custom_reader_singleframe_2(self):
        # Same as before, but only one argument to Universe
        T = MDAnalysis.topology.GROParser.GROParser
        R = MDAnalysis.coordinates.GRO.GROReader
        u = MDAnalysis.Universe(two_water_gro,
                                topology_format=T, format=R)
        assert_equal(len(u.atoms), 6)

    # VALID
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def test_custom_parser(self):
        # topology reader passing works
        u = MDAnalysis.Universe(TRZ_psf, TRZ, topology_format=MDAnalysis.topology.PSFParser.PSFParser)
        assert_equal(len(u.atoms), 8184)

    # VALID
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def test_custom_both(self):
        # use custom for both
        u = MDAnalysis.Universe(TRZ_psf, TRZ, format=MDAnalysis.coordinates.TRZ.TRZReader,
                                topology_format=MDAnalysis.topology.PSFParser.PSFParser)
        assert_equal(len(u.atoms), 8184)


class TestInMemoryUniverse(TestCase):

    @staticmethod
    @dec.skipif(parser_not_found('DCD'),
               'DCD parser not available. Are you using python 3?')
    def test_reader_w_timeseries():
        universe = MDAnalysis.Universe(PSF, DCD, in_memory=True)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 98, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    def test_reader_wo_timeseries():
        universe = MDAnalysis.Universe(GRO, TRR, in_memory=True)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (47681, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    @dec.skipif(parser_not_found('DCD'),
               'DCD parser not available. Are you using python 3?')
    def test_reader_w_timeseries_frame_interval():
        universe = MDAnalysis.Universe(PSF, DCD, in_memory=True,
                                       in_memory_frame_interval=10)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    def test_reader_wo_timeseries_frame_interval():
        universe = MDAnalysis.Universe(GRO, TRR, in_memory=True,
                                       in_memory_frame_interval=3)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (47681, 4, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_existing_universe():
        universe = MDAnalysis.Universe(PDB_small, DCD)
        universe.transfer_to_memory()
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 98, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_frame_interval_convention():
        universe1 = MDAnalysis.Universe(PSF, DCD)
        array1 = universe1.trajectory.timeseries(skip=10)
        universe2 = MDAnalysis.Universe(PSF, DCD, in_memory=True,
                                        in_memory_frame_interval=10)
        array2 = universe2.trajectory.timeseries()
        assert_equal(array1, array2,
                     err_msg="Unexpected differences between arrays.")


class TestWrap(TestCase):
    # VALID
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def setUp(self):
        self.u = MDAnalysis.Universe(TRZ_psf, TRZ)
        self.ag = self.u.atoms[:100]

    # VALID
    def tearDown(self):
        del self.u
        del self.ag

    # VALID
    def test_wrap_comp_fail(self):
        assert_raises(ValueError, self.ag.wrap, compound='strawberries')

    # VALID
    def test_wrap_cent_fail(self):
        assert_raises(ValueError, self.ag.wrap, compound='residues', center='avacado')

    # VALID
    def test_wrap_box_fail(self):
        assert_raises(ValueError, self.ag.wrap, box=np.array([0, 1]))

    def _in_box(self, coords):
        """Check that a set of coordinates are 0.0 <= r <= box"""
        box = self.u.dimensions[:3]

        return (coords >= 0.0).all() and (coords <= box).all()

    # VALID
    def test_wrap_atoms(self):
        ag = self.u.atoms[100:200]
        ag.wrap(compound='atoms')

        assert_equal(self._in_box(ag.positions), True)

    # VALID
    def test_wrap_group(self):
        ag = self.u.atoms[:100]
        ag.wrap(compound='group')

        cen = ag.center_of_mass()

        assert_equal(self._in_box(cen), True)

    # INVALID: must do `r.atoms.center_of_mass()`
    @skip
    def test_wrap_residues(self):
        ag = self.u.atoms[300:400]
        ag.wrap(compound='residues')

        cen = np.vstack([r.center_of_mass() for r in ag.residues])

        assert_equal(self._in_box(cen), True)

    # INVALID: must do `s.atoms.center_of_mass()`
    @skip
    def test_wrap_segments(self):
        ag = self.u.atoms[1000:1200]
        ag.wrap(compound='segments')

        cen = np.vstack([s.center_of_mass() for s in ag.segments])

        assert_equal(self._in_box(cen), True)

    # VALID
    def test_wrap_fragments(self):
        ag = self.u.atoms[:250]
        ag.wrap(compound='fragments')

        cen = np.vstack([f.center_of_mass() for f in ag.fragments])

        assert_equal(self._in_box(cen), True)


class TestGuessBonds(TestCase):
    """Test the AtomGroup methed guess_bonds

    This needs to be done both from Universe creation (via kwarg) and AtomGroup

    It needs to:
     - work if all atoms are in vdwradii table
     - fail properly if not
     - work again if vdwradii are passed.
    """
    # VALID
    def setUp(self):
        self.vdw = {'A':1.05, 'B':0.4}

    # VALID
    def tearDown(self):
        del self.vdw

    # INVALID: Universe has no bonds; AtomGroup does
    def _check_universe(self, u):
        """Verify that the Universe is created correctly"""
        assert_equal(len(u.bonds), 4)
        assert_equal(len(u.angles), 2)
        assert_equal(len(u.dihedrals), 0)
        assert_equal(len(u.atoms[0].bonds), 2)
        assert_equal(len(u.atoms[1].bonds), 1)
        assert_equal(len(u.atoms[2].bonds), 1)
        assert_equal(len(u.atoms[3].bonds), 2)
        assert_equal(len(u.atoms[4].bonds), 1)
        assert_equal(len(u.atoms[5].bonds), 1)
        assert_('guess_bonds' in u.kwargs)

    def test_universe_guess_bonds(self):
        """Test that making a Universe with guess_bonds works"""
        u = MDAnalysis.Universe(two_water_gro, guess_bonds=True)
        self._check_universe(u)
        assert_(u.kwargs['guess_bonds'] is True)

    def test_universe_guess_bonds_no_vdwradii(self):
        """Make a Universe that has atoms with unknown vdwradii."""
        assert_raises(ValueError, MDAnalysis.Universe, two_water_gro_nonames, guess_bonds=True)

    def test_universe_guess_bonds_with_vdwradii(self):
        """Unknown atom types, but with vdw radii here to save the day"""
        u = MDAnalysis.Universe(two_water_gro_nonames, guess_bonds=True,
                                vdwradii=self.vdw)
        self._check_universe(u)
        assert_(u.kwargs['guess_bonds'] is True)
        assert_equal(self.vdw, u.kwargs['vdwradii'])

    def test_universe_guess_bonds_off(self):
        u = MDAnalysis.Universe(two_water_gro_nonames, guess_bonds=False)

        for attr in ('bonds', 'angles', 'dihedrals'):
            assert_(not hasattr(u, attr))
        assert_(u.kwargs['guess_bonds'] is False)

    def _check_atomgroup(self, ag, u):
        """Verify that the AtomGroup made bonds correctly,
        and that the Universe got all this info
        """
        assert_equal(len(ag.bonds), 2)
        assert_equal(len(ag.angles), 1)
        assert_equal(len(ag.dihedrals), 0)
        assert_equal(len(u.bonds), 2)
        assert_equal(len(u.angles), 1)
        assert_equal(len(u.dihedrals), 0)
        assert_equal(len(u.atoms[0].bonds), 2)
        assert_equal(len(u.atoms[1].bonds), 1)
        assert_equal(len(u.atoms[2].bonds), 1)
        assert_equal(len(u.atoms[3].bonds), 0)
        assert_equal(len(u.atoms[4].bonds), 0)
        assert_equal(len(u.atoms[5].bonds), 0)

    def test_atomgroup_guess_bonds(self):
        """Test an atomgroup doing guess bonds"""
        u = MDAnalysis.Universe(two_water_gro)

        ag = u.atoms[:3]
        ag.guess_bonds()
        self._check_atomgroup(ag, u)

    def test_atomgroup_guess_bonds_no_vdwradii(self):
        u = MDAnalysis.Universe(two_water_gro_nonames)

        ag = u.atoms[:3]
        assert_raises(ValueError, ag.guess_bonds)

    def test_atomgroup_guess_bonds_with_vdwradii(self):
        u = MDAnalysis.Universe(two_water_gro_nonames)

        ag = u.atoms[:3]
        ag.guess_bonds(vdwradii=self.vdw)
        self._check_atomgroup(ag, u)


class TestAtomGroupProperties(object):
    """Test working with the properties of Atoms via AtomGroups

    Check that:
    - getting properties from AG matches the Atom values
    - setting properties from AG changes the Atom
    - setting the property on Atom changes AG
    """
    @staticmethod
    def get_new(att_type):
        """Return enough values to change the small g"""
        if att_type == 'string':
            return ['A', 'B', 'C', 'D', 'E', 'F']
        elif att_type == 'float':
            return np.array([0.001, 0.002, 0.003, 0.005, 0.012, 0.025], dtype=np.float32)
        elif att_type == 'int':
            return [4, 6, 8, 1, 5, 4]

    def _check_ag_matches_atom(self, att, atts, ag):
        """Checking Atomgroup property matches Atoms"""
        # Check that accessing via AtomGroup is identical to doing
        # a list comprehension over AG
        ref = [getattr(atom, att) for atom in ag]

        assert_equal(ref, getattr(ag, atts),
                     err_msg="AtomGroup doesn't match Atoms for property: {0}".format(att))

    def _change_atom_check_ag(self, att, atts, vals, ag):
        """Changing Atom, checking AtomGroup matches this"""
        # Set attributes via Atoms
        for atom, val in zip(ag, vals):
            setattr(atom, att, val)
        # Check that AtomGroup returns new values
        other = getattr(ag, atts)

        assert_equal(vals, other,
                     err_msg="Change to Atoms not reflected in AtomGroup for property: {0}".format(att))

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_attributes(self):
        u = make_Universe(('names', 'resids', 'segids', 'types', 'altLocs',
                           'charges', 'masses', 'radii', 'bfactors',
                           'occupancies'))
        u.atoms.occupancies = 1.0
        master = u.atoms
        idx = [0, 1, 4, 7, 11, 14]
        ag = master[idx]

        for att, atts, att_type in (
                ('name', 'names', 'string'),
                ('type', 'types', 'string'),
                ('altLoc', 'altLocs', 'string'),
                ('charge', 'charges', 'float'),
                ('mass', 'masses', 'float'),
                ('radius', 'radii', 'float'),
                ('bfactor', 'bfactors', 'float'),
                ('occupancy', 'occupancies', 'float')
        ):
            vals = self.get_new(att_type)
            yield self._check_ag_matches_atom, att, atts, ag
            yield self._change_atom_check_ag, att, atts, vals, ag


# VALID
class TestOrphans(object):
    """Test moving Universes out of scope and having A/AG persist

    Atoms and AtomGroups from other scopes should work, namely:
      - should have access to Universe
      - should be able to use the Reader (coordinates)
    """
    # VALID
    def test_atom(self):
        u = MDAnalysis.Universe(two_water_gro)

        def getter():
            u2 = MDAnalysis.Universe(two_water_gro)
            return u2.atoms[1]

        atom = getter()

        assert_(atom is not u.atoms[1])
        assert_(len(atom.universe.atoms) == len(u.atoms))
        assert_array_almost_equal(atom.position, u.atoms[1].position)

    # VALID
    def test_atomgroup(self):
        u = MDAnalysis.Universe(two_water_gro)

        def getter():
            u2 = MDAnalysis.Universe(two_water_gro)
            return u2.atoms[:4]

        ag = getter()
        ag2 = u.atoms[:4]
        assert_(ag is not ag2)
        assert_(len(ag.universe.atoms) == len(u.atoms))
        assert_array_almost_equal(ag.positions, ag2.positions)


class TestCrossUniverse(object):
    """Test behaviour when we mix Universes"""

    # VALID
    def _check_badadd(self, a, b):
        def add(x, y):
            return x + y
        assert_raises(ValueError, add, a, b)

    # VALID: currently gives TypeError and unhelpful message
    def test_add_mixed_universes(self):
        # Issue #532
        # Checks that adding objects from different universes
        # doesn't proceed quietly.
        u1 = MDAnalysis.Universe(two_water_gro)
        u2 = MDAnalysis.Universe(two_water_gro)

        A = [u1.atoms[:2], u1.atoms[3]]
        B = [u2.atoms[:3], u2.atoms[0]]

        # Checks Atom to Atom, Atom to AG, AG to Atom and AG to AG
        for x, y in itertools.product(A, B):
            yield self._check_badadd, x, y

    def test_adding_empty_ags(self):
        # Check that empty AtomGroups don't trip up on the Universe check
        u = MDAnalysis.Universe(two_water_gro)

        assert_(len(u.atoms[[]] + u.atoms[:3]) == 3)
        assert_(len(u.atoms[:3] + u.atoms[[]]) == 3)


class TestDihedralSelections(object):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.dih_prec = 2

    def tearDown(self):
        del self.universe
        del self.dih_prec

    def test_phi_selection(self):
        phisel = self.universe.s4AKE.r10.phi_selection()
        assert_equal(phisel.names, ['C', 'N', 'CA', 'C'])
        assert_equal(phisel.residues.resids, [9, 10])
        assert_equal(phisel.residues.resnames, ['PRO', 'GLY'])

    def test_psi_selection(self):
        psisel = self.universe.s4AKE.r10.psi_selection()
        assert_equal(psisel.names, ['N', 'CA', 'C', 'N'])
        assert_equal(psisel.residues.resids, [10, 11])
        assert_equal(psisel.residues.resnames, ['GLY', 'ALA'])

    def test_omega_selection(self):
        osel = self.universe.s4AKE.r8.omega_selection()
        assert_equal(osel.names, ['CA', 'C', 'N', 'CA'])
        assert_equal(osel.residues.resids, [8, 9])
        assert_equal(osel.residues.resnames, ['ALA', 'PRO'])

    def test_chi1_selection(self):
        sel = self.universe.s4AKE.r13.chi1_selection()  # LYS
        assert_equal(sel.names, ['N', 'CA', 'CB', 'CG'])
        assert_equal(sel.residues.resids, [13])
        assert_equal(sel.residues.resnames, ['LYS'])

    def test_phi_sel_fail(self):
        sel = self.universe.residues[0].phi_selection()
        assert_equal(sel, None)

    def test_psi_sel_fail(self):
        sel = self.universe.residues[-1].psi_selection()
        assert_equal(sel, None)

    def test_omega_sel_fail(self):
        sel = self.universe.residues[-1].omega_selection()
        assert_equal(sel, None)

    def test_ch1_sel_fail(self):
        sel = self.universe.s4AKE.r8.chi1_selection()
        assert_equal(sel, None)  # ALA

    def test_dihedral_phi(self):
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        phisel = u.s4AKE.r10.phi_selection()
        assert_almost_equal(phisel.dihedral.value(), -168.57384, self.dih_prec)

    def test_dihedral_psi(self):
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        psisel = u.s4AKE.r10.psi_selection()
        assert_almost_equal(psisel.dihedral.value(), -30.064838, self.dih_prec)

    def test_dihedral_omega(self):
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        osel = u.s4AKE.r8.omega_selection()
        assert_almost_equal(osel.dihedral.value(), -179.93439, self.dih_prec)

    def test_dihedral_chi1(self):
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        sel = u.s4AKE.r13.chi1_selection()  # LYS
        assert_almost_equal(sel.dihedral.value(), -58.428127, self.dih_prec)
