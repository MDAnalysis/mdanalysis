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

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import (PSF, DCD, PDB_small, GRO, TRR,
                                        TRZ, TRZ_psf, PSF_notop,
                                        XYZ_mini, two_water_gro,
                                        two_water_gro_nonames,
                                        COORDINATES_XYZ, COORDINATES_TRR)
import MDAnalysis.core.groups
from MDAnalysis.core.groups import Atom, AtomGroup
from MDAnalysis import NoDataError

import numpy as np
from numpy.testing import (dec, assert_equal,
                           assert_almost_equal, assert_raises, assert_,
                           assert_array_almost_equal, assert_array_equal,
                           assert_allclose)
from nose.plugins.attrib import attr

from unittest import skip

import os
import itertools

from MDAnalysisTests import parser_not_found, tempdir, make_Universe


class TestUniverseSetTopology(object):
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


class _WriteAtoms(object):
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


class TestSequence(object):
    # all tests are done with the AdK system (PSF and DCD) sequence:
    # http://www.uniprot.org/uniprot/P69441.fasta
    # >sp|P69441|KAD_ECOLI Adenylate kinase OS=Escherichia coli (strain K12) GN=adk PE=1 SV=1
    ref_adk_sequence = (
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT"
        "DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI"
        "VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG"
        "YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"
    )

    def setUp(self):
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    # INVALID: AtomGroup no longer gets a sequence method; only in ResidueGroup
    # if resnames in topology
    @skip
    def test_sequence_from_atoms(self):
        p = self.u.select_atoms("protein")
        assert_equal(p.sequence(format="string"),
                     p.residues.sequence(format="string"),
                     err_msg="sequence() yields different results for "
                     "residues and atoms")

    def test_sequence_string(self):
        p = self.u.select_atoms("protein")
        assert_equal(p.residues.sequence(format="string"),
                     self.ref_adk_sequence)

    def test_sequence_SeqRecord(self):
        p = self.u.select_atoms("protein")
        s = p.residues.sequence(format="SeqRecord",
                                id="P69441", name="KAD_ECOLI Adenylate kinase",
                                description="EcAdK from pdb 4AKE")
        assert_equal(s.id, "P69441")
        assert_equal(s.seq.tostring(), self.ref_adk_sequence)

    def test_sequence_SeqRecord_default(self):
        p = self.u.select_atoms("protein")
        s = p.residues.sequence(id="P69441", name="KAD_ECOLI Adenylate kinase",
                                description="EcAdK from pdb 4AKE")
        assert_equal(s.id, "P69441")
        assert_equal(s.seq.tostring(), self.ref_adk_sequence)

    def test_sequence_Seq(self):
        p = self.u.select_atoms("protein")
        s = p.residues.sequence(format="Seq")
        assert_equal(s.tostring(), self.ref_adk_sequence)

    @skip
    def test_sequence_nonIUPACresname(self):
        """test_sequence_nonIUPACresname: non recognized amino acids raise
        ValueError"""
        # fake non-IUPAC residue name for this test
        residues = self.u.select_atoms("resname MET").residues
        residues.resnames = "MSE"

        def wrong_res():
            self.u.atoms.sequence()

        assert_raises(ValueError, wrong_res)
