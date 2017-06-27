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
from __future__ import absolute_import

import numpy as np
from numpy.testing import (
    dec,
    assert_,
    assert_equal,
    assert_raises,
)
from unittest import skip

import MDAnalysis as mda

from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysisTests import parser_not_found


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

    def test_string(self):
        p = self.u.select_atoms("protein")
        assert_equal(p.residues.sequence(format="string"),
                     self.ref_adk_sequence)

    def test_SeqRecord(self):
        p = self.u.select_atoms("protein")
        s = p.residues.sequence(format="SeqRecord",
                                id="P69441", name="KAD_ECOLI Adenylate kinase",
                                description="EcAdK from pdb 4AKE")
        assert_equal(s.id, "P69441")
        assert_equal(s.seq.tostring(), self.ref_adk_sequence)

    def test_SeqRecord_default(self):
        p = self.u.select_atoms("protein")
        s = p.residues.sequence(id="P69441", name="KAD_ECOLI Adenylate kinase",
                                description="EcAdK from pdb 4AKE")
        assert_equal(s.id, "P69441")
        assert_equal(s.seq.tostring(), self.ref_adk_sequence)

    def test_Seq(self):
        p = self.u.select_atoms("protein")
        s = p.residues.sequence(format="Seq")
        assert_equal(s.tostring(), self.ref_adk_sequence)

    def test_nonIUPACresname_VE(self):
        """test_sequence_nonIUPACresname: non recognized amino acids raise
        ValueError"""
        # fake non-IUPAC residue name for this test
        residues = self.u.select_atoms("resname MET").residues
        residues.resnames = "MSE"

        def wrong_res():
            self.u.residues.sequence()

        assert_raises(ValueError, wrong_res)

    def test_format_TE(self):
        assert_raises(TypeError, self.u.residues.sequence, format='chicken')


class TestResidueGroup(object):
    # Legacy tests from before 363
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = mda.Universe(PSF, DCD)
        self.rg = self.universe.residues

    def test_newResidueGroup(self):
        """test that slicing a ResidueGroup returns a new ResidueGroup
        (Issue 135)"""
        rg = self.universe.atoms.residues
        newrg = rg[10:20:2]
        assert_(isinstance(newrg, mda.core.groups.ResidueGroup),
                "Failed to make a new ResidueGroup: type mismatch")

    def test_n_atoms(self):
        assert_equal(self.rg.n_atoms, 3341)

    def test_n_residues(self):
        assert_equal(self.rg.n_residues, 214)

    def test_resids_dim(self):
        assert_equal(len(self.rg.resids), len(self.rg))

    def test_resnums_dim(self):
        assert_equal(len(self.rg.resnums), len(self.rg))

    def test_segids_dim(self):
        assert_equal(len(self.rg.segids), len(self.rg))

    def test_len(self):
        """testing that len(residuegroup) == residuegroup.n_residues"""
        assert_equal(len(self.rg), self.rg.n_residues,
                     "len and n_residues disagree")

    def test_set_resids(self):
        rg = self.universe.select_atoms("bynum 12:42").residues
        resid = 999
        rg.resids = resid
        # check individual atoms
        for at in rg.atoms:
            assert_equal(a.resid, resid,
                         err_msg="failed to set_resid atoms 12:42 to same resid")
        # check residues
        assert_equal(rg.resids, resid * np.ones(rg.n_residues),
                     err_msg="failed to set_resid of residues belonging to "
                     "atoms 12:42 to same resid")

    def test_set_resids(self):
        """test_set_resid: set ResidueGroup resids on a per-residue basis"""
        rg = self.universe.select_atoms("resid 10:18").residues
        resids = np.array(rg.resids) + 1000
        rg.resids = resids
        # check individual atoms
        for r, resid in zip(rg, resids):
            for at in r.atoms:
                assert_equal(at.resid, resid,
                             err_msg="failed to set_resid residues 10:18 to same "
                             "resid in residue {0}\n"
                             "(resids = {1}\nresidues = {2})".format(r, resids, rg))
        assert_equal(rg.resids, resids,
                     err_msg="failed to set_resid of residues belonging to "
                     "residues 10:18 to new resids")

    def test_set_resids_updates_self(self):
        rg = self.universe.select_atoms("resid 10:18").residues
        resids = np.array(rg.resids) + 1000
        rg.resids = resids
        assert_equal(rg.resids, resids,
                     err_msg="old selection was not changed in place "
                     "after set_resid")

    def test_set_resnum_single(self):
        rg = self.universe.residues[:3]
        new = 22
        rg.resnums = new

        assert_equal(all(rg.resnums == new), True)
        for r in rg:
            assert_equal(r.resnum, new)

    def test_set_resnum_many(self):
        rg = self.universe.residues[:3]
        new = [22, 23, 24]
        rg.resnums = new

        assert_equal(all(rg.resnums == new), True)
        for r, v in zip(rg, new):
            assert_equal(r.resnum, v)

    def test_set_resnum_ValueError(self):
        rg = self.universe.residues[:3]
        new = [22, 23, 24, 25]

        assert_raises(ValueError, setattr, rg, 'resnums', new)

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


