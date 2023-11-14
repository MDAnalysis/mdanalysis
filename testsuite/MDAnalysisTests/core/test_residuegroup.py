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
import pickle
from numpy.testing import assert_equal, assert_almost_equal
import pytest

import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import HAS_BIOPYTHON
from MDAnalysisTests.datafiles import PSF, DCD


@pytest.mark.skipif(HAS_BIOPYTHON, reason="biopython is installed")
def test_sequence_import_error():
    p = mda.Universe(PSF, DCD).select_atoms('protein')
    errmsg = "The `sequence_alignment` method requires an installation"
    with pytest.raises(ImportError, match=errmsg):
        _ = p.residues.sequence(format="string")


@pytest.mark.skipif(not HAS_BIOPYTHON, reason='requires biopython')
class TestSequence:
    # all tests are done with the AdK system (PSF and DCD) sequence:
    # http://www.uniprot.org/uniprot/P69441.fasta
    # >sp|P69441|KAD_ECOLI Adenylate kinase OS=Escherichia coli (strain K12) GN=adk PE=1 SV=1
    ref_adk_sequence = (
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT"
        "DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI"
        "VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG"
        "YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"
    )

    @pytest.fixture()
    def u(self):
        return mda.Universe(PSF, DCD)

    def test_string(self, u):
        p = u.select_atoms("protein")
        assert_equal(p.residues.sequence(format="string"),
                     self.ref_adk_sequence)

    def test_SeqRecord(self, u):
        p = u.select_atoms("protein")
        s = p.residues.sequence(format="SeqRecord",
                                id="P69441", name="KAD_ECOLI Adenylate kinase",
                                description="EcAdK from pdb 4AKE")
        assert_equal(s.id, "P69441")
        assert_equal(str(s.seq), self.ref_adk_sequence)

    def test_SeqRecord_default(self, u):
        p = u.select_atoms("protein")
        s = p.residues.sequence(id="P69441", name="KAD_ECOLI Adenylate kinase",
                                description="EcAdK from pdb 4AKE")
        assert_equal(s.id, "P69441")
        assert_equal(str(s.seq), self.ref_adk_sequence)

    def test_Seq(self, u):
        p = u.select_atoms("protein")
        s = p.residues.sequence(format="Seq")
        assert_equal(str(s), self.ref_adk_sequence)

    def test_nonIUPACresname_VE(self, u):
        """test_sequence_nonIUPACresname: non recognized amino acids raise
        ValueError"""
        # fake non-IUPAC residue name for this test
        residues = u.select_atoms("resname MET").residues
        residues.resnames = "MSE"

        def wrong_res():
            u.residues.sequence()

        with pytest.raises(ValueError):
            wrong_res()

    def test_format_TE(self, u):
        with pytest.raises(TypeError):
            u.residues.sequence(format='chicken')


class TestResidueGroup(object):
    # Legacy tests from before 363
    @pytest.fixture()
    def universe(self):
        return mda.Universe(PSF, DCD)

    @pytest.fixture()
    def rg(self, universe):
        return universe.residues

    def test_newResidueGroup(self, universe):
        """test that slicing a ResidueGroup returns a new ResidueGroup
        (Issue 135)"""
        rg = universe.atoms.residues
        newrg = rg[10:20:2]
        assert isinstance(newrg, mda.core.groups.ResidueGroup), \
            "Failed to make a new ResidueGroup: type mismatch"

    def test_n_atoms(self, rg):
        assert_equal(rg.n_atoms, 3341)

    def test_n_residues(self, rg):
        assert_equal(rg.n_residues, 214)

    def test_resids_dim(self, rg):
        assert_equal(len(rg.resids), len(rg))

    def test_resnums_dim(self, rg):
        assert_equal(len(rg.resnums), len(rg))

    def test_segids_dim(self, rg):
        assert_equal(len(rg.segids), len(rg))

    def test_len(self, rg):
        """testing that len(residuegroup) == residuegroup.n_residues"""
        assert_equal(len(rg), rg.n_residues,
                     "len and n_residues disagree")

    def test_set_resids(self, universe):
        rg = universe.select_atoms("bynum 12:42").residues
        resid = 999
        rg.resids = resid
        # check individual atoms
        for at in rg.atoms:
            assert_equal(at.resid, resid,
                         err_msg="failed to set_resid atoms 12:42 to same resid")
        # check residues
        assert_equal(rg.resids, resid * np.ones(rg.n_residues),
                     err_msg="failed to set_resid of residues belonging to "
                             "atoms 12:42 to same resid")

    def test_set_resids(self, universe):
        """test_set_resid: set ResidueGroup resids on a per-residue basis"""
        rg = universe.select_atoms("resid 10:18").residues
        resids = np.array(rg.resids) + 1000
        rg.resids = resids
        # check individual atoms
        for r, resid in zip(rg, resids):
            for at in r.atoms:
                assert_equal(at.resid, resid,
                             err_msg="failed to set_resid residues 10:18 to same "
                                     "resid in residue {0}\n"
                                     "(resids = {1}\nresidues = {2})".format(r,
                                                                             resids,
                                                                             rg))
        assert_equal(rg.resids, resids,
                     err_msg="failed to set_resid of residues belonging to "
                             "residues 10:18 to new resids")

    def test_set_resids_updates_self(self, universe):
        rg = universe.select_atoms("resid 10:18").residues
        resids = np.array(rg.resids) + 1000
        rg.resids = resids
        assert_equal(rg.resids, resids,
                     err_msg="old selection was not changed in place "
                             "after set_resid")

    def test_set_resnum_single(self, universe):
        rg = universe.residues[:3]
        new = 22
        rg.resnums = new

        assert_equal(all(rg.resnums == new), True)
        for r in rg:
            assert_equal(r.resnum, new)

    def test_set_resnum_many(self, universe):
        rg = universe.residues[:3]
        new = [22, 23, 24]
        rg.resnums = new

        assert_equal(all(rg.resnums == new), True)
        for r, v in zip(rg, new):
            assert_equal(r.resnum, v)

    def test_set_resnum_ValueError(self, universe):
        rg = universe.residues[:3]
        new = [22, 23, 24, 25]
        with pytest.raises(ValueError):
            setattr(rg, 'resnums', new)

    # INVALID: no `set_resnames` method; use `resnames` property directly
    @pytest.mark.skip
    def test_set_resname_single(self, universe):
        rg = universe.residues[:3]
        new = 'newname'

        rg.set_resnames(new)
        assert_equal(all(rg.resnames == new), True)
        for r in rg:
            assert_equal(r.name, new)

    # INVALID: no `set_resnames` method; use `resnames` property directly
    @pytest.mark.skip
    def test_set_resname_many(self, universe):
        rg = universe.residues[:3]
        new = ['a', 'b', 'c']
        rg.set_resnames(new)

        assert_equal(all(rg.resnames == new), True)
        for r, v in zip(rg, new):
            assert_equal(r.name, v)

    # INVALID: no `set_resnames` method; use `resnames` property directly
    @pytest.mark.skip
    def test_set_resname_ValueError(self, universe):
        rg = universe.residues[:3]
        new = ['a', 'b', 'c', 'd']

        with pytest.raises(ValueError):
            rg.set_resnames(new)

    # INVALID: no `set_resids` method; also, residues are not mergeable
    # by setting resids; resids are not necessarily unique; atoms must
    # have their resindex set to change residue membership
    @pytest.mark.skip
    def test_merge_residues(self, universe):
        rg = universe.select_atoms("resid 12:14").residues
        nres_old = universe.atoms.n_residues
        natoms_old = rg.n_atoms
        rg.set_resids(12)  # merge all into one with resid 12
        nres_new = universe.atoms.n_residues
        r_merged = universe.select_atoms("resid 12:14").residues
        natoms_new = universe.select_atoms("resid 12").n_atoms
        assert_equal(len(r_merged), 1, err_msg="set_resid failed to merge "
                                               "residues: merged = {0}".format(
            r_merged))
        assert_equal(nres_new, nres_old - 2,
                     err_msg="set_resid failed to merge residues: "
                             "merged = {0}".format(r_merged))
        assert_equal(natoms_new, natoms_old, err_msg="set_resid lost atoms "
                                                     "on merge".format(
            r_merged))

        assert_equal(universe.residues.n_residues,
                     universe.atoms.n_residues,
                     err_msg="Universe.residues and Universe.atoms.n_residues "
                             "do not agree after residue "
                             "merge.")

    # INVALID: no `set_masses` method; use `masses` property directly
    @pytest.mark.skip
    def test_set_masses(self, universe):
        rg = universe.select_atoms("bynum 12:42 and name H*").residues
        mass = 2.0
        rg.set_masses(mass)
        # check individual atoms
        assert_equal([a.mass for a in rg.atoms],
                     mass * np.ones(rg.n_atoms),
                     err_msg="failed to set_mass H* atoms in resid 12:42 to {0}".format(
                         mass))

    # VALID
    def test_atom_order(self, universe):
        assert_equal(universe.residues.atoms.indices,
                     sorted(universe.residues.atoms.indices))

    def test_get_next_residue(self, rg):
        unsorted_rep_res = rg[[0, 1, 8, 3, 4, 0, 3, 1, -1]]
        next_res = unsorted_rep_res._get_next_residues_by_resid()
        resids = list(unsorted_rep_res.resids+1)
        resids[-1] = None
        next_resids = [r.resid if r is not None else None for r in next_res]
        assert_equal(len(next_res), len(unsorted_rep_res))
        assert_equal(next_resids, resids)

    def test_get_prev_residue(self, rg):
        unsorted_rep_res = rg[[0, 1, 8, 3, 4, 0, 3, 1, -1]]
        prev_res = unsorted_rep_res._get_prev_residues_by_resid()
        resids = list(unsorted_rep_res.resids-1)
        resids[0] = resids[5] = None
        prev_resids = [r.resid if r is not None else None for r in prev_res]
        assert_equal(len(prev_res), len(unsorted_rep_res))
        assert_equal(prev_resids, resids)

    @pytest.mark.parametrize("selection", ("name CA", "segid 4AKE"))
    def test_residuegroup_pickle(self, universe, selection):
        seg_res = universe.select_atoms(selection).residues
        seg = pickle.loads(pickle.dumps(seg_res))
        assert_almost_equal(seg_res.atoms.positions, seg.atoms.positions)
