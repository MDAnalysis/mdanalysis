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
from __future__ import absolute_import, print_function

import os
import MDAnalysis
import pytest
from MDAnalysis.tests.datafiles import (
    PSF,
    DCD,
    capping_input,
    capping_output,
    capping_ace,
    capping_nma,
    merge_protein,
    merge_ligand,
    merge_water
)
import MDAnalysis.core.groups
from MDAnalysis.core.groups import AtomGroup

from numpy.testing import assert_equal, assert_array_equal

from MDAnalysis import Merge
from MDAnalysis.analysis.align import alignto


def capping(ref, ace, nma, output):
    resids = ref.select_atoms("all").resids
    resid_min, resid_max = min(resids), max(resids)

    for r in ace.residues:
        r.resid += resid_min - max(ace.atoms.resids)
    for r in nma.residues:
        r.resid = resid_max

    # TODO pick the first residue in the protein (how should we cap the chains?)
    # TODO consider a case when the protein resid is 1 and all peptide has to be shifted by +1, put that in docs as a
    #  post-processing step
    alignto(ace, ref, select={
            "mobile": "resid {0} and backbone".format(resid_min),
            "reference": "resid {0} and backbone".format(resid_min)},
            strict=True)
    alignto(nma, ref, select={
            "mobile": "resid {0} and backbone and not (resname NMA NME)".format(resid_max),
            "reference": "resid {0} and (backbone or name OT2)".format(resid_max)},
            strict=True)

    #  TODO remove the Hydrogen closest to ACE's oxygen
    nma.residues.resids = 16
    u = Merge(ace.select_atoms("resname ACE"),
              ref.select_atoms(
                  "not (resid {0} and name HT*) and not (resid {1} and (name HT* OT1))"
                  "".format(resid_min, resid_max)),
              nma.select_atoms("resname NME NMA"))
    u.trajectory.ts.dimensions = ref.trajectory.ts.dimensions
    u.atoms.write(output)
    return u


class TestCapping(object):

    def test_capping_file(self, tmpdir):
        peptide = MDAnalysis.Universe(capping_input)
        ref = MDAnalysis.Universe(capping_output)
        ace = MDAnalysis.Universe(capping_ace)
        nma = MDAnalysis.Universe(capping_nma)

        outfile = str(tmpdir.join('test.pdb'))
        u = capping(peptide, ace, nma, outfile)

        assert_equal(len(u.select_atoms("not name H*")),
                     len(ref.select_atoms("not name H*")))

        u = MDAnalysis.Universe(outfile)

        ace = u.select_atoms("resname ACE")
        nma = u.select_atoms("resname NMA")

        # Check if the resids were assigned correctly
        assert_equal(ace.resids[0], 1)
        assert_equal(nma.resids[0], 16)

        assert_array_equal(peptide.trajectory.ts.dimensions,
                           u.trajectory.ts.dimensions)

    def test_capping_inmemory(self, tmpdir):
        peptide = MDAnalysis.Universe(capping_input)
        ref = MDAnalysis.Universe(capping_output)
        ace = MDAnalysis.Universe(capping_ace)
        nma = MDAnalysis.Universe(capping_nma)

        outfile = str(tmpdir.join('test.pdb'))
        u = capping(peptide, ace, nma, outfile)
        assert_equal(len(u.select_atoms("not name H*")),
                     len(ref.select_atoms("not name H*")))

        ace = u.select_atoms("resname ACE")
        nma = u.select_atoms("resname NMA")

        # Check if the resids were assigned correctly
        assert_equal(ace.resids[0], 1)
        assert_equal(nma.resids[0], 16)

        assert_array_equal(peptide.trajectory.ts.dimensions,
                           u.trajectory.ts.dimensions)


@pytest.fixture()
def u_protein():
    return MDAnalysis.Universe(merge_protein)


@pytest.fixture()
def u_ligand():
    return MDAnalysis.Universe(merge_ligand)


@pytest.fixture()
def u_water():
    return MDAnalysis.Universe(merge_water)


class TestMerge(object):
    @pytest.mark.xfail(os.name == 'nt',
                       strict=True,
                       reason="Setup fixtures fail on Windows.")
    def test_merge(self, u_protein, u_ligand, u_water, tmpdir):
        ids_before = [a.index for u in [u_protein, u_ligand, u_water] for a in u.atoms]
        # Do the merge
        u0 = MDAnalysis.Merge(u_protein.atoms, u_ligand.atoms, u_water.atoms)
        # Check that the output Universe has the same number of atoms as the
        # starting AtomGroups
        assert_equal(len(u0.atoms), (len(u_protein.atoms) + len(u_ligand.atoms) + len(u_water.atoms)))
        # Check that the output Universe has the same number of residues and
        # segments as the starting AtomGroups
        assert_equal(len(u0.residues), (len(u_protein.residues) + len(u_ligand.residues) +
                                                           len(u_water.residues)))
        assert_equal(len(u0.segments), (len(u_protein.segments) + len(u_ligand.segments) +
                                                           len(u_water.segments)))

        # Make sure that all the atoms in the new universe are assigned to only
        # one, new Universe
        set0 = {a.universe for a in u0.atoms}
        assert_equal(len(set0), 1)
        u = list(set0)[0]
        assert_equal(u, u0)

        # Make sure that the atom ids of the original universes are unchanged,
        # ie we didn't make the original Universes 'dirty'
        ids_after = [a.index for u in [u_protein, u_ligand, u_water] for a in u.atoms]
        assert_equal(len(ids_after), (len(u_protein.atoms) + len(u_ligand.atoms) + len(u_water.atoms)))
        assert_equal(ids_before, ids_after)

        # Test that we have a same number of atoms in a different way
        ids_new = [a.index for a in u0.atoms]
        assert_equal(len(ids_new), len(ids_before))

        outfile = str(tmpdir.join('test.pdb'))

        u0.atoms.write(outfile)
        u = MDAnalysis.Universe(outfile)
        ids_new2 = [a.index for a in u.atoms]
        assert_equal(ids_new, ids_new2)

    @pytest.mark.xfail(os.name == 'nt',
                       strict=True,
                       reason="Setup fixtures fail on Windows.")
    def test_merge_same_universe(self, u_protein):
        u0 = MDAnalysis.Merge(u_protein.atoms, u_protein.atoms, u_protein.atoms)
        assert_equal(len(u0.atoms), 3 * len(u_protein.atoms))
        assert_equal(len(u0.residues), 3 * len(u_protein.residues))
        assert_equal(len(u0.segments), 3 * len(u_protein.segments))

    @pytest.mark.xfail(os.name == 'nt',
                       strict=True,
                       reason="Setup fixtures fail on Windows.")
    def test_residue_references(self, u_protein, u_ligand):
        m = Merge(u_protein.atoms, u_ligand.atoms)
        assert_equal(m.atoms.residues[0].universe, m,
                     "wrong universe reference for residues after Merge()")

    @pytest.mark.xfail(os.name == 'nt',
                       strict=True,
                       reason="Setup fixtures fail on Windows.")
    def test_segment_references(self, u_protein, u_ligand):
        m = Merge(u_protein.atoms, u_ligand.atoms)
        assert_equal(m.atoms.segments[0].universe, m,
                     "wrong universe reference for segments after Merge()")

    def test_empty_ValueError(self):
        with pytest.raises(ValueError):
            Merge()

    def test_nonsense_TypeError(self):
        with pytest.raises(TypeError):
            Merge(['1', 2])

    @pytest.mark.xfail(os.name == 'nt',
                       strict=True,
                       reason="Setup fixtures fail on Windows.")
    def test_emptyAG_ValueError(self, u_protein):
        a = AtomGroup([], u_protein)
        b = AtomGroup([], u_protein)

        with pytest.raises(ValueError):
            Merge(a, b)


class TestMergeTopology(object):
    """Test that Merge correct does topology"""

    @staticmethod
    @pytest.fixture()
    def u():
        return MDAnalysis.Universe(PSF, DCD)

    def test_merge_with_topology(self, u):
        ag1 = u.atoms[:20]
        ag2 = u.atoms[100:110]

        u_merge = MDAnalysis.Merge(ag1, ag2)

        assert(len(u_merge.atoms) == 30)
        assert(len(u_merge.atoms.bonds) == 28)
        assert(len(u_merge.atoms.angles) == 47)
        assert(len(u_merge.atoms.dihedrals) == 53)
        assert(len(u_merge.atoms.impropers) == 1)

        # All these bonds are in the merged Universe
        assert(len(ag1[0].bonds) == len(u_merge.atoms[0].bonds))
        # One of these bonds isn't in the merged Universe
        assert(len(ag2[0].bonds) - 1 == len(u_merge.atoms[20].bonds))

    @pytest.mark.xfail(os.name == 'nt',
                       strict=True,
                       reason="Setup fixtures fail on Windows.")
    def test_merge_with_topology_from_different_universes(self, u, u_ligand):
        u_merge = MDAnalysis.Merge(u.atoms[:110], u_ligand.atoms)

        # merge_protein doesn't contain bond topology, so merged universe
        # shouldn't have one either
        print(u_merge.atoms.bonds)
        # PDB reader yields empty Bonds group, which means bonds from
        # PSF/DCD survive the merge
        # assert(not hasattr(u_merge.atoms, 'bonds') or len(u_merge.atoms.bonds) == 0)
        assert(not hasattr(u_merge.atoms, 'angles') or len(u_merge.atoms.bonds) == 0)
        assert(not hasattr(u_merge.atoms, 'dihedrals') or len(u_merge.atoms.bonds) == 0)
        assert(not hasattr(u_merge.atoms, 'impropers') or len(u_merge.atoms.bonds) == 0)

    def test_merge_without_topology(self, u):
        # This shouldn't have topology as we merged single atoms
        u_merge = MDAnalysis.Merge(u.atoms[0:1], u.atoms[10:11])

        assert(len(u_merge.atoms) == 2)
        assert(len(u_merge.atoms.bonds) == 0)
        assert(len(u_merge.atoms.angles) == 0)
        assert(len(u_merge.atoms.dihedrals) == 0)
        assert(len(u_merge.atoms.impropers) == 0)
