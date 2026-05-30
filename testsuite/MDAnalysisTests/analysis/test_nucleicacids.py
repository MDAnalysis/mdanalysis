# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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

import MDAnalysis as mda
import numpy as np

import pytest

from pytest import approx

from MDAnalysis.analysis.nucleicacids import (
    NucPairDist,
    WatsonCrickDist,
    MajorPairDist,
    MinorPairDist,
)

from MDAnalysisTests.datafiles import RNA_PSF, RNA_PDB

from MDAnalysis.core.groups import ResidueGroup


@pytest.fixture(scope="module")
def u():
    return mda.Universe(RNA_PSF, RNA_PDB)


@pytest.fixture(scope="module")
def strand(unv=u):
    unv = mda.Universe(RNA_PSF, RNA_PDB)
    return unv.select_atoms("segid RNAA")


def test_empty_ag_error(strand):
    strand1 = ResidueGroup([strand.residues[0]])
    strand2 = ResidueGroup([strand.residues[1]])

    with pytest.raises(ValueError, match="returns an empty AtomGroup"):
        NucPairDist.select_strand_atoms(
            strand1, strand2, "UNK1", "O2", g_name="GUA"
        )


@pytest.fixture(scope="module")
def wc_rna(strand, client_NucPairDist):
    strand1 = ResidueGroup([strand.residues[0], strand.residues[21]])
    strand2 = ResidueGroup([strand.residues[1], strand.residues[22]])

    WC = WatsonCrickDist(
        strand1,
        strand2,
        g_name="GUA",
        a_name="ADE",
        c_name="CYT",
        u_name="URA",
    )
    WC.run(**client_NucPairDist)
    return WC


def test_wc_dist_shape(wc_rna):
    assert wc_rna.results.distances.shape == (1, 2)


def test_wc_dist_results_keys(wc_rna):
    assert "distances" in wc_rna.results


def test_wc_dist(wc_rna):
    assert wc_rna.results.distances[0, 0] == approx(4.3874702, rel=1e-3)
    assert wc_rna.results.distances[0, 1] == approx(4.1716404, rel=1e-3)


def test_wc_dist_invalid_residue_types(u):
    strand = u.select_atoms("resid 1-10")
    # residues[0]=GUA, residues[21]=POT (non-nucleic), residues[22]=POT
    strand1 = ResidueGroup([strand.residues[0], strand.residues[21]])
    strand2 = ResidueGroup([strand.residues[2], strand.residues[22]])
    with pytest.raises(ValueError, match="is not a valid nucleic acid"):
        WatsonCrickDist(
            strand1,
            strand2,
            g_name="GUA",
            a_name="ADE",
            c_name="CYT",
            u_name="URA",
        )


def test_selection_length_mismatch(strand):
    sel1 = strand.residues[1:10]
    sel2 = strand.residues[1:9]
    with pytest.raises(ValueError, match="Selections must be same length"):
        NucPairDist(sel1, sel2)


def test_wc_dist_deprecation_warning(strand):
    strand1 = [strand.residues[0], strand.residues[21]]
    strand2 = [strand.residues[1], strand.residues[22]]

    with pytest.deprecated_call():
        WatsonCrickDist(
            strand1,
            strand2,
            g_name="GUA",
            a_name="ADE",
            c_name="CYT",
            u_name="URA",
        )


def test_wc_dist_strand_verification(strand):
    strand1 = [strand.residues[0], strand[0]]
    strand2 = [strand.residues[1], strand.residues[22]]

    with pytest.raises(TypeError, match="contains non-Residue elements"):
        WatsonCrickDist(strand1, strand2)


@pytest.mark.parametrize("key", [0, 1, 2, "parsnips", "time", -1])
def test_wc_dis_results_keyerrs(wc_rna, key):
    with pytest.raises(KeyError, match=f"{key}"):
        wc_rna.results[key]


def test_minor_dist(strand, client_NucPairDist):
    strand1 = ResidueGroup([strand.residues[2], strand.residues[19]])
    strand2 = ResidueGroup([strand.residues[16], strand.residues[4]])

    MI = MinorPairDist(
        strand1,
        strand2,
        g_name="GUA",
        a_name="ADE",
        c_name="CYT",
        u_name="URA",
    )
    MI.run(**client_NucPairDist)

    assert MI.results.distances[0, 0] == approx(15.06506, rel=1e-3)
    assert MI.results.distances[0, 1] == approx(3.219116, rel=1e-3)


def test_major_dist(strand, client_NucPairDist):
    strand1 = ResidueGroup([strand.residues[1], strand.residues[4]])
    strand2 = ResidueGroup([strand.residues[11], strand.residues[8]])

    MA = MajorPairDist(
        strand1,
        strand2,
        g_name="GUA",
        a_name="ADE",
        c_name="CYT",
        u_name="URA",
    )
    MA.run(**client_NucPairDist)

    assert MA.results.distances[0, 0] == approx(26.884272, rel=1e-3)
    assert MA.results.distances[0, 1] == approx(13.578535, rel=1e-3)


@pytest.fixture(scope="module")
def dna_u():
    """Synthetic universe with 2-letter DNA resnames (DG, DC, DA, DT).

    Used to test that WatsonCrickDist matches full resnames rather than only
    the first character.

    Strand layout (one frame, all atoms at origin):
      residue 0 – DG (purine):    atom N1
      residue 1 – DC (pyrimidine): atom N3
      residue 2 – DA (purine):    atom N1
      residue 3 – DT (pyrimidine): atom N3
    """
    n_atoms = 4
    n_residues = 4
    u = mda.Universe.empty(
        n_atoms,
        n_residues=n_residues,
        atom_resindex=[0, 1, 2, 3],
        trajectory=True,
    )
    u.add_TopologyAttr("name", ["N1", "N3", "N1", "N3"])
    u.add_TopologyAttr("resname", ["DG", "DC", "DA", "DT"])
    u.add_TopologyAttr("resid", [1, 2, 3, 4])
    coords = np.zeros((1, n_atoms, 3), dtype=np.float32)
    u.load_new(coords)
    return u


def test_wc_dist_multichar_resnames(dna_u):
    """WatsonCrickDist must match the full resname, not just resname[0]\.

    With 2-letter names like DG/DA, the first character 'D' does not appear in
    the purine/pyrimidine lists, so the buggy resname[0] check previously
    rasied ValueError for every residue.
    """
    strand1 = ResidueGroup([dna_u.residues[0], dna_u.residues[2]])  # DG, DA
    strand2 = ResidueGroup([dna_u.residues[1], dna_u.residues[3]])  # DC, DT

    WC = WatsonCrickDist(
        strand1,
        strand2,
        g_name="DG",
        a_name="DA",
        c_name="DC",
        t_name="DT",
    )
    WC.run()
    assert WC.results.distances.shape == (1, 2)


def test_select_strand_atoms_multichar_resnames(dna_u):
    """select_strand_atoms must recognise multi-character resnames."""
    strand1 = ResidueGroup([dna_u.residues[0], dna_u.residues[2]])  # DG, DA
    strand2 = ResidueGroup([dna_u.residues[1], dna_u.residues[3]])  # DC, DT

    sel1, sel2 = NucPairDist.select_strand_atoms(
        strand1,
        strand2,
        a1_name="N1",
        a2_name="N3",
        g_name="DG",
        a_name="DA",
        c_name="DC",
        t_name="DT",
    )
    assert len(sel1) == 2
    assert len(sel2) == 2
