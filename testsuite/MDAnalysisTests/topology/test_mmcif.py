import glob
import os
from io import StringIO

import MDAnalysis as mda
import numpy as np
import pytest
from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_array_almost_equal,
    assert_equal,
)

from MDAnalysisTests.datafiles import MMCIF as MMCIF_FOLDER

# FIXME: rewrite tests to read trajectories only once


@pytest.mark.parametrize(
    "mmcif_filename,n_chains",
    [
        (f"{MMCIF_FOLDER}/1YJP.cif", 1),
        (f"{MMCIF_FOLDER}/1YJP.cif.gz", 1),
        (f"{MMCIF_FOLDER}/7ETN.cif", 2),
        (f"{MMCIF_FOLDER}/7ETN.cif.gz", 2),
    ],
)
def test_chains(mmcif_filename, n_chains):
    u = mda.Universe(mmcif_filename)
    assert len(u.segments) == n_chains


@pytest.mark.parametrize(
    "mmcif_filename,sequence",
    [
        (f"{MMCIF_FOLDER}/1YJP.cif", ["GLY", "ASN", "ASN", "GLN", "GLN", "ASN", "TYR"]),
        (
            f"{MMCIF_FOLDER}/1YJP.cif.gz",
            ["GLY", "ASN", "ASN", "GLN", "GLN", "ASN", "TYR"],
        ),
        (f"{MMCIF_FOLDER}/7ETN.cif", ["PRO", "PHE", "LEU", "ILE"]),
        (f"{MMCIF_FOLDER}/7ETN.cif.gz", ["PRO", "PHE", "LEU", "ILE"]),
    ],
)
def test_sequence(mmcif_filename, sequence):
    u = mda.Universe(mmcif_filename)
    in_structure = [
        str(res.resname) for res in u.select_atoms("protein and chainid A").residues
    ]
    assert in_structure == sequence, ":".join(in_structure)
