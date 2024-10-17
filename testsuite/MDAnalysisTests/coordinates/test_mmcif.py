import glob

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


@pytest.mark.parametrize("mmcif_filename", glob.glob(f"{MMCIF_FOLDER}/*.cif*"))
def test_works_with_explicit_format(mmcif_filename):
    u = mda.Universe(mmcif_filename, format="MMCIF")
    assert u.trajectory.n_atoms > 0


@pytest.mark.parametrize("mmcif_filename", glob.glob(f"{MMCIF_FOLDER}/*.cif*"))
def test_works_without_explicit_format(mmcif_filename):
    u = mda.Universe(mmcif_filename)
    assert u.trajectory.n_atoms > 0


@pytest.mark.parametrize(
    "mmcif_filename,natoms_protein,natoms_total",
    [
        (f"{MMCIF_FOLDER}/1YJP.cif", 59, 66),
        (f"{MMCIF_FOLDER}/1YJP.cif.gz", 59, 66),
        (f"{MMCIF_FOLDER}/7ETN.cif", 150, 150),
        (f"{MMCIF_FOLDER}/7ETN.cif.gz", 150, 150),
    ],
)
def test_n_atoms(mmcif_filename, natoms_protein, natoms_total):
    u = mda.Universe(mmcif_filename)
    assert len(u.atoms) == natoms_total
    assert len(u.select_atoms("protein").atoms) == natoms_protein


@pytest.mark.parametrize(
    "mmcif_filename,cell",
    [
        (
            f"{MMCIF_FOLDER}/1YJP.cif.gz",
            np.array([21.937, 4.866, 23.477, 90.00, 107.08, 90.00]),
        ),
        (
            f"{MMCIF_FOLDER}/7ETN.cif.gz",
            np.array([5.264, 24.967, 20.736, 90.00, 94.85, 90.00]),
        ),
    ],
)
def test_cell(mmcif_filename, cell):
    assert np.allclose(mda.Universe(mmcif_filename).coord._unitcell, cell)
