import glob
import warnings
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.analysis.dssp import DSSP, assign, translate
from MDAnalysis.analysis.dssp.pydssp_numpy import _get_hydrogen_atom_position

from MDAnalysisTests.datafiles import DSSP as DSSP_FOLDER
from MDAnalysisTests.datafiles import TPR, XTC


# Files that match glob pattern '????.pdb.gz' and matching '????.pdb.dssp' files,
# containing the secondary structure assignment string, will be tested automatically.
@pytest.mark.parametrize(
    "pdb_filename", glob.glob(f"{DSSP_FOLDER}/?????.pdb.gz")
)
def test_file_guess_hydrogens(pdb_filename, client_DSSP):
    u = mda.Universe(pdb_filename)
    with open(f"{pdb_filename.rstrip('.gz')}.dssp", "r") as fin:
        correct_answ = fin.read().strip().split()[0]

    run = DSSP(u, guess_hydrogens=True).run(**client_DSSP)
    answ = "".join(run.results.dssp[0])
    assert answ == correct_answ


def test_trajectory(client_DSSP):
    u = mda.Universe(TPR, XTC).select_atoms("protein").universe
    run = DSSP(u).run(**client_DSSP, stop=10)
    first_frame = "".join(run.results.dssp[0])
    last_frame = "".join(run.results.dssp[-1])
    avg_frame = "".join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert (
        first_frame[:10] != last_frame[:10] == avg_frame[:10] == "-EEEEEE---"
    )


# ensure that we get different assigned results when using and not using the
# donor mask which filters out prolines from being potential HBond donors as they
# are missing hydrogens
def test_donor_mask():
    u = mda.Universe(TPR, XTC).select_atoms("protein").universe
    dssp = DSSP(u)
    coords = dssp._get_coords()
    n_residues = coords.shape[0]

    result_no_mask = assign(coords)
    result_with_mask = assign(coords, donor_mask=dssp._donor_mask)

    assert result_no_mask.shape == (n_residues, 3)
    assert result_with_mask.shape == (n_residues, 3)
    assert not np.array_equal(result_no_mask, result_with_mask)


def test_donor_mask_proline_indices():
    """Test that donor_mask correctly identifies proline residues."""
    u = mda.Universe(TPR, XTC).select_atoms("protein").universe
    dssp = DSSP(u)
    resnames = u.select_atoms("protein").residues.resnames
    is_pro = resnames == "PRO"

    # donor_mask should be False at proline positions and True at non-proline positions
    assert np.all(~dssp._donor_mask[is_pro])
    assert np.all(dssp._donor_mask[~is_pro])


def test_donor_mask_none():
    """Test that assign works with donor_mask=None (pre-2.10.0 behavior)."""
    u = mda.Universe(TPR, XTC).select_atoms("protein").universe
    dssp = DSSP(u)
    coords = dssp._get_coords()
    n_residues = coords.shape[0]

    # None mask should work and return valid shape and each residue should have exactly one
    # True in the one-hot encoding
    result = assign(coords, donor_mask=None)
    assert result.shape == (n_residues, 3)
    assert np.all(result.sum(axis=-1) == 1)


def test_placeholder_box_warnings():
    """Test that placeholder (1,1,1,90,90,90) box dimensions trigger a warning
    and disable PBC."""
    u = mda.Universe(f"{DSSP_FOLDER}/2va0A.pdb.gz")
    u.trajectory.ts.dimensions = np.array([1.0, 1.0, 1.0, 90.0, 90.0, 90.0])

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        dssp = DSSP(u)
        box_warnings = [
            x for x in w if "not using periodic boundary" in str(x.message)
        ]
        assert len(box_warnings) == 1
    assert dssp._box is None


def test_atomgroup(client_DSSP):
    protein = mda.Universe(TPR, XTC).select_atoms("protein")
    run = DSSP(protein).run(**client_DSSP, stop=10)
    first_frame = "".join(run.results.dssp[0])
    last_frame = "".join(run.results.dssp[-1])
    avg_frame = "".join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert (
        first_frame[:10] != last_frame[:10] == avg_frame[:10] == "-EEEEEE---"
    )


def test_trajectory_with_hydrogens(client_DSSP):
    u = mda.Universe(TPR, XTC).select_atoms("protein").universe
    run = DSSP(u, guess_hydrogens=False).run(**client_DSSP, stop=10)
    first_frame = "".join(run.results.dssp[0])
    last_frame = "".join(run.results.dssp[-1])
    avg_frame = "".join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert (
        first_frame[:10] == last_frame[:10] == avg_frame[:10] == "-EEEEEE---"
    )


@pytest.mark.parametrize(
    "pdb_filename", glob.glob(f"{DSSP_FOLDER}/2xdgA.pdb.gz")
)
def test_trajectory_without_hydrogen_fails(pdb_filename, client_DSSP):
    u = mda.Universe(pdb_filename)
    with pytest.raises(ValueError):
        DSSP(u, guess_hydrogens=False).run(**client_DSSP)


@pytest.mark.parametrize(
    "pdb_filename", glob.glob(f"{DSSP_FOLDER}/1mr1D_failing.pdb.gz")
)
def test_trajectory_with_uneven_number_of_atoms_fails(
    pdb_filename, client_DSSP
):
    u = mda.Universe(pdb_filename)
    with pytest.raises(ValueError):
        DSSP(u, guess_hydrogens=True).run(**client_DSSP)


@pytest.mark.parametrize(
    "pdb_filename", glob.glob(f"{DSSP_FOLDER}/wrong_hydrogens.pdb.gz")
)
def test_exception_raises_with_atom_index(pdb_filename, client_DSSP):
    u = mda.Universe(pdb_filename)
    with pytest.raises(
        ValueError,
        match="Residue <Residue SER, 298> contains*",
    ):
        DSSP(u, guess_hydrogens=False).run(**client_DSSP)


def test_insufficient_residues_raises_error(client_DSSP):
    """Test that DSSP raises clear error for insufficient residues."""
    u = mda.Universe(TPR, XTC)

    protein = u.select_atoms("protein")

    with pytest.raises(ValueError, match="DSSP requires at least 6 residues"):
        res2 = protein.residues[:2].atoms
        DSSP(res2)

    with pytest.raises(ValueError, match="DSSP requires at least 6 residues"):
        res4 = protein.residues[:4].atoms
        DSSP(res4)

    res6 = protein.residues[:6].atoms
    dssp = DSSP(res6)
    result = dssp.run(**client_DSSP, stop=1)
    assert result.results.dssp.shape[1] == 6


def _hydrogens_per_residue(ag, hydrogen_name):
    """hydrogens selected one residue at a time (the selection DSSP made
    before it was vectorised): the reference for the selection now"""
    return [
        res.atoms.select_atoms(f"name {hydrogen_name}")
        for res in ag.select_atoms("protein").residues
    ]


@pytest.mark.parametrize(
    "hydrogen_name", ["H HN HT1 HT2 HT3", "H HN H1", "H*"]
)
@pytest.mark.parametrize("selection", ["all", "protein and not name H HN"])
def test_hydrogens_match_per_residue_selection(hydrogen_name, selection):
    # without H and HN in the selection, the hydrogens are still looked up
    # among all atoms of each residue
    ag = mda.Universe(TPR, XTC).select_atoms(selection)
    dssp = DSSP(ag, hydrogen_name=hydrogen_name)
    expected = _hydrogens_per_residue(ag, hydrogen_name)
    np.testing.assert_equal(dssp._has_hydrogen, [bool(h) for h in expected])
    np.testing.assert_equal(
        dssp._first_hydrogens.ix, [h.ix[0] for h in expected if h]
    )


def test_coords_with_hydrogens_match_per_residue_lookup():
    # residues without a hydrogen (here the N-terminal MET, whose hydrogens
    # are H1-H3, and prolines) take the guessed hydrogen of the same index
    u = mda.Universe(TPR, XTC)
    coords = DSSP(u, guess_hydrogens=False)._get_coords()
    guessed = _get_hydrogen_atom_position(coords[:, :4])
    hydrogens = _hydrogens_per_residue(u.atoms, "H HN HT1 HT2 HT3")
    assert not all(hydrogens)
    expected = [
        h.positions[0] if h else guessed[k] for k, h in enumerate(hydrogens)
    ]
    np.testing.assert_equal(coords[:, 4], expected)
