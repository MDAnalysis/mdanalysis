from .base import AnalysisBase
import numpy as np
from einops import repeat, rearrange
from .. import Universe

CONST_Q1Q2 = 0.084
CONST_F = 332
DEFAULT_CUTOFF = -0.5
DEFAULT_MARGIN = 1.0

# code of the off-class functions
# is borrowed from https://github.com/ShintaroMinami/PyDSSP
# in accordance with the MIT License


def _unfold(a: np.ndarray, window: int, axis: int):
    """Unfolds an ndarray if it was batched (pydssp legacy)

    Args:
        a (np.ndarray): input array
        window (int): window
        axis (int): axis to unfold across

    Returns:
        np.ndarray: unfolded array
    """
    idx = np.arange(window)[:, None] + \
        np.arange(a.shape[axis] - window + 1)[None, :]
    unfolded = np.take(a, idx, axis=axis)
    return np.moveaxis(unfolded, axis - 1, -1)


def _check_input(coord):
    """Check if input coordinates have appropriate shape

    Args:
        coord (np.ndarray): input coordinates

    Returns:
        (np.ndarray, tuple): coordinates and appropriate shape
    """
    org_shape = coord.shape
    assert (len(org_shape) == 3) or (
        len(org_shape) == 4
    ), "Shape of input tensor should be [batch, L, atom, xyz] or [L, atom, xyz]"
    coord = repeat(coord, "... -> b ...",
                   b=1) if len(org_shape) == 3 else coord
    return coord, org_shape


def _get_hydrogen_atom_position(coord: np.ndarray) -> np.ndarray:
    """Fills in hydrogen atoms positions is they're abscent

    Args:
        coord (np.ndarray): input coordinates, shape (n_atoms, 4, 3)

    Returns:
        np.ndarray: coordinates of additional hydrogens, shape (n_atoms, 3)
    """
    # A little bit lazy (but should be OK) definition of H position here.
    vec_cn = coord[:, 1:, 0] - coord[:, :-1, 2]
    vec_cn = vec_cn / np.linalg.norm(vec_cn, axis=-1, keepdims=True)
    vec_can = coord[:, 1:, 0] - coord[:, 1:, 1]
    vec_can = vec_can / np.linalg.norm(vec_can, axis=-1, keepdims=True)
    vec_nh = vec_cn + vec_can
    vec_nh = vec_nh / np.linalg.norm(vec_nh, axis=-1, keepdims=True)
    return coord[:, 1:, 0] + 1.01 * vec_nh


def get_hbond_map(
    coord: np.ndarray,
    cutoff: float = DEFAULT_CUTOFF,
    margin: float = DEFAULT_MARGIN,
    return_e: bool = False,
) -> np.ndarray:
    """Returns hydrogen bond map

    Args:
        coord (np.ndarray): input coordinates in either (n, 4, 3) or (n, 5, 3) shape (without or with hydrogens)
        cutoff (float, optional): Defaults to DEFAULT_CUTOFF.
        margin (float, optional): margin. Defaults to DEFAULT_MARGIN.
        return_e (bool, optional): if to return energy instead of hbond map. Defaults to False.

    Returns:
        np.ndarray: output hbond map or energy depending on return_e param
    """
    # check input
    coord, org_shape = _check_input(coord)
    b, l, a, _ = coord.shape
    # add pseudo-H atom if not available
    assert (a == 4) or (
        a == 5
    ), "Number of atoms should be 4 (N,CA,C,O) or 5 (N,CA,C,O,H)"
    h = coord[:, 1:, 4] if a == 5 else _get_hydrogen_atom_position(coord)
    # distance matrix
    nmap = repeat(coord[:, 1:, 0], "... m c -> ... m n c", n=l - 1)
    hmap = repeat(h, "... m c -> ... m n c", n=l - 1)
    cmap = repeat(coord[:, 0:-1, 2], "... n c -> ... m n c", m=l - 1)
    omap = repeat(coord[:, 0:-1, 3], "... n c -> ... m n c", m=l - 1)
    d_on = np.linalg.norm(omap - nmap, axis=-1)
    d_ch = np.linalg.norm(cmap - hmap, axis=-1)
    d_oh = np.linalg.norm(omap - hmap, axis=-1)
    d_cn = np.linalg.norm(cmap - nmap, axis=-1)
    # electrostatic interaction energy
    e = np.pad(
        CONST_Q1Q2 * (1.0 / d_on + 1.0 / d_ch - 1.0 /
                      d_oh - 1.0 / d_cn) * CONST_F,
        [[0, 0], [1, 0], [0, 1]],
    )
    if return_e:
        return e
    # mask for local pairs (i,i), (i,i+1), (i,i+2)
    local_mask = ~np.eye(l, dtype=bool)
    local_mask *= ~np.diag(np.ones(l - 1, dtype=bool), k=-1)
    local_mask *= ~np.diag(np.ones(l - 2, dtype=bool), k=-2)
    # hydrogen bond map (continuous value extension of original definition)
    hbond_map = np.clip(cutoff - margin - e, a_min=-margin, a_max=margin)
    hbond_map = (np.sin(hbond_map / margin * np.pi / 2) + 1.0) / 2
    hbond_map = hbond_map * repeat(local_mask, "l1 l2 -> b l1 l2", b=b)
    # return h-bond map
    hbond_map = np.squeeze(hbond_map, axis=0) if len(
        org_shape) == 3 else hbond_map
    return hbond_map


def assign(coord: np.ndarray) -> np.ndarray:
    """Assigns secondary structure for a given coordinate array

    Args:
        coord (np.ndarray): input coordinates in either (n, 4, 3) or (n, 5, 3) shape -- without or with hydrogens

    Returns:
        np.ndarray: output (n,) array with labels in C3 notation ('-', 'H', 'E')
    """
    # check input
    coord, org_shape = _check_input(coord)
    # get hydrogen bond map
    hbmap = get_hbond_map(coord)
    hbmap = rearrange(
        hbmap, "... l1 l2 -> ... l2 l1"
    )  # convert into "i:C=O, j:N-H" form
    # identify turn 3, 4, 5
    turn3 = np.diagonal(hbmap, axis1=-2, axis2=-1, offset=3) > 0.0
    turn4 = np.diagonal(hbmap, axis1=-2, axis2=-1, offset=4) > 0.0
    turn5 = np.diagonal(hbmap, axis1=-2, axis2=-1, offset=5) > 0.0
    # assignment of helical sses
    h3 = np.pad(turn3[:, :-1] * turn3[:, 1:], [[0, 0], [1, 3]])
    h4 = np.pad(turn4[:, :-1] * turn4[:, 1:], [[0, 0], [1, 4]])
    h5 = np.pad(turn5[:, :-1] * turn5[:, 1:], [[0, 0], [1, 5]])
    # helix4 first
    helix4 = h4 + np.roll(h4, 1, 1) + np.roll(h4, 2, 1) + np.roll(h4, 3, 1)
    h3 = h3 * ~np.roll(helix4, -1, 1) * ~helix4  # helix4 is higher prioritized
    h5 = h5 * ~np.roll(helix4, -1, 1) * ~helix4  # helix4 is higher prioritized
    helix3 = h3 + np.roll(h3, 1, 1) + np.roll(h3, 2, 1)
    helix5 = (
        h5
        + np.roll(h5, 1, 1)
        + np.roll(h5, 2, 1)
        + np.roll(h5, 3, 1)
        + np.roll(h5, 4, 1)
    )
    # identify bridge
    unfoldmap = _unfold(_unfold(hbmap, 3, -2), 3, -2) > 0.0
    unfoldmap_rev = rearrange(unfoldmap, "b l1 l2 ... -> b l2 l1 ...")
    p_bridge = (unfoldmap[:, :, :, 0, 1] * unfoldmap_rev[:, :, :, 1, 2]) + (
        unfoldmap_rev[:, :, :, 0, 1] * unfoldmap[:, :, :, 1, 2]
    )
    p_bridge = np.pad(p_bridge, [[0, 0], [1, 1], [1, 1]])
    a_bridge = (unfoldmap[:, :, :, 1, 1] * unfoldmap_rev[:, :, :, 1, 1]) + (
        unfoldmap[:, :, :, 0, 2] * unfoldmap_rev[:, :, :, 0, 2]
    )
    a_bridge = np.pad(a_bridge, [[0, 0], [1, 1], [1, 1]])
    # ladder
    ladder = (p_bridge + a_bridge).sum(-1) > 0
    # H, E, L of C3
    helix = (helix3 + helix4 + helix5) > 0
    strand = ladder
    loop = ~helix * ~strand
    onehot = np.stack([loop, helix, strand], axis=-1)
    onehot = rearrange(
        onehot, "1 ... -> ...") if len(org_shape) == 3 else onehot
    return onehot


def translate(onehot: np.ndarray) -> np.ndarray:
    C3_ALPHABET = np.array(["-", "H", "E"])
    index = np.argmax(onehot, axis=-1)
    return C3_ALPHABET[index]


class DSSP(AnalysisBase):
    """Assign secondary structure using DSSP algorithm as implemented by pydssp package (v. 0.9.0)
    Here:
     - 'H' represents a 4-helix (alpha-helix)
     - 'E' represents 'extended strand', participating in beta-ladder
     - '-' represents unordered part

    Example:
    >>> from MDAnalysis.analysis.dssp import DSSP
    >>> from MDAnalysisTests.datafiles import PDB
    >>> import MDAnalysis as mda
    >>> u = mda.Universe(PDB)
    >>> run = DSSP(u).run()
    >>> print("".join(run.results.dssp[0]))
    '--EEEEE-----HHHHHHHHHHHH--EEE-HHHHHHHHHHH--HHHHHHHHHHHH-----HHHHHHHHHHH---HHH---EEEE-----HHHHHHHHHH-----EEEEEE--HHHHHHHH--EE--------EE---E------E------E----HHH-HHHHHHHHHHHHHHHHHHHHHHHHHHHH---EEEEEE----HHHHHHHHHHHH-'

    Also, per-frame dssp assignment allows you to build average secondary structure -- `DSSP.results.dssp_ndarray` 
    holds (n_frames, n_residues, 3) shape ndarray with one-hot encoding of loop, helix and sheet, respectively:

    >>> from MDAnalysis.analysis.dssp import translate
    >>> u = mda.Universe(...)
    >>> long_run = DSSP(u).run()
    >>> mean_secondary_structure = translate(long_run.results.dssp_ndarray.mean(axis=0))
    >>> ''.join(mean_secondary_structure)
    ---HHHHHHHH--------------------------HHHHHH-------HHHH---HHHHHHHHHHHHHH------------HHHHH------------------HHHHHHHH---
    """

    def __init__(self, u: Universe, guess_hydrogens: bool = True):
        super().__init__(u.trajectory)
        self._guess_hydrogens = guess_hydrogens

        # define necessary selections
        # O1 is for C-terminal residue
        heavyatom_names = ("N", "CA", "C", "O O1")
        self._heavy_atoms: dict[str, "AtomGroup"] = {
            t: u.select_atoms(f"protein and name {t}") for t in heavyatom_names
        }
        self._hydrogens: list["AtomGroup"] = [
            res.atoms.select_atoms("name H")
            for res in u.select_atoms("protein").residues
        ]  # can't do it the other way because I need missing values to exist so that I could fill them in later

    def _prepare(self):
        self.results.dssp_ndarray = [None for _ in range(self.n_frames)]

    def _single_frame(self):
        coords = np.array(
            [group.positions for group in self._heavy_atoms.values()])

        if not self._guess_hydrogens:
            guessed_h_coords = _get_hydrogen_atom_position(
                coords.swapaxes(0, 1))
            h_coords = np.array(
                [
                    group.positions[0] if group else guessed_h_coords[idx]
                    for idx, group in enumerate(self._hydrogens)
                ]
            )
            h_coords = np.expand_dims(h_coords, axis=0)
            coords = np.vstack([coords, h_coords])

        coords = coords.swapaxes(0, 1)
        dssp = assign(coords)
        self.results.dssp_ndarray[self._frame_index] = dssp

    def _conclude(self):
        self.results.dssp = translate(self.results.dssp_ndarray)
        self.results.dssp_ndarray = np.array(self.results.dssp_ndarray)
        self.results.resids = np.array(
            [at.resid for at in self._heavy_atoms["CA"]])


def report(arr: np.ndarray):
    print(f"{arr.shape=}")
