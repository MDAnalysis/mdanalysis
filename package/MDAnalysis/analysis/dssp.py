"""
Secondary structure assignment (helix, sheet
and loop) --- :mod:`MDAnalysis.analysis.dssp`
==========================================================================

:Author: Egor Marin
:Year: 2024
:Copyright: GNU Public License v2

.. versionadded:: 2.8.0

The module contains code to build hydrogend bond contact map,
and use it to assign protein secondary structure (:class:`DSSP`).

This module uses the python version of the original algorithm by Kabsch & Sander (1983),
re-implemented by @ShintaroMinami and available under MIT license
[here](https://github.com/ShintaroMinami/PyDSSP/tree/master#differences-from-the-original-dssp)
Note that this implementation does not discriminate different types of beta-sheets,
as well as different types of helices, meaning you will get 3_10 helices
and pi-helices labelled as "helix" too.

When using this module in published work please cite [Kabsch1983]_.


Example applications
--------------------

Assigning secondary structure of a PDB file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this example we will simply print a string representing 
protein's secondary structure.

.. code-block:: python
    from MDAnalysis.tests.datafiles import PDB
    from MDAnalysis.analysis.dssp import DSSP
    u = mda.Universe(PDB)
    s = ''.join(DSSP(u).run().results.dssp)
    print(s)


Calculating average secondary structure of a trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here we take a trajectory and calculate its average secondary structure,
i.e. assign a secondary structure label 'X' to a residue if most of the frames
in the trajectory got assigned 'X' label.

.. code-block:: python
    from MDAnalysis.analysis.dssp import DSSP, translate
    from MDAnalysisTests.datafiles import TPR, XTC
    u = mda.Universe(TPR, XTC)
    long_run = DSSP(u).run(stop=20)
    mean_secondary_structure = translate(
        long_run.results.dssp_ndarray.mean(axis=0)
        )
    print(''.join(mean_secondary_structure)[:20])
    # '-EEEEEE------HHHHHHH'

Find parts of the protein that maintain their secondary structure during simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this example, we will find residue groups that maintain their secondary structure
along the simulation, and have some meaningful ('E' or 'H') secondary structure
during more than set `threshold` share of frames. We will call these residues
"persistent", for clarity, and label them according to the structure
that they maintain during the run:

.. code-block:: python
    from MDAnalysis.analysis.dssp import DSSP, translate
    from MDAnalysisTests.datafiles import TPR, XTC
    u = mda.Universe(TPR, XTC)
    threshold = 0.8

    long_run = DSSP(u).run()
    persistent_residues = translate(
        long_run
        .results
        .dssp_ndarray
        .mean(axis=0) > threshold
    )
    print(''.join(persistent_residues)[:20])
    # '--EEEE----------HHHH'


Functions
---------

.. autofunction:: get_hbond_map
.. autofunction:: assign
.. autofunction:: translate

Analysis classes
----------------

.. autoclass:: DSSP
   :members:
   :inherited-members:

   .. attribute:: results.dssp

       Contains the time series of the DSSP assignment
       as chars N×m :class:`numpy.ndarray` array with content
       ``[frame, residue]``, where `structure_type` is one of three characters:
       ['H', 'E', '-'], representing alpha-helix, sheet and loop, respectively.

   .. attribute:: results.dssp_ndarray

       Contains the one-hot encoding of the time series of the DSSP assignment
       as chars N×mx3 :class:`numpy.ndarray` array with content
       ``[frame, residue, encoding]``, where `encoding` is a (3,) shape
       :class:`numpy.ndarray` of booleans with axes representing loop '-',
       helix 'H' and sheet 'E', consequtively.
"""

from .base import AnalysisBase
import numpy as np
from .. import Universe

try:
    from einops import repeat, rearrange
except ModuleNotFoundError:
    HAS_EINOPS = False
else:
    HAS_EINOPS = True

CONST_Q1Q2 = 0.084
CONST_F = 332
DEFAULT_CUTOFF = -0.5
DEFAULT_MARGIN = 1.0


def _upsample(a: np.ndarray, window: int) -> np.ndarray:
    """Performs array upsampling with given window along given axis.
    
    Example
    -------
    .. code-block:: python
        hbmap = np.arange(4*4).reshape(4,4)
        print(hbmap)
        # [[ 0  1  2  3]
        #  [ 4  5  6  7]
        #  [ 8  9 10 11]
        #  [12 13 14 15]]
        
        print(_upsample(hbmap))
        # [[[[ 0  1  2]
        #    [ 4  5  6]
        #    [ 8  9 10]]

        #   [[ 1  2  3]
        #    [ 5  6  7]
        #    [ 9 10 11]]]


        #  [[[ 4  5  6]
        #    [ 8  9 10]
        #    [12 13 14]]

        #   [[ 5  6  7]
        #    [ 9 10 11]
        #    [13 14 15]]]]
        
    Parameters
    ----------
    a : np.ndarray
        input array
    window : int
        upsample window

    Returns
    -------
    np.ndarray
        unfolded array
    """
    return _unfold(_unfold(a, window, -2), window, -2)


def _unfold(a: np.ndarray, window: int, axis: int):
    "Helper function for 2D array upsampling"
    idx = np.arange(window)[:, None] + np.arange(a.shape[axis] - window + 1)[None, :]
    unfolded = np.take(a, idx, axis=axis)
    return  np.moveaxis(unfolded, axis-1, -1)


def _get_hydrogen_atom_position(coord: np.ndarray) -> np.ndarray:
    """Fills in hydrogen atoms positions if they are abscent, under the
    assumption that C-N-H and H-N-CA angles are perfect 120 degrees,
    and N-H bond is 1.01 A.

    Parameters
    ----------
    coord : np.ndarray
        input coordinates in Angstrom, shape (n_atoms, 4, 3),
        where second axes corresponds to (N, CA, C, O) atom coordinates

    Returns
    -------
    np.ndarray
        coordinates of additional hydrogens, shape (n_atoms-1, 3)

    .. versionadded:: 2.8.0
    """
    # C_i, N_i, H_i and CA_{i+1} are all in the peptide bond plane
    # we wanna get C_{i+1} - N_{i} vectors and normalize them
    # ---------
    # v1 = vec(C_i, N_i)
    # v2 = vec(CA_{i+1}, N_i)
    # v3 = vec(N_i, H_i) = ?
    # we use the assumption that all the angles are 120 degrees,
    # and |v3| = 1.01, hence
    # we can derive v3 = (v1/|v1| + v2/|v2|)*|v3|

    # get v1 = vec(C_i, N_i)
    vec_cn = coord[1:, 0] - coord[:-1, 2]
    vec_cn = vec_cn / np.linalg.norm(vec_cn, axis=-1, keepdims=True)

    # get v2 = vec(CA_{i+1}, N_{i})
    vec_can = coord[1:, 0] - coord[1:, 1]
    vec_can = vec_can / np.linalg.norm(vec_can, axis=-1, keepdims=True)

    vec_nh = vec_cn + vec_can
    vec_nh = vec_nh / np.linalg.norm(vec_nh, axis=-1, keepdims=True)

    # vec_(0, H) = vec(0, N) + vec_nh
    return coord[1:, 0] + 1.01 * vec_nh


def get_hbond_map(
    coord: np.ndarray,
    cutoff: float = DEFAULT_CUTOFF,
    margin: float = DEFAULT_MARGIN,
    return_e: bool = False,
) -> np.ndarray:
    """Returns hydrogen bond map

    Parameters
    ----------
    coord : np.ndarray
        input coordinates in either (n, 4, 3) or (n, 5, 3) shape 
        (without or with hydrogens). If hydrogens are not present, then
        ideal positions (see :func:_get_hydrogen_atom_positions) are used.
    cutoff : float, optional
        cutoff, by default DEFAULT_CUTOFF
    margin : float, optional
        margin, by default DEFAULT_MARGIN
    return_e : bool, optional
        if to return energy instead of hbond map, by default False

    Returns
    -------
    np.ndarray
        output hbond map or energy depending on return_e param

    Raises
    ------
    
    ImportError
        if module `einops` is not present

    .. versionadded:: 2.8.0
    """
    if not HAS_EINOPS:
        raise ImportError('DSSP: to use DSSP, please install einops')

    n_atoms, n_atom_types, _ = coord.shape
    assert n_atom_types in (4, 5), "Number of atoms should be 4 (N,CA,C,O) or 5 (N,CA,C,O,H)"

    if n_atom_types == 4:
        h_1 = _get_hydrogen_atom_position(coord)
    elif n_atom_types == 5:
        h_1 = coord[1:, 4]
        coord = coord[:, :4]
    else:
        raise ValueError("Number of atoms should be 4 (N,CA,C,O) or 5 (N,CA,C,O,H)")
    # after this:
    # h.shape == (n_atoms, 3)
    # coord.shape == (n_atoms, 4, 3)

    # distance matrix
    n_1, c_0, o_0 = coord[1:, 0], coord[0:-1, 2], coord[0:-1, 3]

    n = n_atoms - 1
    cmap = np.tile(c_0, (n, 1, 1))
    omap = np.tile(o_0, (n, 1, 1))
    nmap = np.tile(n_1, (1, 1, n)).reshape(n, n, 3)
    hmap = np.tile(h_1, (1, 1, n)).reshape(n, n, 3)

    d_on = np.linalg.norm(omap - nmap, axis=-1)
    d_ch = np.linalg.norm(cmap - hmap, axis=-1)
    d_oh = np.linalg.norm(omap - hmap, axis=-1)
    d_cn = np.linalg.norm(cmap - nmap, axis=-1)

    # electrostatic interaction energy
    # e[i, j] = e(CO_i) - e(NH_j)
    e = np.pad(
        CONST_Q1Q2 * (1.0 / d_on + 1.0 / d_ch - 1.0 / d_oh - 1.0 / d_cn) *
        CONST_F,
        [[1, 0], [0, 1]],
    )

    if return_e:
        return e

    # mask for local pairs (i,i), (i,i+1), (i,i+2)
    local_mask = ~np.eye(n_atoms, dtype=bool)
    local_mask *= ~np.diag(np.ones(n_atoms - 1, dtype=bool), k=-1)
    local_mask *= ~np.diag(np.ones(n_atoms - 2, dtype=bool), k=-2)
    # hydrogen bond map (continuous value extension of original definition)
    hbond_map = np.clip(cutoff - margin - e, a_min=-margin, a_max=margin)
    hbond_map = (np.sin(hbond_map / margin * np.pi / 2) + 1.0) / 2
    hbond_map = hbond_map * local_mask

    return hbond_map


def assign(coord: np.ndarray) -> np.ndarray:
    """Assigns secondary structure for a given coordinate array,
    either with or without assigned hydrogens

    Parameters
    ----------
    coord : np.ndarray
        input coordinates in either (n, 4, 3) or (n, 5, 3) shape,
        without or with hydrogens, respectively

    Returns
    -------
    np.ndarray
        output (n,) array with one-hot labels in C3 notation ('-', 'H', 'E'),
        representing loop, helix and sheet, respectively.

    Raises
    ------
    
        if module `einops` is not present

    .. versionadded:: 2.8.0
    """

    if not HAS_EINOPS:
        raise ImportError('DSSP: to use DSSP, please install einops')

    # get hydrogen bond map
    hbmap = get_hbond_map(coord)
    hbmap = np.swapaxes(hbmap, -1, -2) # convert into "i:C=O, j:N-H" form

    # identify turn 3, 4, 5
    turn3 = np.diagonal(hbmap, offset=3) > 0.0
    turn4 = np.diagonal(hbmap, offset=4) > 0.0
    turn5 = np.diagonal(hbmap, offset=5) > 0.0

    # assignment of helical secondary structures
    h3 = np.pad(turn3[:-1] * turn3[1:], [[1, 3]])
    h4 = np.pad(turn4[:-1] * turn4[1:], [[1, 4]])
    h5 = np.pad(turn5[:-1] * turn5[1:], [[1, 5]])

    # helix4 first, as alpha helix
    helix4 = h4 + np.roll(h4, 1, 0) + np.roll(h4, 2, 0) + np.roll(h4, 3, 0)
    h3 = h3 * ~np.roll(helix4, -1, 0) * ~helix4  # helix4 is higher prioritized
    h5 = h5 * ~np.roll(helix4, -1, 0) * ~helix4  # helix4 is higher prioritized
    helix3 = h3 + np.roll(h3, 1, 0) + np.roll(h3, 2, 0)
    helix5 = (h5 + np.roll(h5, 1, 0) + np.roll(h5, 2, 0) + np.roll(h5, 3, 0) +
              np.roll(h5, 4, 0))

    # identify bridge
    unfoldmap = _upsample(hbmap, 3) > 0.
    unfoldmap_rev = np.swapaxes(unfoldmap, 0, 1)

    p_bridge = (unfoldmap[:, :, 0, 1] * unfoldmap_rev[:, :, 1, 2]) + (
        unfoldmap_rev[:, :, 0, 1] * unfoldmap[:, :, 1, 2])
    p_bridge = np.pad(p_bridge, [[1, 1], [1, 1]])

    a_bridge = (unfoldmap[:, :, 1, 1] * unfoldmap_rev[:, :, 1, 1]) + (
        unfoldmap[:, :, 0, 2] * unfoldmap_rev[:, :, 0, 2])
    a_bridge = np.pad(a_bridge, [[1, 1], [1, 1]])

    # ladder
    ladder = (p_bridge + a_bridge).sum(-1) > 0.

    # H, E, L of C3
    helix = (helix3 + helix4 + helix5) > 0.
    strand = ladder
    loop = ~helix * ~strand

    onehot = np.stack([loop, helix, strand], axis=-1)
    return onehot


def translate(onehot: np.ndarray) -> np.ndarray:
    """Translate a one-hot encoding summary into char-based secondary
    structure assignment. One-hot encoding corresponds to C3 notation:
    '-', 'H', 'E' are loop, helix and sheet, respectively. Input array must
    have its last axis of shape 3: (n_residues, 3) or (n_frames, n_residues, 3)

    Examples
    --------

    >>> from MDAnalysis.analysis.dssp import translate
    >>> import numpy as np
    >>> # encoding 'HE-'
    >>> onehot = np.array([
        [False, True, False],  # 'H'
        [False, False, True],  # 'E'
        [True, False, False]]) # '-'
    >>> ''.join(translate(onehot))
    'HE-'

    Parameters
    ----------
    onehot : np.ndarray
        input array of one-hot encoding in ('-', 'H', 'E') order

    Returns
    -------
    np.ndarray
        array of '-', 'H' and 'E' symbols with secondary structure

    .. versionadded:: 2.8.0
    """
    C3_ALPHABET = np.array(["-", "H", "E"])
    index = np.argmax(onehot, axis=-1)
    return C3_ALPHABET[index]


class DSSP(AnalysisBase):
    """Assign secondary structure using DSSP algorithm as implemented
    by pydssp package (v. 0.9.0): https://github.com/ShintaroMinami/PyDSSP
    Here:
     - 'H' represents a generic helix (alpha-helix, pi-helix or 3-10 helix)
     - 'E' represents 'extended strand', participating in beta-ladder
     (parallel or antiparallel)
     - '-' represents unordered part

    Examples
    --------

    For example, you can assign secondary structure for a single PDB file:

    .. code-block:: python
        from MDAnalysis.analysis.dssp import DSSP
        from MDAnalysisTests.datafiles import PDB
        import MDAnalysis as mda
        u = mda.Universe(PDB)
        run = DSSP(u).run()
        print("".join(run.results.dssp[0][:20]))
        # '--EEEEE-----HHHHHHHH'

    Also, per-frame dssp assignment allows you to build average
    secondary structure -- `DSSP.results.dssp_ndarray` holds
    (n_frames, n_residues, 3) shape ndarray with one-hot encoding of
    loop, helix and sheet, respectively:

    .. code-block:: python
        from MDAnalysis.analysis.dssp import translate
        u = mda.Universe(...)
        long_run = DSSP(u).run()
        mean_secondary_structure = translate(
            long_run.results.dssp_ndarray.mean(axis=0)
            )
        print(''.join(mean_secondary_structure)[:20])
        # '---HHHHHHHH---------'

    .. versionadded:: 2.8.0
    """

    def __init__(self, u: Universe, guess_hydrogens: bool = True):
        if not HAS_EINOPS:
            raise ImportError('DSSP: to use DSSP, please install einops')
        super().__init__(u.trajectory)
        self._guess_hydrogens = guess_hydrogens

        # define necessary selections
        # O1 is for C-terminal residue
        heavyatom_names = ("N", "CA", "C", "O O1")
        self._heavy_atoms: dict[str, "AtomGroup"] = {
            t: u.select_atoms(f"protein and name {t}")
            for t in heavyatom_names
        }
        self._hydrogens: list["AtomGroup"] = [
            res.atoms.select_atoms("name H")
            for res in u.select_atoms("protein").residues
        ]
        # can't do it the other way because I need missing values to exist
        # so that I could fill them in later

    def _prepare(self):
        self.results.dssp_ndarray = []

    def _get_coords(self) -> np.ndarray:
        """Returns coordinates of (N,CA,C,O,H) atoms, as required by
        :func:`get_hbond_map` and :func:`assign` functions.

        Returns
        -------
        np.ndarray
            coordinates of (N,CA,C,O,H) atoms

        Raises
        ------
        ValueError
            if input Universe contains different number of (N,CA,C,O) atoms

        .. versionadded:: 2.8.0
        """
        positions = [group.positions for group in self._heavy_atoms.values()]
        if len(set(map(lambda arr: arr.shape[0], positions))) != 1:
            raise ValueError((
                "Universe contains not equal number of (N,CA,C,O) atoms ('name' field)."
                " Please select appropriate sub-universe manually."))
        coords = np.array(positions)

        if not self._guess_hydrogens:
            guessed_h_coords = _get_hydrogen_atom_position(
                coords.swapaxes(0, 1))
            h_coords = np.array([
                group.positions[0] if group else guessed_h_coords[idx]
                for idx, group in enumerate(self._hydrogens)
            ])
            h_coords = np.expand_dims(h_coords, axis=0)
            coords = np.vstack([coords, h_coords])

        coords = coords.swapaxes(0, 1)
        return coords

    def _single_frame(self):
        coords = self._get_coords()
        dssp = assign(coords)
        self.results.dssp_ndarray.append(dssp)

    def _conclude(self):
        self.results.dssp = translate(np.array(self.results.dssp_ndarray))
        self.results.dssp_ndarray = np.array(self.results.dssp_ndarray)
        self.results.resids = np.array(
            [at.resid for at in self._heavy_atoms["CA"]])
