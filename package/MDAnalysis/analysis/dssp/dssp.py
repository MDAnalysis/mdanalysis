"""
Secondary structure assignment (helix, sheet and loop) --- :mod:`MDAnalysis.analysis.dssp`
=========================================================================================

:Author: Egor Marin
:Year: 2024
:Copyright: LGPL v2.1+

.. versionadded:: 2.8.0

The module contains code to build hydrogend bond contact map,
and use it to assign protein secondary structure (:class:`DSSP`).

This module uses the python version of the original algorithm by Kabsch & Sander (1983),
re-implemented by @ShintaroMinami and available under MIT license
from [ShintaroMinami/PyDSSP](https://github.com/ShintaroMinami/PyDSSP/tree/master#differences-from-the-original-dssp)
Note that this implementation does not discriminate different types of beta-sheets,
as well as different types of helices, meaning you will get 3_10 helices
and pi-helices labelled as "helix" too.


Using original `pydssp`
----------------------

The default implementation uses the original `pydssp` (v.0.9.0) code,
rewritten without usage of `einops` library and hence having no dependencies.
If you want to explicitly use `pydssp` (or its particular version), install it
to your current environment with `python3 -m pip install pydssp`. Please note
that the way MDAnalysis uses pydssp does not support pydssp's capability
for batch processing or its use of the pytorch library.

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
       helix 'H' and sheet 'E', consequently.
"""

import numpy as np
from MDAnalysis import Universe

from ..base import AnalysisBase
from ...due import due, Doi

due.cite(
    Doi("10.1002/bip.360221211"),
    description="DSSP algorithm description",
    path="MDAnalysis.analysis.dssp",
    cite_module=True,
)

del Doi


try:
    from pydssp.pydssp_numpy import (
        assign,
        get_hbond_map,
        _get_hydrogen_atom_position,
    )

    HAS_PYDSSP = True

except ModuleNotFoundError:
    HAS_PYDSSP = False
    from .pydssp_numpy import (
        assign,
        get_hbond_map,
        _get_hydrogen_atom_position,
    )


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
    """Assign secondary structure using DSSP algorithm.

    Analyze a selection containing a protein and assign secondary structure
    using the Kabsch-Sander algorithm [Kabsch1983]_. Only a subset of secondary
    structure categories are implemented:

    - 'H' represents a generic helix (alpha-helix, pi-helix or 3-10 helix)
    - 'E' represents 'extended strand', participating in beta-ladder
      (parallel or antiparallel)
    - '-' represents unordered part

    The implementation was taken from the pydssp package (v. 0.9.0)
    https://github.com/ShintaroMinami/PyDSSP by Shintaro Minami.

    Parameters
    ----------
    u : Universe
        input universe
    guess_hydrogens : bool, optional
        whether you want to guess hydrogens positions, by default True.
        Guessing is made assuming perfect 120 degrees for all bonds that N
        atom makes, and a N-H bond length of 1.01 A.
        If `guess_hydrogens` is False, hydrogen atom positions on N atoms
        will be parsed from the trajectory, except for the "hydrogen" atom
        positions on PRO residues, and an N-terminal residue.

    Raises
    ------
    ValueError
        if `guess_hydrogens` is True but some non-PRO hydrogens are missing.

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
        ]
        # can't do it the other way because I need missing values to exist
        # so that I could fill them in later

        if not self._guess_hydrogens:
            for calpha, hydrogen in zip(
                self._heavy_atoms["CA"][1:], self._hydrogens[1:]
            ):
                if not hydrogen and calpha.resname != "PRO":
                    raise ValueError(
                        (
                            "Universe is missing non-PRO hydrogen "
                            f"on residue {calpha.residue}"
                        )
                    )

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
            raise ValueError(
                (
                    "Universe contains not equal number of (N,CA,C,O) atoms ('name' field)."
                    " Please select appropriate sub-universe manually."
                )
            )
        coords = np.array(positions)

        if not self._guess_hydrogens:
            guessed_h_coords = _get_hydrogen_atom_position(coords.swapaxes(0, 1))

            h_coords = np.array(
                [
                    group.positions[0] if group else guessed_h_coords[idx]
                    for idx, group in enumerate(self._hydrogens)
                ]
            )
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
        self.results.resids = np.array([at.resid for at in self._heavy_atoms["CA"]])
