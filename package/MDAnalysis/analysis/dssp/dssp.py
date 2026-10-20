"""
Secondary structure assignment (helix, sheet and loop) --- :mod:`MDAnalysis.analysis.dssp`
==========================================================================================

:Author: Egor Marin
:Year: 2024
:Copyright: LGPL v2.1+

.. versionadded:: 2.8.0

The module contains code to build hydrogend bond contact map,
and use it to assign protein secondary structure (:class:`DSSP`).

This module uses the python version of the original algorithm :footcite:p:`Kabsch1983`,
re-implemented by @ShintaroMinami and available under MIT license from
`ShintaroMinami/PyDSSP <https://github.com/ShintaroMinami/PyDSSP/tree/master#differences-from-the-original-dssp>`_.


.. Note::
    This implementation does not discriminate different types of
    beta-sheets, as well as different types of helices, meaning you will get
    :math:`3_{10}` helices and π-helices labelled as "helix" too.


.. rubric:: Using original `pydssp`
The default implementation uses the original *pydssp* (v.0.9.0) code, rewritten
without usage of the *einops* library and hence having no dependencies. If you want
to explicitly use *pydssp* (or its particular version), install it to your
current environment with ``python -m pip install pydssp``. Please note that the
way MDAnalysis uses *pydssp* does not support *pydssp* 's capability for batch
processing or its use of the *pytorch* library.

When using this module in published work please cite :footcite:p:`Kabsch1983`.

.. rubric:: References

.. footbibliography::


Example applications
--------------------

Assigning secondary structure of a PDB file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example we will simply print a string representing a protein's secondary
structure.

.. code-block:: python

    from MDAnalysis.tests.datafiles import PDB
    from MDAnalysis.analysis.dssp import DSSP
    u = mda.Universe(PDB)
    s = ''.join(DSSP(u).run().results.dssp[0])
    print(s)


Calculating average secondary structure of a trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we take a trajectory and calculate its average secondary structure, i.e.
assign a secondary structure label 'X' to a residue if most of the frames in the
trajectory got assigned the 'X' label.

.. code-block:: python

    from MDAnalysis.analysis.dssp import DSSP, translate
    from MDAnalysisTests.datafiles import TPR, XTC
    u = mda.Universe(TPR, XTC)
    long_run = DSSP(u).run()
    mean_secondary_structure = translate(long_run.results.dssp_ndarray.mean(axis=0))
    print(''.join(mean_secondary_structure))

Running this code produces ::

   '--EEEE----------HHHHHHH----EE----HHHHH------HHHHHHHHHH------HHHHHHHHHHH---------EEEE-----HHHHHHHHH------EEEEEE--HHHHHH----EE--------EE---E----------------------HHHHHHHHHHHHHHHHHHHHHHHHHHHH----EEEEE------HHHHHHHHH--'

Find parts of the protein that maintain their secondary structure during simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this example, we will find residue groups that maintain their secondary structure
along the simulation, and have some meaningful ('E' or 'H') secondary structure
during more than set `threshold` fraction of frames. We will call these residues
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

Running this code produces ::

    '--EEEE----------HHHH'


Analysis classes
----------------

.. autoclass:: DSSP
   :members:
   :inherited-members:

   .. attribute:: results.dssp

      Contains the time series of the DSSP assignment as a
      :class:`numpy.ndarray` array of shape ``(n_frames, n_residues)`` where each row
      contains the assigned secondary structure character for each residue (whose
      corresponding resid is stored in :attr:`results.resids`). The three characters
      are ['H', 'E', '-'] and representi alpha-helix, sheet and loop, respectively.

   .. attribute:: results.dssp_ndarray

      Contains the one-hot encoding of the time series of the DSSP assignment
      as a :class:`numpy.ndarray` Boolean array of shape ``(n_frames, n_residues, 3)``
      where for each residue the encoding is stored as ``(3,)`` shape
      :class:`numpy.ndarray` of Booleans so that ``True`` at index 0 represents loop
      ('-'), ``True`` at index 1 represents helix ('H'), and ``True`` at index 2
      represents sheet 'E'.

      .. SeeAlso:: :func:`translate`


   .. attribute:: results.resids

      A :class:`numpy.ndarray` of length ``n_residues`` that contains the residue IDs
      (resids) for the protein residues that were assigned a secondary structure.


Functions
---------

.. autofunction:: assign
.. autofunction:: translate
"""

import warnings

import numpy as np

from MDAnalysis import AtomGroup, Universe

from ...due import Doi, due
from ..base import AnalysisBase, ResultsGroup

due.cite(
    Doi("10.1002/bip.360221211"),
    description="DSSP algorithm description",
    path="MDAnalysis.analysis.dssp",
    cite_module=True,
)

del Doi


try:  # pragma: no cover
    from pydssp.pydssp_numpy import (
        _get_hydrogen_atom_position,
        assign,
    )

    HAS_PYDSSP = True

except ModuleNotFoundError:
    HAS_PYDSSP = False
    from .pydssp_numpy import (
        _get_hydrogen_atom_position,
        assign,
    )


class DSSP(AnalysisBase):
    """Assign secondary structure using the DSSP algorithm.

    Analyze a selection containing a protein and assign secondary structure
    using the Kabsch-Sander algorithm :footcite:p:`Kabsch1983`. Only a subset
    of secondary structure categories are implemented:

    - 'H' represents a generic helix (α-helix, π-helix or :math:`3_{10}` helix)
    - 'E' represents 'extended strand', participating in beta-ladder (parallel
      or antiparallel)
    - '-' represents unordered part ("loop")

    The implementation was taken from the pydssp package (v. 0.9.0)
    https://github.com/ShintaroMinami/PyDSSP by Shintaro Minami under the
    MIT license.

    .. Warning::
       For DSSP to work properly, your atoms must represent a protein. The
       hydrogen atom bound to the backbone nitrogen atom is matched by name
       as given by the keyword argument `hydrogen_atom`. There may only be
       a single backbone nitrogen hydrogen atom per residue; the one exception
       is proline, for which there should not exist any such hydrogens.
       The default value of `hydrogen_atom` should handle the common naming
       conventions in the PDB and in force fields but if you encounter an error
       or unusual results during your run, try to figure out how to select the
       correct hydrogen atoms and report an issue in the MDAnalysis
       `issue tracker <https://github.com/MDAnalysis/mdanalysis/issues>`_.

    Parameters
    ----------
    atoms : Union[Universe, AtomGroup]
        input Universe or AtomGroup with at least 6 protein residues. In both cases, only protein residues will
        be chosen prior to the analysis via `select_atoms('protein')`.
        Heavy atoms of the protein are then selected by name
        `heavyatom_names`, and hydrogens are selected by name
        `hydrogen_name`.
    guess_hydrogens : bool, optional
        whether you want to guess hydrogens positions, by default ``True``.
        Guessing is made assuming perfect 120 degrees for all bonds that N
        atom makes, and a N-H bond length of 1.01 A.
        If ``guess_hydrogens`` is False, hydrogen atom positions on N atoms
        will be parsed from the trajectory, except for the "hydrogen" atom
        positions on PRO residues, and an N-terminal residue.
    heavyatom_names : tuple[str], default ("N", "CA", "C", "O O1 OT1")
        selection names that will be used to select "N", "CA", "C" and "O"
        atom coordinates for the secondary structure determination. The last
        string contains multiple values for "O" to account for C-term residues.
    hydrogen_name : str, default "H HN HT1 HT2 HT3"
        This selection should only select a single hydrogen atom in each residue
        (except proline), namely the one bound to the backbone nitrogen.

        .. Note::
           To work with different hydrogen-naming conventions by default, the
           default selection is broad but if hydrogens are incorrectly selected
           (e.g., a :exc:`ValueError` is raised) you must customize `hydrogen_name`
           for your specific case.

    backend : str, optional
        :ref:`Backend for distance calculations<selection-of-acceleration-backend>`, by default ``"serial"``.
        Can be set to ``"distopia"`` if `distopia`_ is installed.

        .. versionadded:: 2.11.0

    .. _distopia: https://www.mdanalysis.org/distopia/


    Raises
    ------
    ValueError
        If fewer than 6 residues are provided in the selection.
    ValueError
        if ``guess_hydrogens`` is True but some non-PRO hydrogens are missing.

    Examples
    --------

    For example, you can assign secondary structure for a single PDB file:

    >>> from MDAnalysis.analysis.dssp import DSSP
    >>> from MDAnalysisTests.datafiles import PDB
    >>> import MDAnalysis as mda
    >>> u = mda.Universe(PDB)
    >>> run = DSSP(u).run()
    >>> print("".join(run.results.dssp[0, :20]))
    --EEEEE-----HHHHHHHH

    The :attr:`results.dssp` holds the time series of assigned secondary
    structure, with one character for each residue.

    (Note that for displaying purposes we only print the first 20 residues
    of frame 0 with ``run.results.dssp[0, :20]`` but one would typically look
    at all residues ``run.results.dssp[0]``.)

    The :attr:`results.dssp_ndarray` attribute holds a
    ``(n_frames, n_residues, 3)`` shape ndarray with a *one-hot encoding*
    of *loop* '-' (index 0), *helix* 'H' (index 1), and *sheet* 'E'
    (index 2), respectively for each frame of the trajectory. It can be
    used to compute, for instance, the **average secondary structure**:

    >>> from MDAnalysis.analysis.dssp import translate, DSSP
    >>> from MDAnalysisTests.datafiles import TPR, XTC
    >>> u = mda.Universe(TPR, XTC)
    >>> run = DSSP(u).run()
    >>> mean_secondary_structure = translate(run.results.dssp_ndarray.mean(axis=0))
    >>> print(''.join(mean_secondary_structure)[:20])
    -EEEEEE------HHHHHHH


    .. versionadded:: 2.8.0

       Enabled **parallel execution** with the ``multiprocessing`` and ``dask``
       backends; use the new method :meth:`get_supported_backends` to see all
       supported backends.

    .. versionchanged:: 2.10.0
       Change treatment of proline and follow pydssp 0.9.1 (prolines are now explicitly
       forbidden to participate in the hydrogen bond network). Previous version could yield
       wrong assignment of prolines.
    """

    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(cls):
        return (
            "serial",
            "multiprocessing",
            "dask",
        )

    def __init__(
        self,
        atoms: Universe | AtomGroup,
        guess_hydrogens: bool = True,
        *,
        heavyatom_names: tuple[str] = ("N", "CA", "C", "O O1 OT1"),
        hydrogen_name: str = "H HN HT1 HT2 HT3",
        backend: str = "serial",
    ):
        self._guess_hydrogens = guess_hydrogens
        self._backend = backend

        ag: AtomGroup = atoms.select_atoms("protein")
        super().__init__(ag.universe.trajectory)

        # define necessary selections
        self._heavy_atoms: dict[str, AtomGroup] = {
            t: ag.atoms[
                np.isin(
                    ag.names, t.split()
                )  # need split() since `np.isin` takes an iterable as second argument
                # and "N".split() -> ["N"]
            ]
            for t in heavyatom_names
        }
        self._donor_mask: np.ndarray | None = ag.residues.resnames != "PRO"
        self._box = self._trajectory.ts.dimensions
        if self._box is not None and np.allclose(
            self._box, [1, 1, 1, 90, 90, 90]
        ):
            self._box = None
            warnings.warn(
                "Box dimensions are (1, 1, 1, 90, 90, 90), not using "
                "periodic boundary conditions in DSSP calculations"
            )
        # hydrogens of each residue (among all of its atoms, not only those in
        # ag), found with one selection: residue positions in ag.residues
        # (sorted by resindex), then the number per residue and the first one
        # of each residue that has any, in atom order
        residues = ag.residues
        hydrogens = residues.atoms.select_atoms(f"name {hydrogen_name}")
        h_residue = np.searchsorted(residues.resindices, hydrogens.resindices)
        n_hydrogens = np.bincount(h_residue, minlength=len(residues))
        self._has_hydrogen: np.ndarray = n_hydrogens > 0
        self._first_hydrogens: AtomGroup = hydrogens[
            np.unique(h_residue, return_index=True)[1]
        ]
        if not self._guess_hydrogens:
            # pairs the k-th CA atom with the k-th residue. This assumes that
            # both are ordered in the same way, which holds as long as ag is
            # sorted (it is, from the hard-coded protein selection above)
            calphas = self._heavy_atoms["CA"]
            n = min(len(calphas), len(residues))
            is_pro = calphas.resnames[1:n] == "PRO"
            counts = n_hydrogens[1:n]
            wrong = np.flatnonzero(np.where(is_pro, counts != 0, counts != 1))
            if wrong.size:
                residue = calphas[wrong[0] + 1].residue
                raise ValueError(
                    f"Residue {residue} contains wrong number of hydrogens: "
                    "exactly 1 hydrogen is expected for non-PRO residues, and "
                    "zero hydrogens for PRO residues."
                )

        positions = [group.positions for group in self._heavy_atoms.values()]
        if len({arr.shape[0] for arr in positions}) != 1:
            raise ValueError(
                "Universe contains unequal numbers of (N,CA,C,O) atoms ('name' field)."
                " Please select appropriate AtomGroup manually."
            )

        if len(ag.residues) < 6:
            raise ValueError(
                "DSSP requires at least 6 residues for secondary structure analysis, "
                f"but only {len(ag.residues)} residue(s) were provided in the selection."
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

        """
        # NOTE: here we explicitly rely on the fact that `self._heavy_atoms`
        # dictionary maintains order of the keys since python 3.7
        positions = [group.positions for group in self._heavy_atoms.values()]
        coords = np.array(positions)

        if not self._guess_hydrogens:
            guessed_h_coords = _get_hydrogen_atom_position(
                coords.swapaxes(0, 1)
            )

            # the first hydrogen of each residue that has one; a residue
            # without one takes the guessed hydrogen of the same index
            has = self._has_hydrogen
            h_coords = np.empty((len(has), 3), dtype=coords.dtype)
            h_coords[has] = self._first_hydrogens.positions
            missing = np.flatnonzero(~has)
            h_coords[missing] = guessed_h_coords[missing]
            coords = np.vstack([coords, h_coords[np.newaxis]])

        coords = coords.swapaxes(0, 1)
        return coords

    def _single_frame(self):
        coords = self._get_coords()
        dssp = assign(
            coords,
            donor_mask=self._donor_mask,
            box=self._box,
            backend=self._backend,
        )
        self.results.dssp_ndarray.append(dssp)

    def _conclude(self):
        self.results.dssp = translate(np.array(self.results.dssp_ndarray))
        self.results.dssp_ndarray = np.array(self.results.dssp_ndarray)
        self.results.resids = self._heavy_atoms["CA"].resids

    def _get_aggregator(self):
        return ResultsGroup(
            lookup={"dssp_ndarray": ResultsGroup.flatten_sequence},
        )


def translate(onehot: np.ndarray) -> np.ndarray:
    """Translate a one-hot encoding summary into char-based secondary structure
    assignment.

    One-hot encoding corresponds to C3 notation:
    '-', 'H', 'E' are loop, helix and sheet, respectively. Input array must
    have its last axis of shape 3: ``(n_residues, 3)`` or ``(n_frames, n_residues, 3)``

    Examples
    --------

    .. code-block:: python

        from MDAnalysis.analysis.dssp import translate
        import numpy as np
        # encoding 'HE-'
        onehot = np.array([[False, True, False],  # 'H'
                           [False, False, True],  # 'E'
                           [True, False, False]]) # '-'
        ''.join(translate(onehot))
        print(''.join(translate(onehot)))

    Running this code produces ::

        HE-

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
