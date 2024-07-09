# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2024 The MDAnalysis Development Team and contributors
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
r"""\
Solvent-Accessible Surface Area --- :mod:`MDAnalysis.analysis.sasa`
===================================================================

:Authors: Jure Cerar
:Year: 2024
:Copyright: GNU Public License v2

Solvent-Accessible Surface Area
-------------------------------

`Solvent-Accessible Surface Area`_ (SASA) or `Accessible Surface Area`_ (ASA)
is the surface area of molecule that is accessible to solvent molecules in a
given environment. In the context of biomolecules such as proteins or nucleic
acids, solvent-exposed surface area quantifies the extent to which specific
regions of the molecule are exposed and interact with the surrounding solvent
and other solutes.

Implementation here uses the Shrake-Rupley algorithm :footcite:p:`Shrake1973`.
Algorithm works by drawing a mesh of equidistant points around each particle.
The points are drawn at particle's van der Waals radius, :math:`r_{VDW}`, with
added radius of solvent probe, :math:`r_{probe}`, usually water molecule (a
typical value is 1.4 Å). All points are checked against the surface of
neighboring atoms to determine whether they are buried or exposed. The number
of points accessible is multiplied by the portion of surface area each point
represents to calculate the accessible surface, :math:`S`:

.. math::

    S = 4 * \pi * (r_{VDW} + r_{probe})^2 * N_{exposed} / N

Relative Solvent-Exposed Surface Area
-------------------------------------

`Relative Solvent-Exposed Surface Area`_ (RSASA) or `Relative Solvent
Accessibility`_ (RSA) is a normalized measure of the solvent accessibility of a
specific region within a molecule, usually a residue. It represents the ratio
of the actual exposed surface area to the maximum possible surface area
accessible, :math:`S_{max}`, in that particular region with value of 0 for
fully buried and 1 for fully exposed residue and is often expressed as a
percentage. By normalizing to the maximum value this allows for a more
standardized comparison between different molecules or regions within the same
molecule, facilitating the analysis of structural changes and functional
implications across various biological contexts. It is calculated as:

.. math::

    S_{relative} = S / S_{max}

To measure the relative solvent accessibility, one usually takes
:math:`S_{max}` values that have been obtained from `Gly-X-Gly` tripeptides,
where `X` is the residue of interest. Commonly pre-calculated :math:`S_{max}`
scales are used e.g. :footcite:p:`Tien2013`. Here a different approach is
implemented. The :math:`S_{max}` is calculated by removing all other residues
except its two next neighbors (if present) for a given structure with no need
for tabulated values. This approach has several advantages. It is more robust,
as it always guaranteed to return value between 0 and 1. It is not limited to
standard amino acid residues, and works for other polymers (DNA, RNA, etc.) as
well. Additionally, it is not dependant on atomic radii and probe radius used.
The results are however still comparable to values from literature.

Examples
--------

For example, we can find how the average area of adenylate kinase (AdK). The
structure is included within the test data files:

.. code-block:: python

    import MDAnalysis as mda
    from MDAnalysis.analysis.sasa import SASA, RSASA

    from MDAnalysisTests.datafiles import TPR, GRO
    from MDAnalysis.transformations import unwrap

    # Load and unwrap the structure
    u = MDAnalysis.Universe(TPR, GRO)
    u.trajectory.add_transformations(unwrap(u.atoms))
    protein = u.select_atoms("protein")

    # Calculate SASA
    S = SASA(protein).run()
    print(f"Total surface: {S.results.area.sum():.2f} Å^2")

    # Calculate SASA more accurately
    S = SASA(protein, n_dots=1024).run()
    print(f"Total surface: {S.results.area.sum():.2f} Å^2")

    # Calculate surface area (no probe)
    S = SASA(protein, probe_radius=0).run()
    print(f"Total surface: {S.results.area.sum():.2f} Å^2")

    # Calculate SASA with custom atom radii for H
    S = SASA(protein, radii_dict={"H": 1.2}).run()
    print(f"Total surface: {S.results.area.sum():.2f} Å^2")

    # Calculate relative SASA
    RS = RSASA(protein).run()
    for res, area in zip(protein.residues, RS.results.relative_area):
        print(res.resnum, res.resname[:3], f"{area:7.2%}")

    # Calculate relative SASA of side-chain atoms
    RS = RSASA(protein, subsele="not backbone").run()
    for res, area in zip(protein.residues, RS.results.relative_area):
        print(res.resnum, res.resname[:3], f"{area:7.2%}")

Notes
_____

There are several things that must be taken into account when calculating
surfaces, like using unwrapped trajectories, etc. When calculating SASA make
sure that the structure fluctuates around the equilibrium state;  If there are
major changes in structure use the `start`, `stop`,  and `step` keywords to
control which frames are used for calculations. Note that computation of SASA
can be intensive for large systems or long trajectories.

References
----------

.. footbibliography::


Classes
-------

.. autoclass:: SASA
.. autoclass:: RSASA


.. _`Solvent-Accessible Surface Area`:
   https://en.wikipedia.org/wiki/Accessible_surface_area
.. _`Accessible Surface Area`:
   https://en.wikipedia.org/wiki/Accessible_surface_area
.. _`Relative Solvent-Exposed Surface Area`:
   https://en.wikipedia.org/wiki/Relative_accessible_surface_area
.. _`Relative Solvent Accessibility`:
   https://en.wikipedia.org/wiki/Relative_accessible_surface_area
"""

import numpy as np
import collections
import logging
import scipy

from ..due import due, Doi
from .base import AnalysisBase
from ..topology.tables import vdwradii

from ..core import groups

logger = logging.getLogger('MDAnalysis.analysis.sasa')

due.cite(
    Doi("10.1016/0022-2836(73)90011-9"),
    description="Environment and exposure to solvent of protein atoms",
    path="MDAnalysis.analysis.sasa",
    cite_module=True,
)
due.cite(
    Doi("10.1371/journal.pone.0080635"),
    description="Max Allowed Solvent Accessibility of Residues in Proteins",
    path="MDAnalysis.analysis.sasa",
    cite_module=True,
)
due.cite(
    Doi("10.12688/f1000research.7931.1"),
    description="FreeSASA: An open source C library for SASA calculations",
    path="MDAnalysis.analysis.sasa",
    cite_module=True,
)
del Doi


class SASA(AnalysisBase):
    """
    Calculate Solvent-Accessible Surface Area for atoms in selection using the
    Shrake-Rupley algorithm.

    Parameters
    ----------
    ag : :class:`AtomGroup`
        An MDAnalysis :class:`AtomGroup`. :class:`UpdatingAtomGroup` instances
        are not accepted.
    n_dots : int
        Resolution of the surface of each atom. A higher number of points
        results in more precise measurements, but slows down the calculation.
        Defaults to 256.
    probe_radius : float
        Radius of the solvent probe in :math:`Angstroms`. Defaults to 1.40,
        roughly the radius of a water molecule.
    radii_dict : dict
        User-provided dictionary of atomic radii to use in the calculation.
        Values will replace or complement those in
        :data:`MDAnalysis.topology.tables.vdwradii`.For unknown particles radii
        defaults to 2.0 :math:`Angstroms`.

    Attributes
    ----------
    results.area : :class:`numpy.ndarray`
        Atom-wise Solvent-Accessible Surface Area in :math:`Angstrom^2`.
    radii : :class:`numpy.ndarray`
        Atomic radii (with probe) used in calculation.
    radii_dict : dict
        Dictionary of atomic radii used for assignment.
    probe_radius : float
        Radius of the solvent probe used in calculation.
    n_dots : int
        Resolution used for calculation.
    n_frames : int
        Number of frames included in calculation.

    Examples
    --------
    >>> R = SASA(u.atoms).run()
    >>> R = SASA(u.atoms, n_dots=256).run()
    >>> R = SASA(u.atoms, probe_radius=1.4).run()
    >>> R = SASA(u.atoms, radii_dict={"H": 1.2}).run()
    >>> R.results.area
    [ 2.999 0.614 19.942 ... 2.830 10.882 40.180 ]
    >>> R.results.area.sum()
    12146.27


    .. versionadded:: 2.8.0
    """

    def __init__(self, ag, probe_radius=1.40, n_dots=256,
                 radii_dict=None, **kwargs):
        """
        Parameters
        ----------
        ag : :class:`AtomGroup`
            An MDAnalysis :class:`AtomGroup`. :class:`UpdatingAtomGroup`
            instances are not accepted.
        n_dots : int
            Resolution of the surface of each atom. A higher number of points
            results in more precise measurements, but slows down the
            calculation. Defaults to 256.
        probe_radius : float
            Radius of the solvent probe in :math:`Angstroms`. Defaults to 1.40,
            roughly the radius of a water molecule.
        radii_dict : dict
            User-provided dictionary of atomic radii to use in the calculation.
            Values will replace or complement those in
            :data:`MDAnalysis.topology.tables.vdwradii`. For unknown particles
            radii defaults to 2.0 :math:`Angstroms`.
        """
        if isinstance(ag, groups.UpdatingAtomGroup):
            raise TypeError(
                f"UpdatingAtomGroups are not valid for surface calculations")

        super(SASA, self).__init__(ag.universe.trajectory, **kwargs)

        # Check input parameters and if AtomGroup has 'elements' property
        if probe_radius < 0.0:
            raise ValueError(
                f"Probe radius must be a positive number: {probe_radius} <= 0")
        if n_dots < 1:
            raise ValueError(
                f"Number of sphere dots must be larger than 1: {n_dots}")
        if not hasattr(ag, "elements"):
            raise ValueError(
                "Cannot assign atomic radii:"
                "Universe has no 'elements' property")

        # Locals
        self.ag = ag
        self.probe_radius = float(probe_radius)
        self.n_dots = int(n_dots)

        # Import internal VDW radii table and update with user values
        self.radii_dict = dict()
        self.radii_dict.update(vdwradii)
        if radii_dict is not None:
            self.radii_dict.update(radii_dict)

        # Assign atoms radii based on elements property
        self.radii = np.array(
            [self.radii_dict.get(e, 2.0) for e in self.ag.elements])
        self.radii += self.probe_radius
        self._max_radii = 2 * np.max(self.radii)

        # Issue a warning if any element is not in radii table.
        if not set(self.ag.elements).issubset(self.radii_dict.keys()):
            logger.warning(
                "Element could not be assigned a radius: Using default radius")

        # Pre-compute Fibonacci sphere
        self._sphere = self._get_sphere(self.n_dots)

    @staticmethod
    def _get_sphere(n_dots):
        """Generate sphere with equidistant points (Fibonacci sphere)"""
        dl = np.pi * (3 - np.sqrt(5))
        dz = 2.0 / n_dots
        longitude = 0
        z = 1 - dz / 2
        xyz = np.zeros((n_dots, 3), dtype=np.float32)
        for i in range(n_dots):
            r = np.sqrt(1 - z * z)
            xyz[i, 0] = np.cos(longitude) * r
            xyz[i, 1] = np.sin(longitude) * r
            xyz[i, 2] = z
            z -= dz
            longitude += dl
        return xyz

    def _prepare(self):
        self.results.area = np.zeros(self.ag.n_atoms)

    def _single_frame(self):
        # Find atom's neighbors using KDTree
        kdt = scipy.spatial.KDTree(self.ag.positions, 10)
        dots_available = set(range(self.n_dots))

        dots = np.zeros(self.ag.n_atoms)
        for i in range(self.ag.n_atoms):
            # Scale sphere and move it to the i-th atom position
            sphere = self._sphere.copy() * self.radii[i]
            sphere += self.ag.positions[i]
            available = dots_available.copy()
            kdt_sphere = scipy.spatial.KDTree(sphere, 10)

            # Iterate over neighbors of atom i
            for j in kdt.query_ball_point(self.ag.positions[i],
                                          self._max_radii, workers=-1):
                if j == i:
                    continue
                if self.radii[j] < (self.radii[i] + self.radii[j]):
                    available -= {
                        n for n in kdt_sphere.query_ball_point(
                            self.ag.positions[j],
                            self.radii[j]
                        )
                    }
            dots[i] = len(available)

        # Convert accessible points to surface area in A^2
        self.results.area += 4 * np.pi * self.radii ** 2 * dots / self.n_dots

    def _conclude(self):
        # Average for number of trajectory frames
        if self.n_frames != 0:
            self.results.area /= self.n_frames


class RSASA(AnalysisBase):
    """
    Calculate Relative Solvent-Accessible Surface Area for residues in
    selection using the Shrake-Rupley algorithm.

    Parameters
    ----------
    ag : :class:`AtomGroup`
        An MDAnalysis :class:`AtomGroup`. :class:`UpdatingAtomGroup` instances
        are not accepted.
    subsele : str
        Calculate surface only for sub-selection within the atomgroup e.g.
        side-chain atoms. Defaults to `None`.
    n_dots : int
        Resolution of the surface of each atom. A higher number of points
        results in more precise measurements, but slows down the calculation.
        Defaults to 256.
    probe_radius : float
        Radius of the solvent probe in :math:`Angstroms`. Defaults to 1.40,
        roughly the radius of a water molecule.
    radii_dict : dict
        User-provided dictionary of atomic radii to use in the calculation.
        Values will replace or complement those in
        :data:`MDAnalysis.topology.tables.vdwradii`. For unknown particles
        radii defaults to 2.0 :math:`Angstroms`.

    Attributes
    ----------
    results.relative_area : :class:`numpy.ndarray`
        Residue-wise Relative Solvent-Accessible Surface Area.
        Ranges from 0 to 1.
    radii : :class:`numpy.ndarray`
        Atomic radii (with probe) used in calculation.
    radii_dict : dict
        Dictionary of atomic radii used for assignment.
    probe_radius : float
        Radius of the solvent probe used in calculation.
    n_dots : int
        Resolution used for calculation.
    n_frames : int
        Number of frames included in calculation.

    Examples
    --------
    >>> R = RSASA(u.atoms).run()
    >>> R = RSASA(u.atoms, n_dots=256).run()
    >>> R = RSASA(u.atoms, probe_radius=1.4).run()
    >>> R = RSASA(u.atoms, radii_dict={"H": 1.2}).run()
    >>> R = RSASA(u.atoms, subsele="not backbone").run()
    >>> R.results.relative_area
    [ 0.215 0.232 0.002 0.000 ... 0.321 0.044 0.605 ]


    .. versionadded:: 2.8.0
    """

    def __init__(self, ag, subsele=None, probe_radius=1.40,
                 n_dots=256, radii_dict=None, **kwargs):
        """
        Parameters
        ----------
        ag : :class:`AtomGroup`
            An MDAnalysis :class:`AtomGroup`. :class:`UpdatingAtomGroup`
            instances are not accepted.
        subsele : str
            Calculate surface only for sub-selection within the atomgroup e.g.
            side-chain atoms. Defaults to `None`.
        n_dots : int
            Resolution of the surface of each atom. A higher number of points
            results in more precise measurements, but slows down the
            calculation. Defaults to 256.
        probe_radius : float
            Radius of the solvent probe in :math:`Angstroms`. Defaults to 1.40,
            roughly the radius of a water molecule.
        radii_dict : dict
            User-provided dictionary of atomic radii to use in the calculation.
            Values will replace or complement those in
            :data:`MDAnalysis.topology.tables.vdwradii`. For unknown particles
            radii defaults to 2.0 :math:`Angstroms`.
        """
        if isinstance(ag, groups.UpdatingAtomGroup):
            raise TypeError(
                f"UpdatingAtomGroups are not valid for surface calculations")

        super(RSASA, self).__init__(ag.universe.trajectory, **kwargs)

        # Check input parameters and if AtomGroup has
        # 'elements' and 'bonds' property
        if probe_radius < 0.0:
            raise ValueError(
                f"Probe radius must be a positive number: {probe_radius} <= 0")
        if n_dots < 1:
            raise ValueError(
                f"Number of sphere dots must be larger than 1: {n_dots}")
        if not hasattr(ag, "elements"):
            raise ValueError("Cannot assign atomic radii:"
                             "Universe has no 'elements' property")
        if not hasattr(ag, "bonds"):
            raise ValueError(
                "Universe has no 'bonds' property")

        # Locals
        self.ag = ag
        self._subsele = subsele if subsele else "all"
        self.probe_radius = float(probe_radius)
        self.n_dots = int(n_dots)

        # Import MDAnalysis VDW radii table and update with user values
        self.radii_dict = dict()
        self.radii_dict.update(vdwradii)
        if radii_dict is not None:
            self.radii_dict.update(radii_dict)

        # Assign atoms radii (for user to see) and issue a warning
        # if any element is not in radii table
        self.radii = np.array(
            [self.radii_dict.get(e, 2.0) for e in self.ag.elements])
        self.radii += self.probe_radius
        if not set(self.ag.elements).issubset(self.radii_dict.keys()):
            logger.warning(
                "Element could not be assigned a radius: Using default radius")

        # Pre-compute Fibonacci sphere
        self._sphere = self._get_sphere(self.n_dots)

    @staticmethod
    def _get_sphere(n_dots):
        """Generate sphere with equidistant points (Fibonacci sphere)"""
        dl = np.pi * (3 - np.sqrt(5))
        dz = 2.0 / n_dots
        longitude = 0
        z = 1 - dz / 2
        xyz = np.zeros((n_dots, 3), dtype=np.float32)
        for i in range(n_dots):
            r = np.sqrt(1 - z * z)
            xyz[i, 0] = np.cos(longitude) * r
            xyz[i, 1] = np.sin(longitude) * r
            xyz[i, 2] = z
            z -= dz
            longitude += dl
        return xyz

    def _get_sasa(self, ag):
        """Calculate SASA for given AtomGroup"""
        # Get radii for current AtomGroup
        radii = np.vectorize(self.radii_dict.get)(ag.elements)
        radii += self.probe_radius
        max_radii = 2 * np.max(radii)

        # Find atom's neighbors using KDTree
        kdt = scipy.spatial.KDTree(ag.positions, 10)
        dots_available = set(range(self.n_dots))

        dots = np.zeros(ag.n_atoms)
        for i in range(ag.n_atoms):
            # Scale sphere and move it to the i-th atom position
            sphere = self._sphere.copy() * radii[i]
            sphere += ag.positions[i]
            available = dots_available.copy()
            kdt_sphere = scipy.spatial.KDTree(sphere, 10)

            # Iterate over neighbors of atom i
            for j in kdt.query_ball_point(ag.positions[i],
                                          max_radii, workers=-1):
                if j == i:
                    continue
                if radii[j] < (radii[i] + radii[j]):
                    available -= {
                        n for n in kdt_sphere.query_ball_point(
                            ag.positions[j],
                            radii[j]
                        )
                    }
            dots[i] = len(available)

        # Convert accessible points to surface area in A^2
        return 4 * np.pi * radii ** 2 * dots / self.n_dots

    def _prepare(self):
        self.results.relative_area = np.zeros(self.ag.n_residues)

    def _single_frame(self):
        # Calculate surface of (sub)selection and accumulate by residues
        sub = self.ag.select_atoms(self._subsele)
        area = self._get_sasa(sub)
        result = collections.defaultdict(float)
        for i, atom in enumerate(sub.atoms):
            result[atom.resid] += area[i]

        # Calculate surface of each isolated "tripeptide"
        for resindex in self.ag.residues.resindices:
            tripep = self.ag.select_atoms(
                f"(byres (bonded resindex {resindex})) and ({self._subsele})")
            if len(tripep) == 0:
                continue
            tripep_area = self._get_sasa(tripep)
            exposed_area = sum([
                a for a, id in zip(tripep_area, tripep.resindices)
                if id == resindex
            ])
            if exposed_area != 0.0:
                result[resindex] /= exposed_area

        # Update the result and account for residues that
        # might have empty selection
        self.results.relative_area += np.array(
            [result[id] for id in self.ag.residues.resids])

    def _conclude(self):
        # Average for number of trajectory frames
        if self.n_frames != 0:
            self.results.relative_area /= self.n_frames
