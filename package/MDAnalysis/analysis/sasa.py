# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2025 The MDAnalysis Development Team and contributors
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
:Year: 2025
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

`Relative Solvent-Accessible Surface Area`_ (RSASA) or `Relative Solvent
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

For example, we can calculate the average area of adenylate kinase (AdK). The
structure is included within the test data files:

.. code-block:: python

    import MDAnalysis
    from MDAnalysis.transformations import unwrap
    from MDAnalysisTests.datafiles import TPR, XTC
    from MDAnalysis.analysis.sasa import SASA, RSASA

    # Load and unwrap the structure
    u = MDAnalysis.Universe(TPR, XTC)
    u.trajectory.add_transformations(unwrap(u.atoms))
    protein = u.select_atoms("protein")

    # Calculate area for protein
    r = SASA(protein).run()
    print("Total surface:", r.results.area.mean(axis=0).sum(), "Å²")

    # Load it into beta-factors property
    if not hasattr(u, 'tempfactors'):
        u.add_TopologyAttribute('tempfactors')
    protein.tempfactors = r.results.area.mean(axis=0)

    # Calculate area for each residue residues
    for res in protein.residues:
        print(res.resnum, res.tempfactors.sum())

    # Increase the accuracy of calculations
    r = SASA(protein, n_dots=1024).run()

    # Calculate protein's accessible surface area (no solvent probe)
    r = SASA(protein, probe_radius=0).run()

    # Calculate residue's relative exposed area
    r = RSASA(protein).run()
    relative_area = r.results.relative_area.mean(axis=0)
    for res, area in zip(protein.residues, relative_area):
        print(res.resname, res.resnum, f"{area:.2%}")

    # Calculate relative area only for side-chain atoms
    r = RSASA(protein, subsele="not backbone").run()

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
.. _`Relative Solvent-Accessible Surface Area`:
   https://en.wikipedia.org/wiki/Relative_accessible_surface_area
.. _`Relative Solvent Accessibility`:
   https://en.wikipedia.org/wiki/Relative_accessible_surface_area
"""

import numpy as np
import collections
import logging
import scipy.spatial

from ..due import due, Doi
from .base import AnalysisBase, ResultsGroup
from ..exceptions import NoDataError
from ..core import groups


logger = logging.getLogger("MDAnalysis.analysis.sasa")

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


def _fib_sphere(n_dots: int) -> np.ndarray:
    """Generate equidistant points on a sphere (Fibonacci sphere)"""
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


def _sasa(positions: np.ndarray, vdw_radii: np.ndarray,
          probe_radius: float, n_dots: int) -> np.ndarray:
    """Calculate Solvent-Accessible Surface
    Area for coordinates with given VdW radii"""
    # Increase radii by probe radius
    radii = vdw_radii.copy() + probe_radius
    max_radii = 2 * np.max(radii)

    # Find atom's neighbors using kd-tree
    kdt = scipy.spatial.KDTree(positions, 10)
    dots_available = set(range(n_dots))

    # Generate Fibonacci sphere
    fib_sphere = _fib_sphere(n_dots)

    dots = np.zeros(positions.shape[0], dtype=np.int64)
    for i in range(positions.shape[0]):
        # Scale fib. sphere and move it to the i-th atom's position
        sphere = np.array(fib_sphere, copy=True) * radii[i]
        sphere += positions[i]
        available = dots_available.copy()

        # Calculate kd-tree for sphere points
        kdt_sphere = scipy.spatial.KDTree(sphere, 10)

        # Iterate over neighbors of the i-th atom
        for j in kdt.query_ball_point(positions[i], max_radii):
            if j == i:
                continue
            if radii[i] < (radii[i] + radii[j]):
                available -= {
                    n for n in kdt_sphere.query_ball_point(
                        positions[j],
                        radii[j]
                    )
                }
        dots[i] = len(available)

    # Calculate fraction of accessible sphere surface
    return 4 * np.pi * radii ** 2 * (dots / n_dots)


class SASA(AnalysisBase):
    """
    Calculate Solvent-Accessible Surface Area for atoms in selection using
    `Shrake-Rupley`_ algorithm.

    Parameters
    ----------
    ag: :class:`AtomGroup`
        An MDAnalysis :class:`AtomGroup`. :class:`UpdatingAtomGroup` instances
        are not accepted.
    n_dots: int
        Resolution of the surface of each atom. A higher number of dots
        results in more precise measurements, but slows down the calculation.
        Defaults to 256.
    probe_radius: float
        Radius of the solvent probe. Defaults to 1.40, roughly the radius
        of a water molecule.

    Attributes
    ----------
    results.area: :class:`numpy.ndarray`
        Atom-wise Solvent-Accessible Surface Area.

    Examples
    --------
    Before analysis you need to unwrap the trajectory. If molecules are split
    over the periodic boundary calculated values will be incorrect:: 

        from MDAnalysisTests.datafiles import TPR, GRO
        from MDAnalysis.transformations import unwrap
        
        u = MDAnalysis.Universe(TPR, GRO)
        u.trajectory.add_transformations(unwrap(u.atoms))
        protein = u.select_atoms('protein')

    Analysis requires the universe to include atomic radii information. If 
    radii are missing, you can infer them from element information. In this 
    example, we also assign a default radius of 2.0 Å for unknown elements::
    
        from MDAnalysis.guesser.tables import vdwradii

        if not hasattr(u, 'radii'):
            u.add_TopologyAttr('radii')
            u.atoms.radii = [vdwradii.get(e, 2.0) for e in u.atoms.elements]

    For most systems, the default parameters are sufficient to run a SASA
    analysis::

        r = SASA(protein).run()
        
        print(f"Total area: {r.results.area.mean(axis=0).sum():.2f} Å²")
    
    To compute the *molecular* accessible surface area only, set the probe
    radius to zero::

        r = SASA(protein, probe_radius=0).run()

    You can improve the accuracy of the calculated area by increasing the 
    number of dots, but this may significantly slow down the calculation (or
    decrease it to speed up the calculation but reduce the accuracy)::

        r = SASA(protein, n_dots=1024).run()

    One neat trick is to store computed area in the `tempfactor` topology 
    attribute making it easy to aggregate values for example summing area by
    residues::

        r = SASA(protein).run()
        if not hasattr(u, 'tempfactors'):
            u.add_TopologyAttr('tempfactors')
        protein.tempfactors = r.results.area.mean(axis=0)
        
        for res in protein.residues:
            print(res.resnum, sum(res.atoms.tempfactors))
    
    .. _Shrake-Rupley::
        https://doi.org/10.1016/0022-2836(73)90011-9

    .. versionadded:: 2.11.0
    """
    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(self):
        return ("serial", "multiprocessing", "dask",)

    def __init__(self, ag, probe_radius=1.40, n_dots=256, **kwargs):
        super(SASA, self).__init__(ag.universe.trajectory, **kwargs)

        if isinstance(ag, groups.UpdatingAtomGroup):
            raise TypeError(
                "UpdatingAtomGroups are not valid for surface calculations")
        if probe_radius < 0.0:
            raise ValueError("Probe radius must be a positive number")
        if n_dots < 1:
            raise ValueError("Number of sphere dots must be larger than 1")
        if not hasattr(ag, "radii"):
            raise NoDataError("Universe has no 'radii' attribute")
        if np.any(ag.radii <= 0.0):
            raise ValueError("Atomic radii must be grater than zero")

        # Locals
        self.ag = ag
        self.probe_radius = float(probe_radius)
        self.n_dots = int(n_dots)

    def _prepare(self):
        self.results.area = []

    def _single_frame(self):
        self.results.area.append(
            _sasa(self.ag.positions, self.ag.radii,
                  self.probe_radius, self.n_dots)
        )

    def _get_aggregator(self):
        return ResultsGroup(lookup={"area": ResultsGroup.ndarray_vstack})

    def _conclude(self):
        self.results.area = np.array(self.results.area)


class RSASA(AnalysisBase):
    """
    Calculate Relative Solvent-Accessible Surface Area for residues in
    selection using `Shrake-Rupley`_ algorithm.

    Parameters
    ----------
    ag: :class:`AtomGroup`
        An MDAnalysis :class:`AtomGroup`. :class:`UpdatingAtomGroup` instances
        are not accepted.
    subsele: str
        Calculate surface only for sub-selection within the atomgroup e.g.,
        side-chain atoms. Defaults to `None`.
    n_dots: int
        Resolution of the surface of each atom. A higher number of dots
        results in more precise measurements, but slows down the calculation.
        Defaults to 256.
    probe_radius: float
        Radius of the solvent probe. Defaults to 1.40, roughly the radius of
        a water molecule.

    Attributes
    ----------
    results.relative_area: :class:`numpy.ndarray`
        Residue-wise Relative Solvent-Accessible Surface Area.
        Ranges from 0 to 1.

    Examples
    --------
    Before analysis you need to unwrap the trajectory. If molecules are split
    over the periodic boundary calculated values will be incorrect:: 

        from MDAnalysisTests.datafiles import TPR, GRO
        from MDAnalysis.transformations import unwrap
        
        u = MDAnalysis.Universe(TPR, GRO)
        u.trajectory.add_transformations(unwrap(u.atoms))
        protein = u.select_atoms('protein')

    Analysis requires the universe to include atomic radii information. If 
    radii are missing, you can infer them from element information. In this 
    example, we also assign a default radius of 2.0 Å for unknown elements::

        from MDAnalysis.guesser.tables import vdwradii

        if not hasattr(u, 'radii'):
            u.add_TopologyAttr('radii')
            u.atoms.radii = [vdwradii.get(e, 2.0) for e in u.atoms.elements]

    Additionally, analysis requires connectivity information i.e. `bonds`
    attribute. Use MDAnalysis built-in guesser if they are missing::

        u.guess_TopologyAttrs(to_guess=('bonds',))

    For most systems, the default parameters are sufficient to run a RSASA
    analysis::

        r = RSASA(protein).run()

        relative_area = r.results.relative_area.mean(axis=0)
        for res, area in zip(protein.residues, relative_area):
            print(res.resname, res.resnum, f"{area:.2%}")
    
    You can improve the accuracy of the calculated area by increasing the 
    number of dots, but this may significantly slow down the calculation (or
    decrease it to speed up the calculation but reduce the accuracy)::

        r = RSASA(protein, n_dots=1024).run()

    You can calculate RSASA for only sub-selection for, example for residue 
    side-chains only:
        
        r = RSASA(protein, subsele="not ").run()

    .. _`Shrake-Rupley`::
        https://doi.org/10.1016/0022-2836(73)90011-9

    .. versionadded:: 2.11.0
    """
    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(self):
        return ("serial", "multiprocessing", "dask",)
    
    def __init__(self, ag, subsele=None, probe_radius=1.40,
                 n_dots=256, **kwargs):
        super(RSASA, self).__init__(ag.universe.trajectory, **kwargs)

        if isinstance(ag, groups.UpdatingAtomGroup):
            raise TypeError(
                "UpdatingAtomGroups are not valid for surface calculations")
        if probe_radius < 0.0:
            raise ValueError("Probe radius must be a positive number")
        if n_dots < 1:
            raise ValueError("Number of sphere dots must be larger than 1")
        if not hasattr(ag, "radii"):
            raise NoDataError("Universe has no 'radii' attribute")
        if np.any(ag.radii <= 0.0):
            raise ValueError("Atomic radii must be grater than zero")
        if not hasattr(ag, "bonds"):
            raise NoDataError("Universe has no 'bonds' attribute")

        # Locals
        self.ag = ag
        self.probe_radius = float(probe_radius)
        self.n_dots = int(n_dots)
        self.subsele = subsele if subsele else "all"

    def _prepare(self):
        self.results.relative_area = []

    def _single_frame(self):
        # Calculate surface of (sub)selection and accumulate by residues
        sub = self.ag.select_atoms(self.subsele)
        area = _sasa(sub.positions, sub.radii, self.probe_radius, self.n_dots)
        result = collections.defaultdict(float)
        for i, atom in enumerate(sub.atoms):
            result[atom.resindex] += area[i]

        # Calculate surface of each isolated "tripeptide"
        for resindex in self.ag.residues.resindices:
            tripep = self.ag.select_atoms(
                f"(byres (bonded resindex {resindex})) and ({self.subsele})")
            if len(tripep) == 0:
                continue
            tripep_area = _sasa(tripep.positions, tripep.radii,
                                self.probe_radius, self.n_dots)
            # Sum only for relevant residue in tripeptide
            exposed_area = sum([
                a for a, i in zip(tripep_area, tripep.resindices)
                if i == resindex
            ])
            if exposed_area != 0.0:
                result[resindex] /= exposed_area

        # Update the result and account for residues that
        # might have empty selection
        self.results.relative_area.append(
            [result[i] for i in self.ag.residues.resindices])

    def _get_aggregator(self):
        return ResultsGroup(
            lookup={"relative_area": ResultsGroup.ndarray_vstack}
        )
    
    def _conclude(self):
        self.results.relative_area = np.array(self.results.relative_area)
