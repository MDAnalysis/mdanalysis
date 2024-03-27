# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
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

"""Hydrogen Bond Analysis --- :mod:`MDAnalysis.analysis.hydrogenbonds.hbond_analysis`
=====================================================================================

:Author: Paul Smith
:Year: 2019
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module provides methods to find and analyse hydrogen bonds in a Universe.

The :class:`HydrogenBondAnalysis` class is a new version of the original
:class:`MDAnalysis.analysis.hbonds.HydrogenBondAnalysis` class from the module
:mod:`MDAnalysis.analysis.hbonds.hbond_analysis`, which itself was modeled after the `VMD
HBONDS plugin`_.

.. _`VMD HBONDS plugin`: http://www.ks.uiuc.edu/Research/vmd/plugins/hbonds/


Input
------

Required:
  - *universe* : an MDAnalysis Universe object

Options:
  - *donors_sel* [None] : Atom selection for donors. If `None`, then will be identified via the topology.
  - *hydrogens_sel* [None] : Atom selection for hydrogens. If `None`, then will be identified via charge and mass.
  - *acceptors_sel* [None] : Atom selection for acceptors. If `None`, then will be identified via charge.
  - *d_h_cutoff* (Å) [1.2] : Distance cutoff used for finding donor-hydrogen pairs
  - *d_a_cutoff* (Å) [3.0] : Distance cutoff for hydrogen bonds. This cutoff refers to the D-A distance.
  - *d_h_a_angle_cutoff* (degrees) [150] : D-H-A angle cutoff for hydrogen bonds.
  - *update_selections* [True] : If true, will update atom selections at each frame.


Output
------

- *frame* : frame at which a hydrogen bond was found
- *donor id* : atom id of the hydrogen bond donor atom
- *hydrogen id* : atom id of the hydrogen bond hydrogen atom
- *acceptor id* : atom id of the hydrogen bond acceptor atom
- *distance* (Å): length of the hydrogen bond
- *angle* (degrees): angle of the hydrogen bond

Hydrogen bond data are returned in a :class:`numpy.ndarray` on a "one line, one observation" basis
and can be accessed via :attr:`HydrogenBondAnalysis.results.hbonds`::

    results = [
        [
            <frame>,
            <donor index (0-based)>,
            <hydrogen index (0-based)>,
            <acceptor index (0-based)>,
            <distance>,
            <angle>
        ],
        ...
    ]

Example use of :class:`HydrogenBondAnalysis`
--------------------------------------------

The simplest use case is to allow :class:`HydrogenBondAnalysis` to guess the acceptor and hydrogen atoms, and to
identify donor-hydrogen pairs via the bonding information in the topology::

  import MDAnalysis
  from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

  u = MDAnalysis.Universe(psf, trajectory)

  hbonds = HBA(universe=u)
  hbonds.run()

It is also possible to specify which hydrogens and acceptors to use in the analysis. For example, to find all hydrogen
bonds in water::

  import MDAnalysis
  from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

  u = MDAnalysis.Universe(psf, trajectory)

  hbonds = HBA(universe=u, hydrogens_sel='resname TIP3 and name H1 H2', acceptors_sel='resname TIP3 and name OH2')
  hbonds.run()

Alternatively, :attr:`hydrogens_sel` and :attr:`acceptors_sel` may be generated via the :attr:`guess_hydrogens` and
:attr:`guess_acceptors`. This selection strings may then be modified prior to calling :attr:`run`, or a subset of
the universe may be used to guess the atoms. For example, find hydrogens and acceptors belonging to a protein::

  import MDAnalysis
  from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

  u = MDAnalysis.Universe(psf, trajectory)

  hbonds = HBA(universe=u)
  hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
  hbonds.acceptors_sel = hbonds.guess_acceptors("protein")
  hbonds.run()

Slightly more complex selection strings are also possible. For example, to find hydrogen bonds involving a protein and
any water molecules within 10 Å of the protein (which may be useful for subsequently finding the lifetime of
protein-water hydrogen bonds or finding water-bridging hydrogen bond paths)::

  import MDAnalysis
  from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

  u = MDAnalysis.Universe(psf, trajectory)

  hbonds = HBA(universe=u)

  protein_hydrogens_sel = hbonds.guess_hydrogens("protein")
  protein_acceptors_sel = hbonds.guess_acceptors("protein")

  water_hydrogens_sel = "resname TIP3 and name H1 H2"
  water_acceptors_sel = "resname TIP3 and name OH2"

  hbonds.hydrogens_sel = f"({protein_hydrogens_sel}) or ({water_hydrogens_sel} and around 10 not resname TIP3})"
  hbonds.acceptors_sel = f"({protein_acceptors_sel}) or ({water_acceptors_sel} and around 10 not resname TIP3})"
  hbonds.run()

To calculate the hydrogen bonds between different groups, for example a
protein and water, one can use the :attr:`between` keyword. The
following will find protein-water hydrogen bonds but not protein-protein
or water-water hydrogen bonds::

  import MDAnalysis
  from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
    HydrogenBondAnalysis as HBA)

  u = MDAnalysis.Universe(psf, trajectory)

  hbonds = HBA(
    universe=u,
    between=['resname TIP3', 'protein']
    )

  protein_hydrogens_sel = hbonds.guess_hydrogens("protein")
  protein_acceptors_sel = hbonds.guess_acceptors("protein")

  water_hydrogens_sel = "resname TIP3 and name H1 H2"
  water_acceptors_sel = "resname TIP3 and name OH2"

  hbonds.hydrogens_sel = f"({protein_hydrogens_sel}) or ({water_hydrogens_sel})"
  hbonds.acceptors_sel = f"({protein_acceptors_sel}) or ({water_acceptors_sel})"

  hbonds.run()

It is further possible to compute hydrogen bonds between several groups with
with use of :attr:`between`. If in the above example,
`between=[['resname TIP3', 'protein'], ['protein', 'protein']]`, all
protein-water and protein-protein hydrogen bonds will be found, but
no water-water hydrogen bonds.

One can also define hydrogen bonds with atom types::

  from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
  hbonds = HBA(
               universe=u,
               donors_sel='type 2',
               hydrogens_sel='type 1',
               acceptors_sel='type 2',
              )

In order to compute the hydrogen bond lifetime, after finding hydrogen bonds
one can use the :attr:`lifetime` function::

    ...
    hbonds.run()
    tau_timeseries, timeseries = hbonds.lifetime()

It is **highly recommended** that a topology with bond information is used to
generate the universe, e.g `PSF`, `TPR`, or `PRMTOP` files. This is the only
method by which it can be guaranteed that donor-hydrogen pairs are correctly
identified. However, if, for example, a `PDB` file is used instead, a
:attr:`donors_sel` may be provided along with a :attr:`hydrogens_sel` and the
donor-hydrogen pairs will be identified via a distance cutoff,
:attr:`d_h_cutoff`::

  import MDAnalysis
  from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
    HydrogenBondAnalysis as HBA)

  u = MDAnalysis.Universe(pdb, trajectory)

  hbonds = HBA(
    universe=u,
    donors_sel='resname TIP3 and name OH2',
    hydrogens_sel='resname TIP3 and name H1 H2',
    acceptors_sel='resname TIP3 and name OH2',
    d_h_cutoff=1.2
  )
  hbonds.run()

The class and its methods
-------------------------

.. autoclass:: HydrogenBondAnalysis
   :members:

   .. attribute:: results.hbonds

      A :class:`numpy.ndarray` which contains a list of all observed hydrogen
      bond interactions. See `Output`_ for more information.

      .. versionadded:: 2.0.0

   .. attribute:: hbonds

      Alias to the :attr:`results.hbonds` attribute.

      .. deprecated:: 2.0.0
         Will be removed in MDAnalysis 3.0.0. Please use
         :attr:`results.hbonds` instead.
"""
import logging
import warnings
from collections.abc import Iterable

import numpy as np

from ..base import AnalysisBase, Results
from MDAnalysis.lib.distances import capped_distance, calc_angles
from MDAnalysis.lib.correlations import autocorrelation, correct_intermittency
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.analysis.hydrogenbonds.hbond_autocorrel import find_hydrogen_donors

from ...due import due, Doi


logger = logging.getLogger(__name__)


due.cite(Doi("10.1039/C9CP01532A"),
         description="Hydrogen bond analysis implementation",
         path="MDAnalysis.analysis.hydrogenbonds.hbond_analysis",
         cite_module=True)

del Doi


class HydrogenBondAnalysis(AnalysisBase):
    """
    Perform an analysis of hydrogen bonds in a Universe.
    """

    def __init__(self, universe,
                 donors_sel=None, hydrogens_sel=None, acceptors_sel=None,
                 between=None, d_h_cutoff=1.2,
                 d_a_cutoff=3.0, d_h_a_angle_cutoff=150,
                 update_selections=True):
        """Set up atom selections and geometric criteria for finding hydrogen
        bonds in a Universe.

        Hydrogen bond selections with `donors_sel` , `hydrogens_sel`, and
        `acceptors_sel` may be achieved with either a *resname*, atom *name* 
        combination, or when those are absent, with atom *type* selections.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        donors_sel : str
            Selection string for the hydrogen bond donor atoms. If the
            universe topology contains bonding information, leave
            :attr:`donors_sel` as `None` so that donor-hydrogen pairs can be
            correctly identified.
        hydrogens_sel :  str
            Selection string for the hydrogen bond hydrogen atoms. Leave as
            `None` to guess which hydrogens to use in the analysis using
            :attr:`guess_hydrogens`. If :attr:`hydrogens_sel` is left as
            `None`, also leave :attr:`donors_sel` as None so that
            donor-hydrogen pairs can be correctly identified.
        acceptors_sel : str
            Selection string for the hydrogen bond acceptor atoms. Leave as
            `None` to guess which atoms to use in the analysis using
            :attr:`guess_acceptors`
        between : List (optional),
            Specify two selection strings for non-updating atom groups between
            which hydrogen bonds will be calculated. For example, if the donor
            and acceptor selections include both protein and water, it is
            possible to find only protein-water hydrogen bonds - and not
            protein-protein or water-water - by specifying
            between=["protein", "SOL"]`. If a two-dimensional list is
            passed, hydrogen bonds between each pair will be found. For
            example, between=[["protein", "SOL"], ["protein", "protein"]]`
            will calculate all protein-water and protein-protein hydrogen
            bonds but not water-water hydrogen bonds. If `None`, hydrogen
            bonds between all donors and acceptors will be calculated.
        d_h_cutoff : float (optional)
            Distance cutoff used for finding donor-hydrogen pairs.
            Only used to find donor-hydrogen pairs if the
            universe topology does not contain bonding information
        d_a_cutoff : float (optional)
            Distance cutoff for hydrogen bonds. This cutoff refers to the D-A distance.
        d_h_a_angle_cutoff : float (optional)
            D-H-A angle cutoff for hydrogen bonds, in degrees.
        update_selections : bool (optional)
            Whether or not to update the acceptor, donor and hydrogen
            lists at each frame.

        Note
        ----

        It is highly recommended that a universe topology with bond
        information is used, as this is the only way that guarantees the
        correct identification of donor-hydrogen pairs.


        .. versionadded:: 2.0.0
            Added `between` keyword
        .. versionchanged:: 2.4.0
            Added use of atom types in selection strings for hydrogen atoms, 
            bond donors, or bond acceptors

        """

        self.u = universe
        self._trajectory = self.u.trajectory

        self.donors_sel = donors_sel.strip() if donors_sel is not None else donors_sel
        self.hydrogens_sel = hydrogens_sel.strip() if hydrogens_sel is not None else hydrogens_sel
        self.acceptors_sel = acceptors_sel.strip() if acceptors_sel is not None else acceptors_sel

        msg = ("{} is an empty selection string - no hydrogen bonds will "
               "be found. This may be intended, but please check your "
               "selection."
               )
        for sel in ['donors_sel', 'hydrogens_sel', 'acceptors_sel']:
            val = getattr(self, sel)
            if isinstance(val, str) and not val:
                warnings.warn(msg.format(sel))

        # If hydrogen bonding groups are selected, then generate
        # corresponding atom groups
        if between is not None:
            if not isinstance(between, Iterable) or len(between) == 0:
                raise ValueError("between must be a non-empty list/iterable")
            if isinstance(between[0], str):
                between = [between]

            between_ags = []
            for group1, group2 in between:
                between_ags.append(
                    [
                        self.u.select_atoms(group1, updating=False),
                        self.u.select_atoms(group2, updating=False)
                    ]
                )

            self.between_ags = between_ags
        else:
            self.between_ags = None


        self.d_h_cutoff = d_h_cutoff
        self.d_a_cutoff = d_a_cutoff
        self.d_h_a_angle = d_h_a_angle_cutoff
        self.update_selections = update_selections
        self.results = Results()
        self.results.hbonds = None

    def guess_hydrogens(self,
                        select='all',
                        max_mass=1.1,
                        min_charge=0.3,
                        min_mass=0.9
                        ):
        """Guesses which hydrogen atoms should be used in the analysis.

        Parameters
        ----------
        select: str (optional)
            :ref:`Selection string <selection-commands-label>` for atom group 
            from which hydrogens will be identified. (e.g., ``(resname X and 
            name H1)`` or ``type 2``)
        max_mass: float (optional)
            The mass of a hydrogen atom must be less than this value.
        min_mass: float (optional)
            The mass of a hydrogen atom must be greater than this value.
        min_charge: float (optional)
            The charge of a hydrogen atom must be greater than this value.

        Returns
        -------
        potential_hydrogens: AtomGroup
            AtomGroup corresponding to all hydrogen atoms potentially capable 
            of forming hydrogen bonds.

        Notes
        -----
        Hydrogen selections may be achieved with either a resname, atom 
        name combination, or when those are absent, atom types.

        This function makes use of atomic masses and atomic charges to identify
        which atoms are hydrogen atoms that are capable of participating in 
        hydrogen bonding. If an atom has a mass less than :attr:`max_mass` and 
        an atomic charge greater than :attr:`min_charge` then it is considered 
        capable of participating in hydrogen bonds.

        If :attr:`hydrogens_sel` is `None`, this function is called to guess 
        the selection.

        Alternatively, this function may be used to quickly generate a 
        :class:`AtomGroup` of potential hydrogen atoms involved in hydrogen 
        bonding. This :class:`AtomGroup` may then be modified before being used
        to set the attribute :attr:`hydrogens_sel`.


        .. versionchanged:: 2.4.0
            Added ability to use atom types

        """

        if min_mass > max_mass:
            raise ValueError("min_mass is higher than (or equal to) max_mass")

        ag = self.u.select_atoms(select)
        hydrogens_ag = ag[
            np.logical_and.reduce((
                ag.masses < max_mass,
                ag.charges > min_charge,
                ag.masses > min_mass,
            ))
        ]

        return hydrogens_ag

    def guess_donors(self, select='all', max_charge=-0.5):
        """Guesses which atoms could be considered donors in the analysis. Only
        use if the universe topology does not contain bonding information, 
        otherwise donor-hydrogen pairs may be incorrectly assigned.

        Parameters
        ----------
        select: str (optional)
            :ref:`Selection string <selection-commands-label>` for atom group 
            from which donors will be identified. (e.g., ``(resname X and name 
            O1)`` or ``type 2``)
        max_charge: float (optional)
            The charge of a donor atom must be less than this value.

        Returns
        -------
        potential_donors: AtomGroup
            AtomGroup corresponding to all atoms potentially capable of forming
            hydrogen bonds.

        Notes
        -----
        Donor selections may be achieved with either a resname, atom 
        name combination, or when those are absent, atom types.

        This function makes use of and atomic charges to identify which atoms 
        could be considered donor atoms in the hydrogen bond analysis. If an 
        atom has an atomic charge less than :attr:`max_charge`, and it is 
        within :attr:`d_h_cutoff` of a hydrogen atom, then it is considered 
        capable of participating in hydrogen bonds.

        If :attr:`donors_sel` is `None`, and the universe topology does not 
        have bonding information, this function is called to guess the 
        selection.

        Alternatively, this function may be used to quickly generate a 
        :class:`AtomGroup` of potential donor atoms involved in hydrogen 
        bonding. This :class:`AtomGroup` may then be modified before being used
        to set the attribute :attr:`donors_sel`.


        .. versionchanged:: 2.4.0
            Added ability to use atom types

        """

        # We need to know `hydrogens_sel` before we can find donors
        # Use a new variable `hydrogens_sel` so that we do not set 
        # `self.hydrogens_sel` if it is currently `None`
        hydrogens_ag = self.guess_hydrogens()
        if self.hydrogens_sel is None:
            hydrogens_sel = self._group_categories(hydrogens_ag)
        else:
            hydrogens_sel = self.hydrogens_sel

        # We're using u._topology.bonds rather than u.bonds as it is a million
        # times faster to access. This is because u.bonds also calculates
        # properties of each bond (e.g bond length). See:
        # https://github.com/MDAnalysis/mdanalysis/issues/2396#issuecomment-596251787
        if (hasattr(self.u._topology, 'bonds') 
           and len(self.u._topology.bonds.values) != 0):
            donors_ag = find_hydrogen_donors(hydrogens_ag)
            donors_ag = donors_ag.intersection(self.u.select_atoms(select))
        else:
            donors_ag = hydrogens_ag.residues.atoms.select_atoms(
                "({donors_sel}) and around {d_h_cutoff} {hydrogens_sel}".format(
                    donors_sel=select,
                    d_h_cutoff=self.d_h_cutoff,
                    hydrogens_sel=hydrogens_sel
                )
            )

        return donors_ag[donors_ag.charges < max_charge]

    def guess_acceptors(self, select='all', max_charge=-0.5):
        """Guesses which atoms could be considered acceptors in the analysis.

        Acceptor selections may be achieved with either a resname, atom 
        name combination, or when those are absent, atom types.

        Parameters
        ----------
        select: str (optional)
            :ref:`Selection string <selection-commands-label>` for atom group
            from which acceptors will be identified. (e.g., ``(resname X and 
            name O1)`` or ``type 2``)
        max_charge: float (optional)
            The charge of an acceptor atom must be less than this value.

        Returns
        -------
        potential_acceptors: AtomGroup
            AtomGroup corresponding to all atoms potentially capable of forming
            hydrogen bonds.

        Notes
        -----
        Acceptor selections may be achieved with either a resname, atom 
        name combination, or when those are absent, atom types.

        This function makes use of and atomic charges to identify which atoms 
        could be considered acceptor atoms in the hydrogen bond analysis. If 
        an atom has an atomic charge less than :attr:`max_charge` then it is 
        considered capable of participating in hydrogen bonds.

        If :attr:`acceptors_sel` is `None`, this function is called to guess 
        the selection.

        Alternatively, this function may be used to quickly generate a 
        :class:`AtomGroup` of potential acceptor atoms involved in hydrogen 
        bonding. This :class:`AtomGroup` may then be modified before being used
        to set the attribute :attr:`acceptors_sel`.


        .. versionchanged:: 2.4.0
            Added ability to use atom types

        """

        ag = self.u.select_atoms(select, updating=self.update_selections)

        return ag[ag.charges < max_charge]

    @staticmethod
    def _group_categories(group):
        """ Find categories according to universe constraints
        
        Parameters
        ----------
        group : AtomGroup
            AtomGroups corresponding to either hydrogen bond acceptors, 
            donors, or hydrogen atoms that meet their respective charge
            and mass constraints. 

        Returns
        -------
        select : str
            String for each hydrogen bond acceptor/donor/hydrogen atom category.


        .. versionadded:: 2.4.0

        """

        if hasattr(group, "resnames") and hasattr(group, "names"):
            group_list = np.unique([
                '(resname {} and name {})'.format(r, 
                    p) for r, p in zip(group.resnames, group.names)
            ])
        else:
            group_list = np.unique(
                [
                    'type {}'.format(tp) for tp in group.types
                ]
            )

        return " or ".join(group_list)

    def _get_dh_pairs(self):
        """Finds donor-hydrogen pairs.

        Returns
        -------
        donors, hydrogens: AtomGroup, AtomGroup
            AtomGroups corresponding to all donors and all hydrogens. 
            AtomGroups are ordered such that, if zipped, will
            produce a list of donor-hydrogen pairs.
        """

        # # If donors_sel is not provided, use topology to find d-h pairs
        # if self.donors_sel is None:

        #     # We're using u._topology.bonds rather than u.bonds as it is a million times faster to access.
        #     # This is because u.bonds also calculates properties of each bond (e.g bond length).
        #     # See https://github.com/MDAnalysis/mdanalysis/issues/2396#issuecomment-596251787
        #     if not (hasattr(self.u._topology, 'bonds') and len(self.u._topology.bonds.values) != 0):
        #         raise NoDataError('Cannot assign donor-hydrogen pairs via topology as no bond information is present. '
        #                           'Please either: load a topology file with bond information; use the guess_bonds() '
        #                           'topology guesser; or set HydrogenBondAnalysis.donors_sel so that a distance cutoff '
        #                           'can be used.')

        #     hydrogens = self.u.select_atoms(self.hydrogens_sel)
        #     donors = sum(h.bonded_atoms[0] for h in hydrogens) if hydrogens \
        #         else AtomGroup([], self.u)

        # # Otherwise, use d_h_cutoff as a cutoff distance
        # else:

        #     hydrogens = self.guess_hydrogens()
        #     donors = self.guess_donors()
        #     donors_indices, hydrogen_indices = capped_distance(
        #         donors.positions,
        #         hydrogens.positions,
        #         max_cutoff=self.d_h_cutoff,
        #         box=self.u.dimensions,
        #         return_distances=False
        #     ).T

        #     donors = donors[donors_indices]
        #     hydrogens = hydrogens[hydrogen_indices]

        hydrogens = self.hydrogens_ag #NOTE: do I need this?
        donors = self.donors_ag #NOTE: do I need this?
        donors_indices, hydrogen_indices = capped_distance(
            donors.positions,
            hydrogens.positions,
            max_cutoff=self.d_h_cutoff,
            box=self.u.dimensions,
            return_distances=False
        ).T

        donors = donors[donors_indices]
        hydrogens = hydrogens[hydrogen_indices]

        return donors, hydrogens

    def _filter_atoms(self, donors, acceptors):
        """Create a mask to filter donor, hydrogen and acceptor atoms.

        This can be used to consider only hydrogen bonds between two or more
        specified groups.

        Groups are specified with the `between` keyword when creating the
        HydrogenBondAnalysis object.

           Returns
           -------
           mask: np.ndarray


        .. versionchanged:: 2.5.0
           Change return value to a mask instead of separate AtomGroups.
``
        """

        mask = np.full(donors.n_atoms, fill_value=False)
        for group1, group2 in self.between_ags:

            # Find donors in G1 and acceptors in G2
            mask[
                    np.logical_and(
                        np.isin(donors.indices, group1.indices),
                        np.isin(acceptors.indices, group2.indices)
                    )
            ] = True

            # Find acceptors in G1 and donors in G2
            mask[
                np.logical_and(
                    np.isin(acceptors.indices, group1.indices),
                    np.isin(donors.indices, group2.indices)
                )
            ] = True

        return mask


    def _prepare(self):
        self.results.hbonds = [[], [], [], [], [], []]

        # Set unpaired acceptor and hydrogen AtomGroups
        self.acceptors_ag = self.guess_acceptors() #do i need this?
        self.hydrogens_ag = self.guess_hydrogens()
        self.donors_ag = self.guess_donors()

        # Set atom selections if they have not been provided
        if self.acceptors_sel is None:
            self.acceptors_sel = self._group_categories(self.acceptors_ag) #NOTE: do I need this?
        if self.hydrogens_sel is None:
            self.hydrogens_sel = self._group_categories(self.hydrogens_ag)
        #NOTE: add in self.donors_sel? Is this not done for a reason?
        if self.donors_sel is None:
            self.donors_sel = self._group_categories(self.donors_ag)

        # Select atom groups
        self._acceptors = self.guess_acceptors()
        self._donors, self._hydrogens = self._get_dh_pairs()

    def _single_frame(self):

        box = self._ts.dimensions

        # Update donor-hydrogen pairs if necessary
        if self.update_selections:
            self._donors, self._hydrogens = self._get_dh_pairs()

        # find D and A within cutoff distance of one another
        # min_cutoff = 1.0 as an atom cannot form a hydrogen bond with itself
        d_a_indices, d_a_distances = capped_distance(
            self._donors.positions,
            self._acceptors.positions,
            max_cutoff=self.d_a_cutoff,
            min_cutoff=1.0,
            box=box,
            return_distances=True,
        )

        if np.size(d_a_indices) == 0:
            warnings.warn(
                "No hydrogen bonds were found given d-a cutoff of "
                f"{self.d_a_cutoff} between Donor, {self.donors_sel}, and "
                f"Acceptor, {self.acceptors_sel}."
            )

        # Remove D-A pairs more than d_a_cutoff away from one another
        tmp_donors = self._donors[d_a_indices.T[0]]
        tmp_hydrogens = self._hydrogens[d_a_indices.T[0]]
        tmp_acceptors = self._acceptors[d_a_indices.T[1]]

        # Remove donor-acceptor pairs between pairs of AtomGroups we are not
        # interested in
        if self.between_ags is not None:
            between_mask = self._filter_atoms(tmp_donors, tmp_acceptors)
            tmp_donors = tmp_donors[between_mask]
            tmp_hydrogens = tmp_hydrogens[between_mask]
            tmp_acceptors = tmp_acceptors[between_mask]
            d_a_distances = d_a_distances[between_mask]

        # Find D-H-A angles greater than d_h_a_angle_cutoff
        d_h_a_angles = np.rad2deg(
            calc_angles(
                tmp_donors.positions,
                tmp_hydrogens.positions,
                tmp_acceptors.positions,
                box=box
            )
        )
        hbond_indices = np.where(d_h_a_angles > self.d_h_a_angle)[0]

        if np.size(hbond_indices) == 0:
            warnings.warn(
                "No hydrogen bonds were found given angle of "
                f"{self.d_h_a_angle} between Donor, {self.donors_sel}, and "
                f"Acceptor, {self.acceptors_sel}."
            )

        # Retrieve atoms, distances and angles of hydrogen bonds
        hbond_donors = tmp_donors[hbond_indices]
        hbond_hydrogens = tmp_hydrogens[hbond_indices]
        hbond_acceptors = tmp_acceptors[hbond_indices]
        hbond_distances = d_a_distances[hbond_indices]
        hbond_angles = d_h_a_angles[hbond_indices]

        # Store data on hydrogen bonds found at this frame
        self.results.hbonds[0].extend(np.full_like(hbond_donors,
                                      self._ts.frame))
        self.results.hbonds[1].extend(hbond_donors.indices)
        self.results.hbonds[2].extend(hbond_hydrogens.indices)
        self.results.hbonds[3].extend(hbond_acceptors.indices)
        self.results.hbonds[4].extend(hbond_distances)
        self.results.hbonds[5].extend(hbond_angles)

    def _conclude(self):

        self.results.hbonds = np.asarray(self.results.hbonds).T

    @property
    def hbonds(self):
        wmsg = ("The `hbonds` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.hbonds` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.hbonds

    def lifetime(self, tau_max=20, window_step=1, intermittency=0):
        """Computes and returns the time-autocorrelation
        (HydrogenBondLifetimes) of hydrogen bonds.

        Before calling this method, the hydrogen bonds must first be computed
        with the `run()` function. The same `start`, `stop` and `step`
        parameters used in finding hydrogen bonds will be used here for
        calculating hydrogen bond lifetimes. That is, the same frames will be
        used in the analysis.

        Unique hydrogen bonds are identified using hydrogen-acceptor pairs.
        This means an acceptor switching to a different hydrogen atom - with
        the same donor - from one frame to the next is considered a different
        hydrogen bond.

        Please see :func:`MDAnalysis.lib.correlations.autocorrelation` and
        :func:`MDAnalysis.lib.correlations.intermittency` functions for more
        details.


        Parameters
        ----------
        window_step : int, optional
            The number of frames between each t(0).
        tau_max : int, optional
            Hydrogen bond lifetime is calculated for frames in the range
            1 <= `tau` <= `tau_max`
        intermittency : int, optional
            The maximum number of consecutive frames for which a bond can
            disappear but be counted as present if it returns at the next
            frame. An intermittency of `0` is equivalent to a continuous
            autocorrelation, which does not allow for hydrogen bond
            disappearance. For example, for `intermittency=2`, any given
            hydrogen bond may disappear for up to two consecutive frames yet
            be treated as being present at all frames. The default is
            continuous (intermittency=0).

        Returns
        -------
        tau_timeseries : np.array
            tau from 1 to `tau_max`
        timeseries : np.array
            autcorrelation value for each value of `tau`
        """

        if self.results.hbonds is None:
            logger.error(
                "Autocorrelation analysis of hydrogen bonds cannot be done"
                "before the hydrogen bonds are found"
            )
            logger.error(
                "Autocorrelation: Please use the .run() before calling this"
                "function"
            )
            raise NoDataError(".hbonds attribute is None: use .run() first")

        if self.step != 1:
            logger.warning(
                "Autocorrelation: Hydrogen bonds were computed with step > 1."
            )
            logger.warning(
                "Autocorrelation: We recommend recomputing hydrogen bonds with"
                " step = 1."
            )
            logger.warning(
                "Autocorrelation: if you would like to allow bonds to break"
                " and reform, please use 'intermittency'"
            )

        # Extract the hydrogen bonds IDs only in the format
        # [set(superset(x1,x2), superset(x3,x4)), ..]
        found_hydrogen_bonds = [set() for _ in self.frames]
        for frame_index, frame in enumerate(self.frames):
            for hbond in self.results.hbonds[self.results.hbonds[:, 0] == frame]:
                found_hydrogen_bonds[frame_index].add(frozenset(hbond[2:4]))

        intermittent_hbonds = correct_intermittency(
            found_hydrogen_bonds,
            intermittency=intermittency
        )
        tau_timeseries, timeseries, timeseries_data = autocorrelation(
            intermittent_hbonds,
            tau_max,
            window_step=window_step
        )

        return np.vstack([tau_timeseries, timeseries])

    def count_by_time(self):
        """Counts the number of hydrogen bonds per timestep.

        Returns
        -------
        counts : numpy.ndarray
             Contains the total number of hydrogen bonds found at each timestep.
             Can be used along with :attr:`HydrogenBondAnalysis.times` to plot
             the number of hydrogen bonds over time.
        """

        indices, tmp_counts = np.unique(self.results.hbonds[:, 0], axis=0,
                                        return_counts=True)

        indices -= self.start
        indices /= self.step

        counts = np.zeros_like(self.frames)
        counts[indices.astype(np.intp)] = tmp_counts

        return counts

    def count_by_type(self):
        """Counts the total number of each unique type of hydrogen bond.

        Returns
        -------
        counts : numpy.ndarray
             Each row of the array contains the donor resname, donor atom type, 
             acceptor resname, acceptor atom type and the total number of times 
             the hydrogen bond was found.

        Note
        ----
        Unique hydrogen bonds are determined through a consideration of the 
        resname and atom type of the donor and acceptor atoms in a hydrogen bond.
        """

        d = self.u.atoms[self.results.hbonds[:, 1].astype(np.intp)]
        a = self.u.atoms[self.results.hbonds[:, 3].astype(np.intp)]

        if hasattr(d, "resnames"):
            d_res = d.resnames
            a_res = a.resnames
        else:
            d_res = len(d.types) * ["None"]
            a_res = len(a.types) * ["None"]

        tmp_hbonds = np.array([d_res, d.types, a_res, a.types], dtype=str).T
        hbond_type, type_counts = np.unique(
            tmp_hbonds, axis=0, return_counts=True)
        hbond_type_list = []
        for hb_type, hb_count in zip(hbond_type, type_counts):
            hbond_type_list.append([":".join(hb_type[:2]),
                                    ":".join(hb_type[2:4]), hb_count])

        return np.array(hbond_type_list)

    def count_by_ids(self):
        """Counts the total number hydrogen bonds formed by unique combinations of donor, hydrogen and acceptor atoms.

        Returns
        -------
        counts : numpy.ndarray
             Each row of the array contains the donor atom id, hydrogen atom id, acceptor atom id and the total number
             of times the hydrogen bond was observed. The array is sorted by frequency of occurrence.

        Note
        ----
        Unique hydrogen bonds are determined through a consideration of the hydrogen atom id and acceptor atom id
        in a hydrogen bond.
        """

        d = self.u.atoms[self.results.hbonds[:, 1].astype(np.intp)]
        h = self.u.atoms[self.results.hbonds[:, 2].astype(np.intp)]
        a = self.u.atoms[self.results.hbonds[:, 3].astype(np.intp)]

        tmp_hbonds = np.array([d.ids, h.ids, a.ids]).T
        hbond_ids, ids_counts = np.unique(tmp_hbonds, axis=0,
                                          return_counts=True)

        # Find unique hbonds and sort rows so that most frequent observed bonds are at the top of the array
        unique_hbonds = np.concatenate((hbond_ids, ids_counts[:, None]),
                                       axis=1)
        unique_hbonds = unique_hbonds[unique_hbonds[:, 3].argsort()[::-1]]

        return unique_hbonds
