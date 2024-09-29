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
r"""Dihedral angles analysis --- :mod:`MDAnalysis.analysis.dihedrals`
=================================================================

:Author: Henry Mull
:Year: 2018
:Copyright: GNU Public License v2

.. versionadded:: 0.19.0

This module contains classes for calculating dihedral angles for a given set of
atoms or residues. This can be done for selected frames or whole trajectories.

A list of time steps that contain angles of interest is generated and can be
easily plotted if desired. For the :class:`~MDAnalysis.analysis.dihedrals.Ramachandran`
and :class:`~MDAnalysis.analysis.dihedrals.Janin` classes, basic plots can be
generated using the method :meth:`Ramachandran.plot` or :meth:`Janin.plot`.
These plots are best used as references, but they also allow for user customization.


See Also
--------
:func:`MDAnalysis.lib.distances.calc_dihedrals()`
   function to calculate dihedral angles from atom positions


Example applications
--------------------

General dihedral analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`~MDAnalysis.analysis.dihedrals.Dihedral` class is useful for calculating
angles for many dihedrals of interest. For example, we can find the phi angles
for residues 5-10 of adenylate kinase (AdK). The trajectory is included within
the test data files::

   import MDAnalysis as mda
   from MDAnalysisTests.datafiles import GRO, XTC
   u = mda.Universe(GRO, XTC)

   # selection of atomgroups
   ags = [res.phi_selection() for res in u.residues[4:9]]

   from MDAnalysis.analysis.dihedrals import Dihedral
   R = Dihedral(ags).run()

The angles can then be accessed with :attr:`Dihedral.results.angles`.

Ramachandran analysis
~~~~~~~~~~~~~~~~~~~~~

The :class:`~MDAnalysis.analysis.dihedrals.Ramachandran` class allows for the
quick calculation of classical Ramachandran plots :footcite:p:`Ramachandran1963` in
the backbone :math:`phi` and :math:`psi` angles. Unlike the
:class:`~MDanalysis.analysis.dihedrals.Dihedral` class which takes a list of
`atomgroups`, this class only needs a list of residues or atoms from those
residues. The previous example can repeated with::

   u = mda.Universe(GRO, XTC)
   r = u.select_atoms("resid 5-10")

   R = Ramachandran(r).run()

Then it can be plotted using the built-in plotting method :meth:`Ramachandran.plot()`::

   import matplotlib.pyplot as plt
   fig, ax = plt.subplots(figsize=plt.figaspect(1))
   R.plot(ax=ax, color='k', marker='o', ref=True)
   fig.tight_layout()

as shown in the example :ref:`Ramachandran plot figure <figure-ramachandran>`.

.. _figure-ramachandran:

.. figure:: /images/rama_demo_plot.png
   :scale: 50 %
   :alt: Ramachandran plot

   Ramachandran plot for residues 5 to 10 of AdK, sampled from the AdK test
   trajectory (XTC). The contours in the background are the "allowed region"
   and the "marginally allowed" regions.

To plot the data yourself, the angles can be accessed using
:attr:`Ramachandran.results.angles`.

.. Note::

   The Ramachandran analysis is prone to errors if the topology contains
   duplicate or missing atoms (e.g. atoms with `altloc` or incomplete
   residues). If the topology has as an `altloc` attribute, you must specify
   only one `altloc` for the atoms with more than one (``"protein and not
   altloc B"``).


Janin analysis
~~~~~~~~~~~~~~

Janin plots :footcite:p:`Janin1978` for side chain conformations (:math:`\chi_1`
and :math:`chi_2` angles) can be created with the
:class:`~MDAnalysis.analysis.dihedrals.Janin` class. It works in the same way,
only needing a list of residues; see the :ref:`Janin plot figure
<figure-janin>` as an example.

The data for the angles can be accessed in the attribute
:attr:`Janin.results.angles`.

.. _figure-janin:

.. figure:: /images/janin_demo_plot.png
   :scale: 50 %
   :alt: Janin plot

   Janin plot for all residues of AdK, sampled from the AdK test trajectory
   (XTC). The contours in the background are the "allowed region" and the
   "marginally allowed" regions for all possible residues.

.. Note::

   The Janin analysis is prone to errors if the topology contains duplicate or
   missing atoms (e.g. atoms with `altloc` or incomplete residues). If the
   topology has as an `altloc` attribute, you must specify only one `altloc`
   for the atoms with more than one (``"protein and not altloc B"``).

   Furthermore, many residues do not have a :math:`\chi_2` dihedral and if the
   selections of residues is not carefully filtered to only include those
   residues with *both* sidechain dihedrals then a :exc:`ValueError` with the
   message *Too many or too few atoms selected* is raised.


Reference plots
~~~~~~~~~~~~~~~

Reference plots can be added to the axes for both the Ramachandran and Janin
classes using the kwarg ``ref=True`` for the :meth:`Ramachandran.plot`
and :meth:`Janin.plot` methods. The Ramachandran reference data
(:data:`~MDAnalysis.analysis.data.filenames.Rama_ref`) and Janin reference data
(:data:`~MDAnalysis.analysis.data.filenames.Janin_ref`) were made using data
obtained from a large selection of 500 PDB files, and were analyzed using these
classes :footcite:p:`Mull2018`. The allowed and marginally allowed regions of the
Ramachandran reference plot have cutoffs set to include 90% and 99% of the data
points, and the Janin reference plot has cutoffs for 90% and 98% of the data
points. The list of PDB files used for the reference plots was taken from
:footcite:p:`Lovell2003` and information about general Janin regions was taken from
:footcite:p:`Janin1978`.



Analysis Classes
----------------

.. autoclass:: Dihedral
   :members:
   :inherited-members:

   .. attribute:: results.angles

       Contains the time steps of the angles for each atomgroup in the list as
       an ``n_frames×len(atomgroups)`` :class:`numpy.ndarray` with content
       ``[[angle 1, angle 2, ...], [time step 2], ...]``.

       .. versionadded:: 2.0.0

   .. attribute:: angles

       Alias to the :attr:`results.angles` attribute.

       .. deprecated:: 2.0.0
          Will be removed in MDAnalysis 3.0.0. Please use
          :attr:`results.angles` instead.


.. autoclass:: Ramachandran
   :members:
   :inherited-members:

   .. attribute:: results.angles

       Contains the time steps of the :math:`\phi` and :math:`\psi` angles for
       each residue as an ``n_frames×n_residues×2`` :class:`numpy.ndarray` with
       content ``[[[phi, psi], [residue 2], ...], [time step 2], ...]``.

       .. versionadded:: 2.0.0

   .. attribute:: angles

       Alias to the :attr:`results.angles` attribute.

       .. deprecated:: 2.0.0
          Will be removed in MDAnalysis 3.0.0. Please use
          :attr:`results.angles` instead.


.. autoclass:: Janin
   :members:
   :inherited-members:

   .. attribute:: results.angles

       Contains the time steps of the :math:`\chi_1` and :math:`\chi_2` angles
       for each residue as an ``n_frames×n_residues×2`` :class:`numpy.ndarray`
       with content ``[[[chi1, chi2], [residue 2], ...], [time step 2], ...]``.

       .. versionadded:: 2.0.0

   .. attribute:: angles

       Alias to the :attr:`results.angles` attribute.

       .. deprecated:: 2.0.0
          Will be removed in MDAnalysis 3.0.0. Please use
          :attr:`results.angles` instead.


References
----------

.. footbibliography::

"""
import numpy as np
import matplotlib.pyplot as plt

import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase, ResultsGroup
from MDAnalysis.lib.distances import calc_dihedrals
from MDAnalysis.analysis.data.filenames import Rama_ref, Janin_ref


class Dihedral(AnalysisBase):
    """Calculate dihedral angles for specified atomgroups.

    Dihedral angles will be calculated for each atomgroup that is given for
    each step in the trajectory. Each :class:`~MDAnalysis.core.groups.AtomGroup`
    must contain 4 atoms.

    Note
    ----
    This class takes a list as an input and is most useful for a large
    selection of atomgroups. If there is only one atomgroup of interest, then
    it must be given as a list of one atomgroup.


    .. versionchanged:: 2.0.0
       :attr:`angles` results are now stored in a
       :class:`MDAnalysis.analysis.base.Results` instance.
    .. versionchanged:: 2.8.0
       introduced :meth:`get_supported_backends` allowing for parallel
       execution on ``multiprocessing`` and ``dask`` backends.
    """
    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(cls):
        return ('serial', 'multiprocessing', 'dask',)


    def __init__(self, atomgroups, **kwargs):
        """Parameters
        ----------
        atomgroups : list[AtomGroup]
            a list of :class:`~MDAnalysis.core.groups.AtomGroup` for which
            the dihedral angles are calculated

        Raises
        ------
        ValueError
            If any atomgroups do not contain 4 atoms

        """
        super(Dihedral, self).__init__(
            atomgroups[0].universe.trajectory, **kwargs)
        self.atomgroups = atomgroups

        if any([len(ag) != 4 for ag in atomgroups]):
            raise ValueError("All AtomGroups must contain 4 atoms")

        self.ag1 = mda.AtomGroup([ag[0] for ag in atomgroups])
        self.ag2 = mda.AtomGroup([ag[1] for ag in atomgroups])
        self.ag3 = mda.AtomGroup([ag[2] for ag in atomgroups])
        self.ag4 = mda.AtomGroup([ag[3] for ag in atomgroups])

    def _prepare(self):
        self.results.angles = []

    def _get_aggregator(self):
        return ResultsGroup(lookup={'angles': ResultsGroup.ndarray_vstack})

    def _single_frame(self):
        angle = calc_dihedrals(self.ag1.positions, self.ag2.positions,
                               self.ag3.positions, self.ag4.positions,
                               box=self.ag1.dimensions)
        self.results.angles.append(angle)

    def _conclude(self):
        self.results.angles = np.rad2deg(np.array(self.results.angles))

    @property
    def angles(self):
        wmsg = ("The `angle` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.angles` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.angles


class Ramachandran(AnalysisBase):
    r"""Calculate :math:`\phi` and :math:`\psi` dihedral angles of selected
    residues.

    :math:`\phi` and :math:`\psi` angles will be calculated for each residue
    corresponding to `atomgroup` for each time step in the trajectory. A
    :class:`~MDAnalysis.ResidueGroup` is generated from `atomgroup` which is
    compared to the protein to determine if it is a legitimate selection.

    Parameters
    ----------
    atomgroup : AtomGroup or ResidueGroup
        atoms for residues for which :math:`\phi` and :math:`\psi` are
        calculated
    c_name : str (optional)
        name for the backbone C atom
    n_name : str (optional)
        name for the backbone N atom
    ca_name : str (optional)
        name for the alpha-carbon atom
    check_protein : bool (optional)
        whether to raise an error if the provided atomgroup is not a
        subset of protein atoms

    Example
    -------
    For standard proteins, the default arguments will suffice to run a
    Ramachandran analysis::

        r = Ramachandran(u.select_atoms('protein')).run()

    For proteins with non-standard residues, or for calculating dihedral
    angles for other linear polymers, you can switch off the protein checking
    and provide your own atom names in place of the typical peptide backbone
    atoms::

        r = Ramachandran(u.atoms, c_name='CX', n_name='NT', ca_name='S',
                         check_protein=False).run()

    The above analysis will calculate angles from a "phi" selection of
    CX'-NT-S-CX and "psi" selections of NT-S-CX-NT'.

    Raises
    ------
    ValueError
        If the selection of residues is not contained within the protein
        and ``check_protein`` is ``True``

    Note
    ----
    If ``check_protein`` is ``True`` and the residue selection is beyond
    the scope of the protein and, then an error will be raised.
    If the residue selection includes the first or last residue,
    then a warning will be raised and they will be removed from the list of
    residues, but the analysis will still run. If a :math:`\phi` or :math:`\psi`
    selection cannot be made, that residue will be removed from the analysis.


    .. versionchanged:: 1.0.0
        added c_name, n_name, ca_name, and check_protein keyword arguments
    .. versionchanged:: 2.0.0
       :attr:`angles` results are now stored in a
       :class:`MDAnalysis.analysis.base.Results` instance.
    .. versionchanged:: 2.8.0
       introduced :meth:`get_supported_backends` allowing for parallel
       execution on ``multiprocessing`` and ``dask`` backends.
    """
    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(cls):
        return ('serial', 'multiprocessing', 'dask',)

    def __init__(self, atomgroup, c_name='C', n_name='N', ca_name='CA',
                 check_protein=True, **kwargs):
        super(Ramachandran, self).__init__(
            atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup
        residues = self.atomgroup.residues

        if check_protein:
            protein = self.atomgroup.universe.select_atoms("protein").residues

            if not residues.issubset(protein):
                raise ValueError("Found atoms outside of protein. Only atoms "
                                "inside of a 'protein' selection can be used to "
                                "calculate dihedrals.")
            elif not residues.isdisjoint(protein[[0, -1]]):
                warnings.warn("Cannot determine phi and psi angles for the first "
                            "or last residues")
                residues = residues.difference(protein[[0, -1]])

        prev = residues._get_prev_residues_by_resid()
        nxt = residues._get_next_residues_by_resid()
        keep = np.array([r is not None for r in prev])
        keep = keep & np.array([r is not None for r in nxt])

        if not np.all(keep):
            warnings.warn("Some residues in selection do not have "
                          "phi or psi selections")
        prev = sum(prev[keep])
        nxt = sum(nxt[keep])
        residues = residues[keep]

        # find n, c, ca
        keep_prev = [sum(r.atoms.names==c_name)==1 for r in prev]
        rnames = [n_name, c_name, ca_name]
        keep_res = [all(sum(r.atoms.names == n) == 1 for n in rnames)
                    for r in residues]
        keep_next = [sum(r.atoms.names == n_name) == 1 for r in nxt]

        # alright we'll keep these
        keep = np.array(keep_prev) & np.array(keep_res) & np.array(keep_next)
        prev = prev[keep]
        res = residues[keep]
        nxt = nxt[keep]

        rnames = res.atoms.names
        self.ag1 = prev.atoms[prev.atoms.names == c_name]
        self.ag2 = res.atoms[rnames == n_name]
        self.ag3 = res.atoms[rnames == ca_name]
        self.ag4 = res.atoms[rnames == c_name]
        self.ag5 = nxt.atoms[nxt.atoms.names == n_name]


    def _prepare(self):
        self.results.angles = []

    def _get_aggregator(self):
        return ResultsGroup(lookup={'angles': ResultsGroup.ndarray_vstack})

    def _single_frame(self):
        phi_angles = calc_dihedrals(self.ag1.positions, self.ag2.positions,
                                    self.ag3.positions, self.ag4.positions,
                                    box=self.ag1.dimensions)
        psi_angles = calc_dihedrals(self.ag2.positions, self.ag3.positions,
                                    self.ag4.positions, self.ag5.positions,
                                    box=self.ag1.dimensions)
        phi_psi = [(phi, psi) for phi, psi in zip(phi_angles, psi_angles)]
        self.results.angles.append(phi_psi)

    def _conclude(self):
        self.results.angles = np.rad2deg(np.array(self.results.angles))

    def plot(self, ax=None, ref=False, **kwargs):
        """Plots data into standard Ramachandran plot.

        Each time step in :attr:`Ramachandran.results.angles` is plotted onto
        the same graph.

        Parameters
        ----------
        ax : :class:`matplotlib.axes.Axes`
              If no `ax` is supplied or set to ``None`` then the plot will
              be added to the current active axes.

        ref : bool, optional
              Adds a general Ramachandran plot which shows allowed and
              marginally allowed regions

        kwargs : optional
              All other kwargs are passed to :func:`matplotlib.pyplot.scatter`.

        Returns
        -------
        ax : :class:`matplotlib.axes.Axes`
             Axes with the plot, either `ax` or the current axes.

        """
        if ax is None:
            ax = plt.gca()
        ax.axis([-180, 180, -180, 180])
        ax.axhline(0, color='k', lw=1)
        ax.axvline(0, color='k', lw=1)
        ax.set(xticks=range(-180, 181, 60), yticks=range(-180, 181, 60),
               xlabel=r"$\phi$", ylabel=r"$\psi$")
        degree_formatter = plt.matplotlib.ticker.StrMethodFormatter(
            r"{x:g}$\degree$")
        ax.xaxis.set_major_formatter(degree_formatter)
        ax.yaxis.set_major_formatter(degree_formatter)

        if ref:
            X, Y = np.meshgrid(np.arange(-180, 180, 4),
                               np.arange(-180, 180, 4))
            levels = [1, 17, 15000]
            colors = ['#A1D4FF', '#35A1FF']
            ax.contourf(X, Y, np.load(Rama_ref), levels=levels, colors=colors)
        a = self.results.angles.reshape(
                np.prod(self.results.angles.shape[:2]), 2)
        ax.scatter(a[:, 0], a[:, 1], **kwargs)
        return ax

    @property
    def angles(self):
        wmsg = ("The `angle` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.angles` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.angles


class Janin(Ramachandran):
    r"""Calculate :math:`\chi_1` and :math:`\chi_2` dihedral angles of selected
    residues.

    :math:`\chi_1` and :math:`\chi_2` angles will be calculated for each residue
    corresponding to `atomgroup` for each time step in the trajectory. A
    :class:`~MDAnalysis.ResidueGroup` is generated from `atomgroup` which is
    compared to the protein to determine if it is a legitimate selection.

    Note
    ----
    If the residue selection is beyond the scope of the protein, then an error
    will be raised. If the residue selection includes the residues ALA, CYS*,
    GLY, PRO, SER, THR, or VAL (the default of the `select_remove` keyword
    argument) then a warning will be raised and they will be removed from the
    list of residues, but the analysis will still run. Some topologies have
    altloc attributes which can add duplicate atoms to the selection and must
    be removed.

    """

    def __init__(self, atomgroup,
                 select_remove="resname ALA CYS* GLY PRO SER THR VAL",
                 select_protein="protein",
                 **kwargs):
        r"""Parameters
        ----------
        atomgroup : AtomGroup or ResidueGroup
            atoms for residues for which :math:`\chi_1` and :math:`\chi_2` are
            calculated

        select_remove : str
            selection string to remove residues that do not have :math:`chi_2`
            angles

        select_protein : str
            selection string to subselect protein-only residues from
            `atomgroup` to check that only amino acids are selected; if you
            have non-standard amino acids then adjust this selection to include
            them

        Raises
        ------
        ValueError
             if the final selection of residues is not contained within the
             protein (as determined by
             ``atomgroup.select_atoms(select_protein)``)

        ValueError
             if not enough or too many atoms are found for a residue in the
             selection, usually due to missing atoms or alternative locations,
             or due to non-standard residues


        .. versionchanged:: 2.0.0
           `select_remove` and `select_protein` keywords were added.
           :attr:`angles` results are now stored in a
           :class:`MDAnalysis.analysis.base.Results` instance.
        """
        super(Ramachandran, self).__init__(
            atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup
        residues = atomgroup.residues
        protein = atomgroup.select_atoms(select_protein).residues
        remove = residues.atoms.select_atoms(select_remove).residues

        if not residues.issubset(protein):
            raise ValueError("Found atoms outside of protein. Only atoms "
                             "inside of a protein "
                             f"(select_protein='{select_protein}') can be "
                             "used to calculate dihedrals.")
        elif len(remove) != 0:
            warnings.warn(f"All residues selected with '{select_remove}' "
                          "have been removed from the selection.")
            residues = residues.difference(remove)

        self.ag1 = residues.atoms.select_atoms("name N")
        self.ag2 = residues.atoms.select_atoms("name CA")
        self.ag3 = residues.atoms.select_atoms("name CB")
        self.ag4 = residues.atoms.select_atoms("name CG CG1")
        self.ag5 = residues.atoms.select_atoms("name CD CD1 OD1 ND1 SD")

        # if there is an altloc attribute, too many atoms will be selected which
        # must be removed before using the class, or the file is missing atoms
        # for some residues which must also be removed
        if any(len(self.ag1) != len(ag) for ag in [self.ag2, self.ag3,
                                                   self.ag4, self.ag5]):
            raise ValueError("Too many or too few atoms selected. Check for "
                             "missing or duplicate atoms in topology.")

    def _conclude(self):
        self.results.angles = (np.rad2deg(np.array(
            self.results.angles)) + 360) % 360

    def plot(self, ax=None, ref=False, **kwargs):
        """Plots data into standard Janin plot.

        Each time step in :attr:`Janin.results.angles` is plotted onto the
        same graph.

        Parameters
        ----------
        ax : :class:`matplotlib.axes.Axes`
              If no `ax` is supplied or set to ``None`` then the plot will
              be added to the current active axes.

        ref : bool, optional
              Adds a general Janin plot which shows allowed and marginally
              allowed regions

        kwargs : optional
              All other kwargs are passed to :func:`matplotlib.pyplot.scatter`.

        Returns
        -------
        ax : :class:`matplotlib.axes.Axes`
             Axes with the plot, either `ax` or the current axes.

        """
        if ax is None:
            ax = plt.gca()
        ax.axis([0, 360, 0, 360])
        ax.axhline(180, color='k', lw=1)
        ax.axvline(180, color='k', lw=1)
        ax.set(xticks=range(0, 361, 60), yticks=range(0, 361, 60),
               xlabel=r"$\chi_1$", ylabel=r"$\chi_2$")
        degree_formatter = plt.matplotlib.ticker.StrMethodFormatter(
            r"{x:g}$\degree$")
        ax.xaxis.set_major_formatter(degree_formatter)
        ax.yaxis.set_major_formatter(degree_formatter)

        if ref:
            X, Y = np.meshgrid(np.arange(0, 360, 6), np.arange(0, 360, 6))
            levels = [1, 6, 600]
            colors = ['#A1D4FF', '#35A1FF']
            ax.contourf(X, Y, np.load(Janin_ref), levels=levels, colors=colors)
        a = self.results.angles.reshape(np.prod(
            self.results.angles.shape[:2]), 2)
        ax.scatter(a[:, 0], a[:, 1], **kwargs)
        return ax
