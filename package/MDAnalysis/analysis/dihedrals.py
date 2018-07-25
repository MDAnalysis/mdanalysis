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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
Dihedral and Ramachandran analysis --- :mod:`MDAnalysis.analysis.dihedrals`
===========================================================================

:Author: Henry Mull
:Year: 2018
:Copyright: GNU Public License v2

.. versionadded:: 0.18.1

This module calculates the dihedral angles phi and psi in degrees for a given
list of residues (or atoms corresponding to the residues). This can be done for
selected frames or whole trajectories.

A list of time steps that contain phi and psi angles for each residue is
generated, and a basic Ramachandran plot can be generated using the method
:meth:`Ramachandran.plot()`. This plot is best used as a reference, but it also
allows for user customization.


See Also
--------
:func:`MDAnalysis.lib.distances.calc_dihedrals()`
   function to calculate dihedral angles from atom positions


Example application
-------------------
This example will show how to calculate the phi and psi angles of adenylate
kinase (AdK) and generate a basic Ramachandran plot. The trajectory is included
within the test data files::

   import MDAnalysis as mda
   from MDAnalysisTests.datafiles import GRO, XTC
   u = mda.Universe(GRO, XTC)
   r = u.select_atoms("protein")    # selection of residues

   from MDAnalysis.analysis.dihedrals import Ramachandran
   R = Ramachandran(r).run()

   import matplotlib.pyplot as plt
   fig, ax = plt.subplots(figsize=plt.figaspect(1))
   ax.set_title("Ramachandran Plot (AdK)")
   R.plot(ax=ax, color='k', marker='s')

Alternatively, if you wanted to plot the data yourself, the angles themselves
can be accessed using :attr:`Ramachandran.angles`::

   fig, ax = plt.subplots(figsize=plt.figaspect(1))
   ax.axis([-180,180,-180,180])
   ax.axhline(0, color='k', lw=1)
   ax.axvline(0, color='k', lw=1)
   ax.set(xticks=range(-180,181,60), yticks=range(-180,181,60),
          xlabel=r"$\phi$ (deg)", ylabel=r"$\psi$ (deg)")
   for ts in R.angles:
       ax.scatter(ts[:,0], ts[:,1], color='k', marker='s')

Analysis Class
--------------

.. autoclass:: Ramachandran
   :members:
   :inherited-members:

   .. attribute:: angles

       Contains the time steps of the phi and psi angles for each residue as
       an n_frames×n_residues×2 :class:`numpy.ndarray` with content
       ``[[[phi, psi], [residue 2], ...], [time step 2], ...]``.

"""
from __future__ import absolute_import
from six.moves import zip, range

import numpy as np
import matplotlib.pyplot as plt

import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import calc_dihedrals


class Ramachandran(AnalysisBase):
    """Calculate phi and psi dihedral angles of selected residues.

    Phi and psi angles will be calculated for each residue corresponding to
    `atomgroup` for each time step in the trajectory. A :class:`ResidueGroup`
    is generated from `atomgroup` which is compared to the protein to determine
    if it is a legitimate selection.

    Note
    ----
    If the residue selection is beyond the scope of the protein, then an error
    will be raised. If the residue selection includes the first or last residue,
    then a warning will be raised and they will be removed from the list of
    residues, but the analysis will still run.

    """
    def __init__(self, atomgroup, **kwargs):
        """Parameters
        ----------
        atomgroup : AtomGroup or ResidueGroup
            atoms for residues for which phi and psi are calculated

        Raises
        ------
        ValueError
             If the selection of residues is not contained within the protein

        """
        super(Ramachandran, self).__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup
        residues = self.atomgroup.residues
        protein = self.atomgroup.universe.select_atoms("protein").residues

        if not residues.issubset(protein):
            raise ValueError("Found atoms outside of protein. Only atoms "
                             "inside of a 'protein' selection can be used to "
                             "calculate dihedrals.")
        elif not residues.isdisjoint(protein[[0, -1]]):
            warnings.warn("Cannot determine phi and psi angles for the first "
                          "or last residues")
            residues = residues.difference(protein[[0, -1]])

        phi_sel = [res.phi_selection() for res in residues]
        psi_sel = [res.psi_selection() for res in residues]
        self.ag1 = mda.AtomGroup([atoms[0] for atoms in phi_sel])
        self.ag2 = mda.AtomGroup([atoms[1] for atoms in phi_sel])
        self.ag3 = mda.AtomGroup([atoms[2] for atoms in phi_sel])
        self.ag4 = mda.AtomGroup([atoms[3] for atoms in phi_sel])
        self.ag5 = mda.AtomGroup([atoms[3] for atoms in psi_sel])

    def _prepare(self):
        self.angles = []

    def _single_frame(self):
        phi_angles = calc_dihedrals(self.ag1.positions, self.ag2.positions,
                                    self.ag3.positions, self.ag4.positions,
                                    box=self.ag1.dimensions)
        psi_angles = calc_dihedrals(self.ag2.positions, self.ag3.positions,
                                    self.ag4.positions, self.ag5.positions,
                                    box=self.ag1.dimensions)
        phi_psi = [(phi, psi) for phi, psi in zip(phi_angles, psi_angles)]
        self.angles.append(phi_psi)

    def _conclude(self):
        self.angles = np.rad2deg(np.array(self.angles))


    def plot(self, ax=None, **kwargs):
        """Plots data into standard ramachandran plot. Each time step in
        :attr:`Ramachandran.angles` is plotted onto the same graph.

        Parameters
        ----------
        ax : :class:`matplotlib.axes.Axes`
              If no `ax` is supplied or set to ``None`` then the plot will
              be added to the current active axes.

        Returns
        -------
        ax : :class:`matplotlib.axes.Axes`
             Axes with the plot, either `ax` or the current axes.

        """
        if ax is None:
            ax = plt.gca()
        ax.axis([-180,180,-180,180])
        ax.axhline(0, color='k', lw=1)
        ax.axvline(0, color='k', lw=1)
        ax.set(xticks=range(-180,181,60), yticks=range(-180,181,60),
               xlabel=r"$\phi$ (deg)", ylabel=r"$\psi$ (deg)")
        a = self.angles.reshape(np.prod(self.angles.shape[:2]), 2)
        ax.scatter(a[:,0], a[:,1], **kwargs)
        return ax
