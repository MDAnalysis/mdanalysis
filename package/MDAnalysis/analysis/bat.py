# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v21 or any higher version
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
r"""Bond-Angle-Torsion coordinates analysis --- :mod:`MDAnalysis.analysis.bat`
===========================================================================

:Author: David Minh
:Year: 2020
:Copyright: GNU Public License v21

.. versionadded:: N/A

This module contains classes for interconverting between Cartesian and an
internal coordinate system, Bond-Angle-Torsion (BAT) coordinates [Chang2003]_,
for a given set of atoms or residues. This coordinate system is designed
to be complete, non-redundant, and minimize correlations between degrees
of freedom. Complete and non-redundant means that for N atoms there will
be 3N Cartesian coordinates and 3N BAT coordinates. Correlations are
minimized by using improper torsions [Hikiri2016]_.

More specifically, bond refers to the bond length, or distance between
a pair of bonded atoms. Angle refers to the bond angle, the angle between
a pair of bonds to a central atom. Torsion refers to the torsion angle.
For a set of four atoms a, b, c, and d, a torsion requires bonds between
a and b, b and c, and c and d. The torsion is the angle between a plane
containing atoms a, b, and c and another plane containing b, c, and d.
For a set of torsions that share atoms b and c, one torsion is defined as
the primary torsion. The others are defined as improper torsions, differences
between the raw torsion angle and the primary torsion. This definition reduces
the correlation between the torsion angles.

Each molecule also has six external coordinates that define its translation and
rotation in space. The three Cartesian coordinates of the first atom are the
molecule's translational degrees of freedom. Rotational degreres of freedom are
specified by the axis-angle convention. The rotation axis is a normalized vector
pointing from the first to second atom. It is described by the polar angle, phi,
and azimuthal angle, theta. omega is a third angle that describes the rotation
of the third atom about the axis.

This module was adapted from AlGDock [Minh2020]_.


See Also
--------
:class:`~MDAnalysis.analysis.dihedrals.Dihedral`
   class to calculate dihedral angles for a given set of atoms or residues
:func:`MDAnalysis.lib.distances.calc_dihedrals()`
   function to calculate dihedral angles from atom positions


Example applications
--------------------

The :class:`~MDAnalysis.analysis.bat.BAT` class defines bond-angle-torsion
coordinates based on the topology of an atom group and interconverts between
Cartesian and BAT coordinate systems. For example, we can determine internal
coordinates for residues 5-10 of adenylate kinase (AdK). The trajectory is
included within the test data files::

   import MDAnalysis as mda
   from MDAnalysisTests.datafiles import PSF, DCD
   import numpy as np

   u = mda.Universe(PSF, DCD)

   # selection of atomgroups
   selected_residues = u.select_atoms("resid 5-10")

   from MDAnalysis.analysis.bat import BAT
   R = BAT(selected_residues)

   # Calculate BAT coordinates for a trajectory
   R.run()

   # Reconstruct Cartesian coordinates from BAT coordinates
   bat = R.bat[0]
   XYZ = R.Cartesian(bat)

   # The difference between the original and reconstructed coordinates
   # should be zero.
   print(np.sum(np.abs(XYZ - selected_residues.positions)>1E-6))

After R.run(), the coordinates can be accessed with :attr:`R.bat`.


References
----------

.. [Chang2003] Chang, Chia-En, Michael J Potter, and Michael K Gilson (2003).
   "Calculation of Molecular Configuration Integrals". *Journal of Physical
   Chemistry B* 107(4): 1048–55. doi:`10.1021/jp027149c
   <https://doi.org/10.1021/jp027149c>`_

.. [Hikiri2016] Hikiri, Simon, Takashi Yoshidome, and Mitsunori Ikeguchi (2016).
   "Computational Methods for Configurational Entropy Using Internal and
   Cartesian Coordinates." *Journal of Chemical Theory and Computation*
   12(12): 5990–6000. doi:`10.1021/acs.jctc.6b00563
   <https://doi.org/10.1021/acs.jctc.6b00563>`_

.. [Minh2020] Minh, David D L (2020). "Alchemical Grid Dock (AlGDock): Binding
   Free Energy Calculations between Flexible Ligands and Rigid Receptors."
   *Journal of Computational Chemistry* 41(7): 715–30.
   doi:`10.1002/jcc.26036 <https://doi.org/10.1002/jcc.26036>`_

"""
from __future__ import absolute_import
# from six.moves import zip, range

import numpy as np

import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

from MDAnalysis.lib.distances import calc_bonds, calc_angles, calc_dihedrals


def _sort_atoms_by_mass(atoms, reverse=False):
    r""" Sorts a list of atoms by name and then by mass

    Parameters
    ----------
    ag_o : list of Atoms
      List to sort
    reverse : bool
      Atoms will be in descending order of mass

    Returns
    -------
    ag_n : list of Atoms
      Sorted list
    """
    return sorted(sorted(atoms, key = lambda atom:atom.name), \
        key = lambda atom:atom.mass, reverse = reverse)


def _find_torsion(selected_atoms, allowed_atoms):
    """ Finds a torsion angle adjacent to the selected atoms

    The torsion angle includes an atom that is allowed_atoms and is not in
    selected_atoms.

    Parameters
    ----------
    selected_atoms : AtomGroup
      Atoms that have already been selected
    allowed_atoms : AtomGroup
      Atoms that are allowed to be part of the torsion angle

    Returns
    -------
    new_torsion : AtomGroup
      an AtomGroup that defines the torsion angle
    """
    for a1 in selected_atoms:
        # Loop over new atoms connected to the selected atom
        for a0 in _sort_atoms_by_mass([a for a in a1.bonded_atoms \
            if (a in allowed_atoms) and (not a in selected_atoms)]):
            # Find the third atom
            for a2 in _sort_atoms_by_mass([a for a in a1.bonded_atoms \
                if (a in allowed_atoms) and (a in selected_atoms) and (a!=a0)]):
                # Find the fourth atom
                for a3 in _sort_atoms_by_mass([a for a in a2.bonded_atoms \
                    if (a in allowed_atoms) and (a in selected_atoms) and (a!=a1)]):
                    return mda.AtomGroup([a0, a1, a2, a3])
    raise Exception('No new torsion angle found!')


class BAT(AnalysisBase):
    """Calculate BAT coordinates for the specified AtomGroup.

    BAT coordinates will be computed for the group of atoms and all frames
    in the trajectory belonging to `ag'.`

    """
    def __init__(self, ag, initial_atom=None, **kwargs):
        r"""Parameters
        ----------
        ag : AtomGroup or Universe
            Group of atoms for which the BAT coordinates are calculated.
            ag must have a bonds attribute. If not available,
            bonds may be guessed using
            :meth:`AtomGroup.guess_bonds <MDAnalysis.core.groups.AtomGroup.guess_bonds>`.
            ag must only include one molecule.
            If a trajectory is associated with the atoms then the computation
            iterates over the trajectory.
        initial_atom : Atom
            The atom whose Cartesian coordinates define the translation
            of the molecule. If not specified, the heaviest terminal atom
            will be selected.

        Raises
        ------
        AttributeError
            If ag does not contain a bonds attribute
        ValueError
            If ag contains more than one molecule

        """
        super(BAT, self).__init__(ag.universe.trajectory, **kwargs)
        self._ag = ag

        # Check that the ag contains bonds
        if not hasattr(ag, 'bonds'):
            raise AttributeError('AtomGroup has no attribute bonds')

        # Determine the root
        # The initial atom must be a terminal atom
        terminal_atoms = _sort_atoms_by_mass(\
          [a for a in self._ag.atoms if len(a.bonds)==1])
        if (initial_atom is None):
            # Select the heaviest root atoms from the heaviest terminal atom
            initial_atom = terminal_atoms[-1]
        elif (not initial_atom in terminal_atoms):
            raise Exception('Initial atom is not a terminal atom')
        # The next atom in the root is bonded to the initial atom
        second_atom = initial_atom.bonded_atoms[0]
        # The last atom in the root is the heaviest atom
        # bonded to the second atom
        third_atom = [a for a in second_atom.bonded_atoms \
          if (a in self._ag) and (a!=initial_atom)][-1]
        self._root = mda.AtomGroup([initial_atom, second_atom, third_atom])

        # Construct a list of torsion angles
        self._torsions = []
        selected_atoms = mda.AtomGroup(self._root)
        while len(selected_atoms) < self._ag.n_atoms:
            try:
                torsion = _find_torsion(selected_atoms, self._ag)
            except:
                raise ValueError('AtomGroup is more than one molecule')

            self._torsions.append(torsion)
            selected_atoms += torsion[0]

        # Get indices of the root and torsion atoms
        # in a Cartesian positions array that matches the AtomGroup
        ag_top_inds = list(self._ag.indices)
        self._root_XYZ_inds = [ag_top_inds.index(a.index) for a in self._root]
        self._torsion_XYZ_inds = [[ag_top_inds.index(a.index) for a in t] \
          for t in self._torsions]

        # The primary torsion is the first torsion on the list
        # with the same central atoms
        prior_atoms = [sorted([a1, a2]) for (a0, a1, a2, a3) in self._torsions]
        self._primary_torsion_indices = [prior_atoms.index(prior_atoms[n]) \
          for n in range(len(prior_atoms))]

        self._ag1 = mda.AtomGroup([ag[0] for ag in self._torsions])
        self._ag2 = mda.AtomGroup([ag[1] for ag in self._torsions])
        self._ag3 = mda.AtomGroup([ag[2] for ag in self._torsions])
        self._ag4 = mda.AtomGroup([ag[3] for ag in self._torsions])

    def _prepare(self):
        self.bat = []

    def _single_frame(self):
        # Calculate coordinates based on the root atoms
        # The rotation axis is a normalized vector pointing from atom 0 to 1
        # It is described in two degrees of freedom
        # by the polar angle and azimuth
        (p0, p1, p2) = self._root.positions
        v01 = p1 - p0
        v21 = p1 - p2
        # Internal coordinates
        r01 = np.sqrt(np.sum(v01 *
                             v01))  # Distance between first two root atoms
        r12 = np.sqrt(np.sum(v21 *
                             v21))  # Distance between second two root atoms
        a012 = np.arccos(max(-1.,min(1.,np.sum(v01*v21)/\
          np.sqrt(np.sum(v01*v01)*np.sum(v21*v21))))) # Angle between root atoms
        # Exernal coordinates
        e = v01 / r01
        phi = np.arctan2(e[1], e[0])  # Polar angle
        theta = np.arccos(e[2])  # Azimuthal angle
        # Rotation to the z axis
        cp = np.cos(phi)
        sp = np.sin(phi)
        ct = np.cos(theta)
        st = np.sin(theta)
        Rz = np.array([[cp * ct, ct * sp, -st], [-sp, cp, 0],
                       [cp * st, sp * st, ct]])
        pos2 = Rz.dot(p2 - p1)
        # Angle about the rotation axis
        omega = np.arctan2(pos2[1], pos2[0])
        root_based = np.concatenate((p0, [phi, theta, omega, r01, r12, a012]))

        # Calculate internal coordinates from the torsion list
        bonds = calc_bonds(self._ag1.positions,
                           self._ag2.positions,
                           box=self._ag1.dimensions)
        angles = calc_angles(self._ag1.positions,
                             self._ag2.positions,
                             self._ag3.positions,
                             box=self._ag1.dimensions)
        torsions = calc_dihedrals(self._ag1.positions,
                                  self._ag2.positions,
                                  self._ag3.positions,
                                  self._ag4.positions,
                                  box=self._ag1.dimensions)
        # When appropriate, calculate improper torsions
        torsions = np.array([\
          torsions[n] - torsions[self._primary_torsion_indices[n]] \
            if self._primary_torsion_indices[n]!=n else torsions[n] \
              for n in range(len(torsions))])
        # Wrap torsions to between -np.pi and np.pi
        torsions = ((torsions + np.pi) % (2 * np.pi)) - np.pi

        self.bat.append(np.concatenate((root_based, bonds, angles, torsions)))

    def Cartesian(self, bat):
        """Conversion from BAT to Cartesian coordinates

        Parameters
        ----------
        bat : np.array
          an array with dimensions (3N,) array with external then internal
          degrees of freedom based on the root atoms, followed by the bond,
          angle, and (proper and improper) torsion coordinates.

        Returns
        -------
        xyz : np.array
          an array with dimensions (N,3) with Cartesian coordinates. The first
          dimension has the same ordering as the AtomGroup used to initialize
          the class.
        """
        # Split the bat vector into more convenient variables
        origin = bat[:3]
        (phi, theta, omega) = bat[3:6]
        (r01, r12, a012) = bat[6:9]
        n_torsions = (self._ag.n_atoms - 3)
        bonds = bat[9:n_torsions + 9]
        angles = bat[n_torsions + 9:2 * n_torsions + 9]
        torsions = bat[2 * n_torsions + 9:]
        # When appropriate, convert improper to proper torsions
        torsions = np.array([\
          torsions[n] + torsions[self._primary_torsion_indices[n]] \
            if self._primary_torsion_indices[n]!=n else torsions[n] \
              for n in range(len(torsions))])
        # Wrap torsions to between -np.pi and np.pi
        torsions = ((torsions + np.pi) % (2 * np.pi)) - np.pi

        # Set initial root atom positions based on internal coordinates
        p0 = np.array([0., 0., 0.])
        p1 = np.array([0., 0., r01])
        p2 = np.array([r12 * np.sin(a012), 0., r01 - r12 * np.cos(a012)])

        # Rotate the third atom by the appropriate value
        co = np.cos(omega)
        so = np.sin(omega)
        Romega = np.array([[co, -so, 0], [so, co, 0], [0, 0, 1]])
        p2 = Romega.dot(p2)
        # Rotate the second two atoms to point in the right direction
        cp = np.cos(phi)
        sp = np.sin(phi)
        ct = np.cos(theta)
        st = np.sin(theta)
        Re = np.array([[cp * ct, -sp, cp * st], [ct * sp, cp, sp * st],
                       [-st, 0, ct]])
        p1 = Re.dot(p1)
        p2 = Re.dot(p2)
        # Translate the first three atoms by the origin
        p0 += origin
        p1 += origin
        p2 += origin

        XYZ = np.zeros((self._ag.n_atoms, 3))
        XYZ[self._root_XYZ_inds[0]] = p0
        XYZ[self._root_XYZ_inds[1]] = p1
        XYZ[self._root_XYZ_inds[2]] = p2

        # Place the remaining atoms
        for ((a0,a1,a2,a3), r01, angle, torsion) \
            in zip(self._torsion_XYZ_inds, bonds, angles, torsions):

            p1 = XYZ[a1]
            p3 = XYZ[a3]
            p2 = XYZ[a2]

            sn_ang = np.sin(angle)
            cs_ang = np.cos(angle)
            sn_tor = np.sin(torsion)
            cs_tor = np.cos(torsion)

            v21 = p1 - p2
            v21 /= np.sqrt(np.sum(v21 * v21))
            v32 = p2 - p3
            v32 /= np.sqrt(np.sum(v32 * v32))

            vp = np.cross(v32, v21)
            cs = np.sum(v21 * v32)
            if abs(cs) > 1:
                print('cos ', cs)

            sn = np.sqrt(max(1.0 - cs * cs, 0.0000000001))
            vp = vp / sn
            vu = np.cross(vp, v21)

            XYZ[a0] = p1 + \
              r01*(vu*sn_ang*cs_tor + vp*sn_ang*sn_tor - v21*cs_ang)
        return XYZ
