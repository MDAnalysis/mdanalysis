# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2 or any higher version
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

r"""
Default Guesser
================
.. _DefaultGuesser:

DefaultGuesser is a generic guesser class that has basic guessing methods.
This class is a general purpose guesser that can be used with most topologies,
but being generic makes it the less accurate among all guessers.





Classes
-------

.. autoclass:: DefaultGuesser
   :members:
   :inherited-members:

"""
from .base import GuesserBase
import numpy as np
import warnings
import math

import re

from ..lib import distances
from . import tables


class DefaultGuesser(GuesserBase):
    """
    This guesser holds generic methods (not directed to specific contexts) for
    guessing different topology attribute. It has the same methods which where
    originally found in Topology.guesser.py. The attributes that can be
    guessed by this class are:
    * masses
    * types
    * elements
    * angles
    * dihedrals
    * bonds
    * improper dihedrals
    * aromaticities

    You can use this guesser either directly through an instance, or through
    the :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttrs` method.

    Examples
    --------
    to guess bonds for a universe::

        import MDAnalysis as mda
        from MDAnalysisTests.datafiles import two_water_gro

        u = mda.Universe(two_water_gro, context='default', to_guess=['bonds'])

    .. versionadded:: 2.8.0

    """
    context = 'default'

    def __init__(self, universe, **kwargs):
        super().__init__(universe, **kwargs)
        self._guesser_methods = {
            'masses': self.guess_masses,
            'types': self.guess_types,
            'elements': self.guess_types,
            'bonds': self.guess_bonds,
            'angles': self.guess_angles,
            'dihedrals': self.guess_dihedrals,
            'impropers': self.guess_improper_dihedrals,
            'aromaticities': self.guess_aromaticities,
        }

    def guess_masses(self, atom_types=None, indices_to_guess=None):
        """Guess the mass of many atoms based upon their type.
        For guessing masses through Universe.guess_TopologyAttrs():
        First try to guess masses from atom elements, if not available,
        try to guess masses from types and if not available, try to guess
        types.

        Parameters
        ----------
        atom_types
          Atom types/elements to guess masses from
        indices_to_guess (optional)
          Mask array for partially guess masses for certain atoms

        Returns
        -------
        atom_masses : np.ndarray dtype float64

        Raises
        ------
        :exc:`ValueError`
            If there are no atom types or elements to guess mass from.

        """
        if atom_types is None:
            try:
                atom_types = self._universe.atoms.elements
            except AttributeError:
                try:
                    atom_types = self._universe.atoms.types
                except AttributeError:
                    try:
                        atom_types = self.guess_types(
                            atom_types=self._universe.atoms.names)
                    except ValueError:
                        raise ValueError(
                            "there is no reference attributes"
                            " (elements, types, or names)"
                            " in this universe to guess mass from")

        if indices_to_guess is not None:
            atom_types = atom_types[indices_to_guess]

        masses = np.array([self.get_atom_mass(atom)
                           for atom in atom_types], dtype=np.float64)
        return masses

    def get_atom_mass(self, element):
        """Return the atomic mass in u for *element*.
        Masses are looked up in :data:`MDAnalysis.guesser.tables.masses`.

        .. Warning:: Until version 3.0.0 unknown masses are set to 0.0

        """
        try:
            return tables.masses[element]
        except KeyError:
            try:
                return tables.masses[element.upper()]
            except KeyError:
                warnings.warn(
                    "Unknown masses are set to 0.0 for current version, "
                    "this will be deprecated in version 3.0.0 and replaced by"
                    " Masse's no_value_label (np.nan)",
                    PendingDeprecationWarning)
                return 0.0

    def guess_atom_mass(self, atomname):
        """Guess a mass based on the atom name.

        :func:`guess_atom_element` is used to determine the kind of atom.

        .. warning:: Until version 3.0.0 anything not recognized is simply
           set to 0.0; if you rely on the masses you might want to double-check.
        """
        return self.get_atom_mass(self.guess_atom_element(atomname))

    def guess_types(self, atom_types=None, masses=None, indices_to_guess=None):
        """Guess the atom type of many atoms based on atom name

        Parameters
        ----------
        atom_types (optional)
          atoms names if types guessing is desired to be from names
        masses (optional)
          atoms masses if types guessing is desired to be from masses
        indices_to_guess (optional)
          Mask array for partially guess types for certain atoms

        Returns
        -------
        atom_types : np.ndarray dtype object

        Raises
        ------
        :exc:`ValueError`
           If there is no names or masses to guess types from.

        """
        if atom_types is None and masses is None:
            try:
                atom_types = self._universe.atoms.names
            except AttributeError:
                try:
                    # if names is not available we can guess types from masses
                    masses = self._universe.atoms.masses
                except AttributeError:
                    raise ValueError(
                        "there is no reference attributes in this universe"
                        "to guess types from")

        if masses is not None:
            if indices_to_guess is not None:
                masses = masses[indices_to_guess]

            return np.array([self.guess_element_from_mass(mass)
                             for mass in masses], dtype=object)

        else:
            if indices_to_guess is not None:
                atom_types = atom_types[indices_to_guess]

            return np.array([self.guess_atom_element(atom)
                            for atom in atom_types], dtype=object)

    def guess_atom_element(self, atomname):
        """Guess the element of the atom from the name.

        Looks in dict to see if element is found, otherwise it uses the first
        character in the atomname. The table comes from CHARMM and AMBER atom
        types, where the first character is not sufficient to determine the
        atom type. Some GROMOS ions have also been added.

        .. Warning: The translation table is incomplete.
           This will probably result in some mistakes,
           but it still better than nothing!

        See Also
        --------
        :func:`guess_atom_type`
        :mod:`MDAnalysis.guesser.tables`
        """
        NUMBERS = re.compile(r'[0-9]')  # match numbers
        SYMBOLS = re.compile(r'[*+-]')  # match *, +, -
        if atomname == '':
            return ''
        try:
            return tables.atomelements[atomname.upper()]
        except KeyError:
            # strip symbols and numbers
            no_symbols = re.sub(SYMBOLS, '', atomname)
            name = re.sub(NUMBERS, '', no_symbols).upper()

            # just in case
            if name in tables.atomelements:
                return tables.atomelements[name]

            while name:
                if name in tables.elements:
                    return name
                if name[:-1] in tables.elements:
                    return name[:-1]
                if name[1:] in tables.elements:
                    return name[1:]
                if len(name) <= 2:
                    return name[0]
                name = name[:-1]  # probably element is on left not right

            # if it's numbers
            return no_symbols

    def guess_element_from_mass(self, atom_mass):
        """Guess the element of the atom from the mass.

        Compare mass to the known element masses with
        a percision of five decimal places. If masses
        matched, return the corresponding element,
        else return empty string.
        NB: guessing elements from masses fail
        with elements that doesn't have unique mass.

        See Also
        --------
        :func:`guess_atom_type`
        :mod:`MDAnalysis.guesser.tables`
        """
        for element, mass in tables.masses.items():
            if math.isclose(mass, atom_mass, rel_tol=1e-5):
                # skipping masses that represent more than one element
                # neglect it better that guessing wrong type
                if atom_mass == 247 or atom_mass == 262:
                    return ''
                else:
                    return element.upper()
        return ''

    def guess_bonds(self, atoms=None, coords=None):
        r"""Guess if bonds exist between two atoms based on their distance.

        Bond between two atoms is created, if the two atoms are within

        .. math::

              d < f \cdot (R_1 + R_2)

        of each other, where :math:`R_1` and :math:`R_2` are the VdW radii
        of the atoms and :math:`f` is an ad-hoc *fudge_factor*. This is
        the `same algorithm that VMD uses`_.

        Parameters
        ----------
        atoms : AtomGroup
             atoms for which bonds should be guessed
        fudge_factor : float, optional
            The factor by which atoms must overlap eachother to be considered a
            bond. Larger values will increase the number of bonds found. [0.55]
        vdwradii : dict, optional
            To supply custom vdwradii for atoms in the algorithm. Must be a
            dict of format {type:radii}. The default table of van der Waals
            radii is hard-coded as :data:`MDAnalysis.guesser.tables.vdwradii`.
            Any user defined vdwradii passed as an argument will supercede the
            table values. [``None``]
        lower_bound : float, optional
            The minimum bond length. All bonds found shorter than this length
            will be ignored. This is useful for parsing PDB with altloc records
            where atoms with altloc A and B maybe very close together and
            there should be no chemical bond between them. [0.1]
        box : array_like, optional
            Bonds are found using a distance search, if unit cell information
            is given, periodic boundary conditions will be considered in the
            distance search. [``None``]

        Returns
        -------
        list
            List of tuples suitable for use in Universe topology building.

        Warnings
        --------
        No check is done after the bonds are guessed to see if Lewis
        structure is correct. This is wrong and will burn somebody.

        Raises
        ------
        :exc:`ValueError` 
           If inputs are malformed or `vdwradii` data is missing.


        .. _`same algorithm that VMD uses`:
           http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node26.html

        """
        if atoms is None:
            atoms = self._universe.atoms

        if coords is None:
            coords = self._universe.atoms.positions

        if len(atoms) != len(coords):
            raise ValueError("'atoms' and 'coord' must be the same length")

        fudge_factor = self._kwargs.get('fudge_factor', 0.55)

        # so I don't permanently change it
        vdwradii = tables.vdwradii.copy()
        user_vdwradii = self._kwargs.get('vdwradii', None)
        # this should make algo use their values over defaults
        if user_vdwradii:
            vdwradii.update(user_vdwradii)

        # Try using types, then elements
        if hasattr(atoms, 'types'):
            atomtypes = atoms.types
        else:
            atomtypes = self.guess_types(atom_types=atoms.names)

        # check that all types have a defined vdw
        if not all(val in vdwradii for val in set(atomtypes)):
            raise ValueError(("vdw radii for types: " +
                              ", ".join([t for t in set(atomtypes) if
                                         t not in vdwradii]) +
                              ". These can be defined manually using the" +
                              f" keyword 'vdwradii'"))

        lower_bound = self._kwargs.get('lower_bound', 0.1)

        box = self._kwargs.get('box', None)

        if box is not None:
            box = np.asarray(box)

        # to speed up checking, calculate what the largest possible bond
        # atom that would warrant attention.
        # then use this to quickly mask distance results later
        max_vdw = max([vdwradii[t] for t in atomtypes])

        bonds = []

        pairs, dist = distances.self_capped_distance(coords,
                                                     max_cutoff=2.0 * max_vdw,
                                                     min_cutoff=lower_bound,
                                                     box=box)
        for idx, (i, j) in enumerate(pairs):
            d = (vdwradii[atomtypes[i]] +
                 vdwradii[atomtypes[j]]) * fudge_factor
            if (dist[idx] < d):
                bonds.append((atoms[i].index, atoms[j].index))
        return tuple(bonds)

    def guess_angles(self, bonds=None):
        """Given a list of Bonds, find all angles that exist between atoms.

        Works by assuming that if atoms 1 & 2 are bonded, and 2 & 3 are bonded,
        then (1,2,3) must be an angle.

        Parameters
        ----------
        bonds : Bonds
             from which angles should be guessed

        Returns
        -------
        list of tuples
            List of tuples defining the angles.
            Suitable for use in u._topology


        See Also
        --------
        :meth:`guess_bonds`

      """
        from ..core.universe import Universe

        angles_found = set()
   
        if bonds is None:
            if hasattr(self._universe.atoms, 'bonds'):
                bonds = self._universe.atoms.bonds
            else:
                temp_u = Universe.empty(n_atoms=len(self._universe.atoms))
                temp_u.add_bonds(self.guess_bonds(
                    self._universe.atoms, self._universe.atoms.positions))
                bonds = temp_u.atoms.bonds

        for b in bonds:
            for atom in b:
                other_a = b.partner(atom)  # who's my friend currently in Bond
                for other_b in atom.bonds:
                    if other_b != b:  # if not the same bond I start as
                        third_a = other_b.partner(atom)
                        desc = tuple(
                            [other_a.index, atom.index, third_a.index])
                        # first index always less than last
                        if desc[0] > desc[-1]:
                            desc = desc[::-1]
                        angles_found.add(desc)

        return tuple(angles_found)

    def guess_dihedrals(self, angles=None):
        """Given a list of Angles, find all dihedrals that exist between atoms.

        Works by assuming that if (1,2,3) is an angle, and 3 & 4 are bonded,
        then (1,2,3,4) must be a dihedral.

        Parameters
        ----------
        angles : Angles
             from which dihedrals should be guessed

        Returns
        -------
        list of tuples
            List of tuples defining the dihedrals.
            Suitable for use in u._topology

        """
        from ..core.universe import Universe

        if angles is None:
            if hasattr(self._universe.atoms, 'angles'):
                angles = self._universe.atoms.angles

            else:
                temp_u = Universe.empty(n_atoms=len(self._universe.atoms))

                temp_u.add_bonds(self.guess_bonds(
                    self._universe.atoms, self._universe.atoms.positions))
     
                temp_u.add_angles(self.guess_angles(temp_u.atoms.bonds))

                angles = temp_u.atoms.angles

        dihedrals_found = set()

        for b in angles:
            a_tup = tuple([a.index for a in b])  # angle as tuple of numbers
            # if searching with b[0], want tuple of (b[2], b[1], b[0], +new)
            # search the first and last atom of each angle
            for atom, prefix in zip([b.atoms[0], b.atoms[-1]],
                                    [a_tup[::-1], a_tup]):
                for other_b in atom.bonds:
                    if not other_b.partner(atom) in b:
                        third_a = other_b.partner(atom)
                        desc = prefix + (third_a.index,)
                        if desc[0] > desc[-1]:
                            desc = desc[::-1]
                        dihedrals_found.add(desc)

        return tuple(dihedrals_found)

    def guess_improper_dihedrals(self, angles=None):
        """Given a list of Angles, find all improper dihedrals
        that exist between atoms.

        Works by assuming that if (1,2,3) is an angle, and 2 & 4 are bonded,
        then (2, 1, 3, 4) must be an improper dihedral.
        ie the improper dihedral is the angle between the planes formed by
        (1, 2, 3) and (1, 3, 4)

        Returns
        -------
            List of tuples defining the improper dihedrals.
            Suitable for use in u._topology

        """

        from ..core.universe import Universe

        if angles is None:
            if hasattr(self._universe.atoms, 'angles'):
                angles = self._universe.atoms.angles

            else:
                temp_u = Universe.empty(n_atoms=len(self._universe.atoms))

                temp_u.add_bonds(self.guess_bonds(
                    self._universe.atoms, self._universe.atoms.positions))

                temp_u.add_angles(self.guess_angles(temp_u.atoms.bonds))

                angles = temp_u.atoms.angles

        dihedrals_found = set()

        for b in angles:
            atom = b[1]  # select middle atom in angle
            # start of improper tuple
            a_tup = tuple([b[a].index for a in [1, 2, 0]])
            # if searching with b[1], want tuple of (b[1], b[2], b[0], +new)
            # search the first and last atom of each angle
            for other_b in atom.bonds:
                other_atom = other_b.partner(atom)
                # if this atom isn't in the angle I started with
                if other_atom not in b:
                    desc = a_tup + (other_atom.index,)
                    if desc[0] > desc[-1]:
                        desc = desc[::-1]
                    dihedrals_found.add(desc)

        return tuple(dihedrals_found)

    def guess_atom_charge(self, atoms):
        """Guess atom charge from the name.

        .. Warning:: Not implemented; simply returns 0.
        """
        # TODO: do something slightly smarter, at least use name/element
        return 0.0

    def guess_aromaticities(self, atomgroup=None):
        """Guess aromaticity of atoms using RDKit

        Returns
        -------
        aromaticities : numpy.ndarray
            Array of boolean values for the aromaticity of each atom

        """
        if atomgroup is None:
            atomgroup = self._universe.atoms

        mol = atomgroup.convert_to("RDKIT")
        return np.array([atom.GetIsAromatic() for atom in mol.GetAtoms()])

    def guess_gasteiger_charges(self, atomgroup):
        """Guess Gasteiger partial charges using RDKit

        Parameters
        ----------
        atomgroup : mda.core.groups.AtomGroup
            Atoms for which the charges will be guessed

        Returns
        -------
        charges : numpy.ndarray
            Array of float values representing the charge of each atom

        """

        mol = atomgroup.convert_to("RDKIT")
        from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
        ComputeGasteigerCharges(mol, throwOnParamFailure=True)
        return np.array([atom.GetDoubleProp("_GasteigerCharge")
                         for atom in mol.GetAtoms()],
                        dtype=np.float32)
