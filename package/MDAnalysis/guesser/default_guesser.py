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
but being generic makes it the least accurate among all guessers.


Guessing behavior
-----------------

This section describes how each attribute is guessed by the DefaultGuesser.

Masses
~~~~~~

We first attempt to look up the mass of an atom based on its element if the
element TopologyAttr is available. If not, we attempt to lookup the mass based
on the atom type (``type``) TopologyAttr. If neither of these is available, we
attempt to guess the atom type based on the atom name (``name``) and then
lookup the mass based on the guessed atom type.


Types
~~~~~

We attempt to guess the atom type based on the atom name (``name``).
The name is first stripped of any numbers and symbols, and then looked up in
the :data:`MDAnalysis.guesser.tables.atomelements` table. If the name is not
found, we continue checking variations of the name following the logic in
:meth:`DefaultGuesser.guess_atom_element`. Ultimately, if no match is found,
the first character of the stripped name is returned.

Elements
~~~~~~~~

This follows the same method as guessing atom types.


Bonds
~~~~~

Bonds are guessed based on the distance between atoms.
See :meth:`DefaultGuesser.guess_bonds` for more details.

Angles
~~~~~~

Angles are guessed based on the bonds between atoms.
See :meth:`DefaultGuesser.guess_angles` for more details.

Dihedrals
~~~~~~~~~

Dihedrals are guessed based on the angles between atoms.
See :meth:`DefaultGuesser.guess_dihedrals` for more details.

Improper Dihedrals
~~~~~~~~~~~~~~~~~~

Improper dihedrals are guessed based on the angles between atoms.
See :meth:`DefaultGuesser.guess_improper_dihedrals` for more details.

Aromaticities
~~~~~~~~~~~~~

Aromaticity is guessed using RDKit's GetIsAromatic method.
See :meth:`DefaultGuesser.guess_aromaticities` for more details.




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

from ..exceptions import NoDataError
from ..lib import distances
from . import tables


class DefaultGuesser(GuesserBase):
    """
    This guesser holds generic methods (not directed to specific contexts) for
    guessing different topology attribute. It has the same methods which where
    originally found in Topology.guesser.py. The attributes that can be
    guessed by this class are:
     - masses
     - types
     - elements
     - angles
     - dihedrals
     - bonds
     - improper dihedrals
     - aromaticities

    You can use this guesser either directly through an instance, or through
    the :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttrs` method.

    Parameters
    ----------
    universe : Universe
        The Universe to apply the guesser on
    box : np.ndarray, optional
        The box of the Universe. This is used for bond guessing.
    vdwradii : dict, optional
        Dict relating atom types: vdw radii. This is used for bond guessing
    fudge_factor : float, optional
        The factor by which atoms must overlap each other to be considered
        a bond.  Larger values will increase the number of bonds found. [0.55]
    lower_bound : float, optional
        The minimum bond length. All bonds found shorter than this length
        will be ignored. This is useful for parsing PDB with altloc records
        where atoms with altloc A and B may be very close together and
        there should be no chemical bond between them. [0.1]

    Examples
    --------
    to guess bonds for a universe::

        import MDAnalysis as mda
        from MDAnalysisTests.datafiles import two_water_gro

        u = mda.Universe(two_water_gro, context='default', to_guess=['bonds'])

    .. versionadded:: 2.8.0

    """

    context = "default"

    def __init__(
        self,
        universe,
        box=None,
        vdwradii=None,
        fudge_factor=0.55,
        lower_bound=0.1,
        **kwargs,
    ):
        super().__init__(
            universe,
            box=box,
            vdwradii=vdwradii,
            fudge_factor=fudge_factor,
            lower_bound=lower_bound,
            **kwargs,
        )
        self._guesser_methods = {
            "masses": self.guess_masses,
            "types": self.guess_types,
            "elements": self.guess_types,
            "bonds": self.guess_bonds,
            "angles": self.guess_angles,
            "dihedrals": self.guess_dihedrals,
            "impropers": self.guess_improper_dihedrals,
            "aromaticities": self.guess_aromaticities,
        }

    def guess_masses(self, atom_types=None, indices_to_guess=None):
        """Guess the mass of many atoms based upon their type.
        For guessing masses through :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttrs`:

        First try to guess masses from atom elements, if not available,
        try to guess masses from types and if not available, try to guess
        types.

        Parameters
        ----------
        atom_types : Optional[np.ndarray]
          Atom types/elements to guess masses from
        indices_to_guess : Optional[np.ndarray]
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
            except NoDataError:
                try:
                    atom_types = self._universe.atoms.types
                except NoDataError:
                    try:
                        atom_types = self.guess_types(
                            atom_types=self._universe.atoms.names
                        )
                    except NoDataError:
                        raise NoDataError(
                            "there is no reference attributes"
                            " (elements, types, or names)"
                            " in this universe to guess mass from"
                        ) from None

        if indices_to_guess is not None:
            atom_types = atom_types[indices_to_guess]

        masses = np.array(
            [self.get_atom_mass(atom) for atom in atom_types], dtype=np.float64
        )
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
                    PendingDeprecationWarning,
                )
                return 0.0

    def guess_atom_mass(self, atomname):
        """Guess a mass based on the atom name.

        :func:`guess_atom_element` is used to determine the kind of atom.

        .. warning:: Until version 3.0.0 anything not recognized is simply
           set to 0.0; if you rely on the masses you might want to double-check.
        """
        return self.get_atom_mass(self.guess_atom_element(atomname))

    def guess_types(self, atom_types=None, indices_to_guess=None):
        """Guess the atom type of many atoms based on atom name

        Parameters
        ----------
        atom_types (optional)
          atoms names if types guessing is desired to be from names
        indices_to_guess (optional)
          Mask array for partially guess types for certain atoms

        Returns
        -------
        atom_types : np.ndarray dtype object

        Raises
        ------
        :exc:`ValueError`
           If there is no names to guess types from.

        """
        if atom_types is None:
            try:
                atom_types = self._universe.atoms.names
            except NoDataError:
                raise NoDataError(
                    "there is no reference attributes in this universe "
                    "to guess types from"
                ) from None

        if indices_to_guess is not None:
            atom_types = atom_types[indices_to_guess]

        # guess each distinct name once
        guesses = {}
        for atom in atom_types:
            if atom not in guesses:
                guesses[atom] = self.guess_atom_element(atom)
        return np.array([guesses[atom] for atom in atom_types], dtype=object)

    def guess_atom_element(self, atomname):
        """Guess the element of the atom from the name.

        First all numbers and symbols are stripped from the name.
        Then the name is looked up in the
        :data:`MDAnalysis.guesser.tables.atomelements` table.
        If the name is not found, we remove the last character or
        first character from the name and check the table for both,
        with a preference for removing the last character. If the name is
        still not found, we iteratively continue to remove the last character
        or first character until we find a match. If ultimately no match
        is found, the first character of the stripped name is returned.

        If the input name is an empty string, an empty string is returned.

        The table comes from CHARMM and AMBER atom
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
        NUMBERS = re.compile(r"[0-9]")  # match numbers
        SYMBOLS = re.compile(r"[*+-]")  # match *, +, -
        if atomname == "":
            return ""
        try:
            return tables.atomelements[atomname.upper()]
        except KeyError:
            # strip symbols and numbers
            no_symbols = re.sub(SYMBOLS, "", atomname)
            name = re.sub(NUMBERS, "", no_symbols).upper()

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
        atoms: AtomGroup
            atoms for which bonds should be guessed
        coords: np.ndarray, optional
            coordinates of the atoms. If not provided, the coordinates
            of the ``atoms`` in the universe are used.

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

        fudge_factor = self._kwargs.get("fudge_factor", 0.55)

        # so I don't permanently change it
        vdwradii = tables.vdwradii.copy()
        user_vdwradii = self._kwargs.get("vdwradii", None)
        # this should make algo use their values over defaults
        if user_vdwradii:
            vdwradii.update(user_vdwradii)

        # Try using types, then elements
        if hasattr(atoms, "types"):
            atomtypes = atoms.types
        else:
            atomtypes = self.guess_types(atom_types=atoms.names)

        # vdwradii keys are uppercase, so normalize atomtypes to match
        atomtypes = np.char.upper(np.asarray(atomtypes, dtype=str))

        # check that all types have a defined vdw
        if not all(val in vdwradii for val in set(atomtypes)):
            raise ValueError(
                (
                    "vdw radii for types: "
                    + ", ".join(
                        [t for t in set(atomtypes) if t not in vdwradii]
                    )
                    + ". These can be defined manually using the"
                    + f" keyword 'vdwradii'"
                )
            )

        lower_bound = self._kwargs.get("lower_bound", 0.1)

        box = self._kwargs.get("box", None)

        if box is not None:
            box = np.asarray(box)

        # to speed up checking, calculate what the largest possible bond
        # atom that would warrant attention.
        # then use this to quickly mask distance results later
        unique_types, type_idx = np.unique(atomtypes, return_inverse=True)
        radii = np.array([vdwradii[t] for t in unique_types])[type_idx]
        max_vdw = radii.max()

        pairs, dist = distances.self_capped_distance(
            coords, max_cutoff=2.0 * max_vdw, min_cutoff=lower_bound, box=box
        )
        bonded = (
            dist < (radii[pairs[:, 0]] + radii[pairs[:, 1]]) * fudge_factor
        )
        return tuple(map(tuple, atoms.ix[pairs[bonded]].tolist()))

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

        if bonds is None:
            if hasattr(self._universe.atoms, "bonds"):
                bonds = self._universe.atoms.bonds
            else:
                temp_u = Universe.empty(n_atoms=len(self._universe.atoms))
                temp_u.add_bonds(
                    self.guess_bonds(
                        self._universe.atoms, self._universe.atoms.positions
                    )
                )
                bonds = temp_u.atoms.bonds

        bix, neighbours = _connection_indices(bonds, 2)
        if len(bix) == 0:
            return ()
        # each bond (i, j) is extended at both ends: (j, i, k) and (i, j, k)
        centre = np.concatenate([bix[:, 0], bix[:, 1]])
        other = np.concatenate([bix[:, 1], bix[:, 0]])
        row, third = neighbours(centre)
        keep = third != other[row]  # not the bond we started from
        row, third = row[keep], third[keep]
        found = np.column_stack([other[row], centre[row], third])
        return _unique_connections(found)

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
            if hasattr(self._universe.atoms, "angles"):
                angles = self._universe.atoms.angles

            else:
                temp_u = Universe.empty(n_atoms=len(self._universe.atoms))

                temp_u.add_bonds(
                    self.guess_bonds(
                        self._universe.atoms, self._universe.atoms.positions
                    )
                )

                temp_u.add_angles(self.guess_angles(temp_u.atoms.bonds))

                angles = temp_u.atoms.angles

        aix, neighbours = _connection_indices(angles, 3)
        if len(aix) == 0:
            return ()
        # extend each angle at its first atom, (a2, a1, a0, x), and at its
        # last atom, (a0, a1, a2, x)
        prefix = np.concatenate([aix[:, ::-1], aix])
        row, fourth = neighbours(prefix[:, -1])
        prefix = prefix[row]
        keep = (prefix != fourth[:, None]).all(axis=1)
        found = np.column_stack([prefix[keep], fourth[keep]])
        return _unique_connections(found)

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
            if hasattr(self._universe.atoms, "angles"):
                angles = self._universe.atoms.angles

            else:
                temp_u = Universe.empty(n_atoms=len(self._universe.atoms))

                temp_u.add_bonds(
                    self.guess_bonds(
                        self._universe.atoms, self._universe.atoms.positions
                    )
                )

                temp_u.add_angles(self.guess_angles(temp_u.atoms.bonds))

                angles = temp_u.atoms.angles

        aix, neighbours = _connection_indices(angles, 3)
        if len(aix) == 0:
            return ()
        # extend each angle (a0, a1, a2) at its middle atom: (a1, a2, a0, x)
        prefix = aix[:, [1, 2, 0]]
        row, fourth = neighbours(prefix[:, 0])
        prefix = prefix[row]
        keep = (prefix != fourth[:, None]).all(axis=1)
        found = np.column_stack([prefix[keep], fourth[keep]])
        return _unique_connections(found)

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
        return np.array(
            [
                atom.GetDoubleProp("_GasteigerCharge")
                for atom in mol.GetAtoms()
            ],
            dtype=np.float32,
        )


def _connection_indices(connections, n_atoms):
    """Atom indices of `connections` and a bonded-neighbour lookup.

    Parameters
    ----------
    connections : TopologyGroup or iterable of TopologyObject
        bonds or angles to extend
    n_atoms : int
        number of atoms per connection

    Returns
    -------
    indices : numpy.ndarray
        ``(n_connections, n_atoms)`` atom indices
    neighbours : callable
        ``neighbours(atoms)`` returns ``(row, neighbour)`` arrays: for every
        bond of ``atoms[row]`` in the whole universe, the other atom
    """
    try:
        indices, universe = connections.indices, connections.universe
    except AttributeError:
        connections = list(connections)
        if not connections:
            return np.empty((0, n_atoms), dtype=np.intp), None
        indices = np.array([c.indices for c in connections])
        universe = connections[0].universe
    if len(indices) == 0:
        return np.empty((0, n_atoms), dtype=np.intp), None
    indices = np.asarray(indices, dtype=np.intp)

    try:
        values = universe._topology.bonds.values
    except AttributeError:
        raise NoDataError(
            "This Universe does not contain bonds information"
        ) from None
    values = np.array(values, dtype=np.intp).reshape(-1, 2)
    # sparse (CSR) neighbour list over all bonds in the universe
    src = np.concatenate([values[:, 0], values[:, 1]])
    dst = np.concatenate([values[:, 1], values[:, 0]])
    nbrs = dst[np.argsort(src, kind="stable")]
    degree = np.bincount(src, minlength=universe._topology.n_atoms)
    start = np.cumsum(degree) - degree

    def neighbours(atoms):
        deg = degree[atoms]
        row = np.repeat(np.arange(len(atoms)), deg)
        offset = np.arange(len(row)) - np.repeat(np.cumsum(deg) - deg, deg)
        return row, nbrs[start[atoms][row] + offset]

    return indices, neighbours


def _unique_connections(found):
    """Orient each connection so the first index is lower than the last
    and return the distinct connections as a tuple of tuples."""
    flip = found[:, 0] > found[:, -1]
    found[flip] = found[flip, ::-1]
    return tuple(set(map(tuple, found.tolist())))
