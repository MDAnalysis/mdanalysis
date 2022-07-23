from .base import GuesserBase
import numpy as np
import warnings

import re

from ..lib import distances
from . import tables


class DefaultGuesser(GuesserBase):
    context = 'default'
 
    def __init__(self, atoms):
        super().__init__(atoms)
        self._guess = {'mass': self.guess_masses,
                       'type': self.guess_types}
        self._rank = {'mass': 2,
                      'type': 1
                      }

    def guess_masses(self):
        """Guess the mass of many atoms based upon their type

        Returns
        -------
        atom_masses : np.ndarray dtype float64
        """
        if hasattr(self._atoms, 'elements'):
            atom_types = self._atoms.elements
        else:
            atom_types = self._atoms.types

        self.validate_atom_types(atom_types)
        masses = np.array([self.get_atom_mass(atom_t)
                           for atom_t in atom_types], dtype=np.float64)
        return masses

    def validate_atom_types(self, atom_types):
        """Vaildates the atom types based on whether they are available
        in our tables

        Parameters
        ----------
        atom_types
          Type of each atom

        Returns
        -------
        None

        .. versionchanged:: 0.20.0
           Try uppercase atom type name as well
        """
        for atom_type in np.unique(atom_types):
            try:
                tables.masses[atom_type]
            except KeyError:
                try:
                    tables.masses[atom_type.upper()]
                except KeyError:
                    warnings.warn
                    ("Failed to guess the mass for the following atoms: {}"
                     .format(atom_type))

    def get_atom_mass(self, element):
        """Return the atomic mass in u for *element*.

        Masses are looked up in :data:`MDAnalysis.topology.tables.masses`.

        .. Warning:: Unknown masses are set to 0.0

        .. versionchanged:: 0.20.0
           Try uppercase atom type name as well
        """
        try:
            return tables.masses[element]
        except KeyError:
            try:
                return tables.masses[element.upper()]
            except KeyError:
                return 0.0

    def guess_atom_mass(self, atomname):
        """Guess a mass based on the atom name.

        :func:`guess_atom_element` is used to determine the kind of atom.

        .. warning:: Anything not recognized is simply set to 0;
        if you rely on the masses you might want to double check.
        """
        return self.get_atom_mass(self.guess_atom_element(atomname))
    
    def guess_types(self):
        """Guess the atom type of many atoms based on atom name

        Parameters
        ----------
        atom_names
          Name of each atom

        Returns
        -------
        atom_types : np.ndarray dtype object
        """
        names = self._atoms.names
        return
        np.array([self.guess_atom_element(n) for n in names], dtype=object)

    NUMBERS = re.compile(r'[0-9]')  # match numbers
    SYMBOLS = re.compile(r'[*+-]')  # match *, +, -

    def guess_atom_element(self, atomname):
        """Guess the element of the atom from the name.

        Looks in dict to see if element is found, otherwise it uses the first
        character in the atomname. The table comes from CHARMM and AMBER atom
        types, where the first character is not sufficient to determine the
        atom type. Some GROMOS ions have also been added.

        .. Warning: The translation table is incomplete.
           This will probably result in some mistakes,
           but it still bbetter than nothing!

        See Also
        --------
        :func:`guess_atom_type`
        :mod:`MDAnalysis.topology.tables`
        """
        if atomname == '':
            return ''
        try:
            return tables.atomelements[atomname.upper()]
        except KeyError:
            # strip symbols and numbers
            no_symbols = re.sub(self.SYMBOLS, '', atomname)
            name = re.sub(self.NUMBERS, '', no_symbols).upper()

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
