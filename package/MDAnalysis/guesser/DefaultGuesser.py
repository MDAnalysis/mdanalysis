from .base import GuesserBase
import numpy as np
import warnings
from . import tables


class DefaultGuesser(GuesserBase):
    context = 'default'

    def __init__(self):
        self._guess = {'mass': self.guess_masses}

    def guess_masses(self):
        """Guess the mass of many atoms based upon their type

        Returns
        -------
        atom_masses : np.ndarray dtype float64
        """
        atom_types = self._atoms.types
        self.validate_atom_types(atom_types)
        masses = np.array([self.get_atom_mass(atom_t)
                           for atom_t in atom_types], dtype=np.float64)
        return masses

    def validate_atom_types(self, atom_types):
        """Vaildates the atom types based on whether they are available in our tables

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
                    warnings.warn("Failed to guess the mass for the following atom types: {}".format(atom_type))

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

        .. warning:: Anything not recognized is simply set to 0; if you rely on the
                     masses you might want to double check.
        """
        return self.get_atom_mass(self.guess_atom_element(atomname))
