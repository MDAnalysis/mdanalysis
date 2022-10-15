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
from .base import GuesserBase
import re
import numpy as np
from . import tables, pdb_tables
import warnings


class PDBGuesser(GuesserBase):
    context = 'pdb'

    def __init__(self, universe, **kwargs):
        super().__init__(universe, **kwargs)
        self._guesser_methods = {'masses': self.guess_masses,
                                 'types': self.guess_types,
                                 'elements': self.guess_types,

                                 }

    def guess_types(self, atoms_index=None):          

        names = []
        residues =[]

        try:
            names = self._universe.atoms.names
        except AttributeError:
            raise ValueError("there is no name attribute in this universe"
                                 "to guess types from")
             

        try:
            residues = self._universe.atoms.residues.resnames
        except AttributeError:
            raise ValueError("there is no residue name attribute in this universe"
                                 "to guess types from")

        if not atoms_index:
           atoms_index = list(range(len(names)))

        elements = []
        failed_guessing = set()
        failiur_count = 0
        SYMBOLS = re.compile(r'[*\'\s\d\"]')  # match *, ', " space charachter, and digits
        for n in atoms_index:
            if(SYMBOLS.match(names[n][0])):
                if names[n][0] == ' ' and names[n][2].islower():
                    if names[n][1:2].upper() in tables.masses:
                        # concord file ??
                        elements.append(names[n][1:2].upper())
                    elif names[n][1].upper() in tables.masses:
                        elements.append(names[n][1].upper())
                    else:
                        elements.append('')
                        failed_guessing.add(names[n])
                        failiur_count += 1
                elif names[n][1] in tables.masses:
                    elements.append(names[n][1])

                elif names[n][2] in ['H', 'C', 'N', 'O', 'P', 'S']:
                    elements.append(names[n][2])
                else:
                    elements.append('')
                    failed_guessing.add(names[n])
                    failiur_count += 1
            else:
                if names[n][0] == 'H' and residues[n] in pdb_tables.standard_groups:
                    elements.append(names[n][0])
                elif names[n][:1].upper() in tables.masses:
                    elements.append(names[n][:1].upper())
                elif names[n][0] == 'H':
                    elements.append(names[n][0])
                else:
                    elements.append('')
                    failed_guessing.add(names[n])
                    failiur_count += 1

        if failed_guessing:
            sucessful = len(elements) - failiur_count
            warnings.warn(f'Sucessfuly guessed types for {sucessful}/{len(elements)} atoms. '
                          f'Failed to guess types for the following atoms: {failed_guessing}')

        return np.array(elements, dtype=object)

    def guess_masses(self, atoms=None):
        """Guess the mass of many atoms based upon their type

        Parameters
        ----------
        atoms
          Atom types/elements to guess masses from

        Returns
        -------
        atom_masses : np.ndarray dtype float64
        """
        atom_types = None
        if atoms is not None:
            atom_types = atoms
        elif hasattr(self._universe.atoms, 'elements'):
            atom_types = self._universe.atoms.elements
        elif hasattr(self._universe.atoms, 'types'):
            atom_types = self._universe.atoms.types

        else:
            try:
                self._universe.guess_TopologyAttributes(
                    self.context, ['elements'])
                atom_types = self._universe.atoms.types
            except AttributeError:
                raise ValueError("there is no reference attributes in this universe"
                                 "to guess masses from")

        masses = []
        failed = set()
        failiur_count = 0
        for a in atom_types:
            try:
                masses.append(tables.masses[a.upper()])
            except KeyError:
                masses.append(0)
                failed.add(a)
                failiur_count += 1

        if failed:
            sucessful = len(masses) - failiur_count
            warnings.warn(f'Sucessfuly guessed masses for {sucessful}/{len(masses)} atoms. '
                          f'Failed to guess masses for the following atoms: {failed}')

            warnings.warn("Unknown masses are set to 0.0 for current version, this will be"
                          "depracated in version 3.0.0 and replaced by Masse's no_value_label (np.nan)", PendingDeprecationWarning)

        return np.array(masses, dtype=np.float64)
