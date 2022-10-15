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

"""
PDB Guesser
=========================================================================

This guesser class is tailored to be used with PDB standard files data. 
to be continued after settling the funtionalities

Elements guessing
PDB has a well-defined formatt for names, which from it we can get the atomic symbol easily.

Names always found in columns 13-16. The first two charachters represents the atom symbol. If
the symbol consists of one characters then the first character is blank. Then at the third character
comes the remoteness indicator code for amino acids residues [A, B, G, D, E, Z, H]. The last character
is a branching factor if needed.

The above rules is the standard rules but there is some exceptions to it:
    1- CONCORD generated PDB files shift the atom name columns and generate the second character of the atom symbol
    as a lower case, so if the third charachter is a lower case, then the atom symbols is taken from
    the second and third characters
    2- If the first character is blank and the second character is not recognized as atomic symbol,
    we check if the third character contains "H", "C", "N", "O", "P" or "S", then it is considered
    the atomic symbol
    3-  If the first character is a digit, \”, \’, or *, then the second character is the atomic symbol
    4- If the fisrt character in 'H' and the residudes is a standard amino acid, nucliec acid or known
    hetero groups (found in pdb_tables.py), then the atom element is 'H'
    5- If the first two characters are not recognized as atomic symmbol and the first character is 'H',
    then the element is considered to be H


Classes
-------

.. autoclass:: PDBGuesser
   :members:
   :inherited-members:

"""


class PDBGuesser(GuesserBase):
    context = 'pdb'

    def __init__(self, universe, **kwargs):
        super().__init__(universe, **kwargs)
        self._guesser_methods = {'masses': self.guess_masses,
                                 'types': self.guess_types,
                                 'elements': self.guess_types,

                                 }

    def guess_types(self, atoms_index=None, partial_guess=None):

        names = []
        residues = []

        try:
            names = self._universe.atoms.names
        except AttributeError:
            raise ValueError("there is no name attribute in this universe"
                             "to guess types from")

        if not hasattr(self._universe.atoms.residues, 'resnames'):
            raise ValueError("there is no residue name attribute in this universe"
                             "to guess types from")

        atoms_indices = partial_guess if partial_guess else list(
            range(len(names)))

        elements = []
        failed_guessing = set()
        failiur_count = 0
        # match *, ', " space charachter, and digits
        SYMBOLS = re.compile(r'[*\'\s\d\"]')
        for n in atoms_indices:
            if(SYMBOLS.match(names[n][0])):
                # concord files special case
                if names[n][0] == ' ' and names[n][2].islower():
                    if names[n][1:3].upper() in tables.masses:
                        elements.append(names[n][1:3].upper())
                elif names[n][1] in tables.masses:
                    elements.append(names[n][1])

                elif names[n][2] in ['H', 'C', 'N', 'O', 'P', 'S']:
                    elements.append(names[n][2])
                else:
                    elements.append('')
                    failed_guessing.add(names[n])
                    failiur_count += 1
            else:
                if names[n][0] == 'H' and self._universe.atoms[n].residue.resname in pdb_tables.standard_groups:
                    elements.append(names[n][0])
                elif names[n][:2].upper() in tables.masses:
                    elements.append(names[n][:2].upper())
                elif names[n][1].upper() in tables.masses:
                    elements.append(names[n][1])
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
