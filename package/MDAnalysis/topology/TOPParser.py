# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.mdanalysis.org
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
AMBER PRMTOP topology parser
============================

Reads an AMBER top file to build the system.

Amber keywords are turned into the following attributes:

+-----------------+----------------------+
| AMBER flag      | MDAnalysis attribute |
+-----------------+----------------------+
| ATOM_NAME       | names                |
+-----------------+----------------------+
| CHARGE          | charges              |
+-----------------+----------------------+
| ATOMIC_NUMBER   | elements             |
+-----------------+----------------------+
| MASS            | masses               |
+-----------------+----------------------+
| ATOM_TYPE_INDEX | type_indices         |
+-----------------+----------------------+
| AMBER_ATOM_TYPE | types                |
+-----------------+----------------------+
| RESIDUE_LABEL   | resnames             |
+-----------------+----------------------+
| RESIDUE_POINTER | residues             |
+-----------------+----------------------+

TODO:
  Add reading of bonds etc

.. Note::

   The Amber charge is converted to electron charges as used in
   MDAnalysis and other packages. To get back Amber charges, multiply
   by 18.2223.

.. _`PARM parameter/topology file specification`:
   http://ambermd.org/formats.html#topology

Classes
-------

.. autoclass:: TOPParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import, division

from six.moves import range, zip
import numpy as np
import functools
from math import ceil

from . import guessers
from .tables import NUMBER_TO_ELEMENT
from ..lib.util import openany, FORTRANReader
from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomtypes,
    Atomids,
    Charges,
    Elements,
    Masses,
    Resnames,
    Resids,
    Resnums,
    Segids,
    AtomAttr,
)


class TypeIndices(AtomAttr):
    """Numerical type of each Atom"""
    attrname = 'type_indices'
    singular = 'type_index'
    level = 'atom'


class TOPParser(TopologyReaderBase):
    """Reads topology information from an AMBER top file.

    Reads the following Attributes if in topology:
    - Atomnames
    - Charges
    - Masses
    - Elements
    - Atomtypes
    - Resnames
    - Type_indices

    Guesses the following attributes:
     - Elements (if not included in topology)

    The format is defined in `PARM parameter/topology file
    specification`_.  The reader tries to detect if it is a newer
    (AMBER 12?) file format by looking for the flag "ATOMIC_NUMBER".

    .. _`PARM parameter/topology file specification`:
       http://ambermd.org/formats.html#topology

   .. versionchanged:: 0.7.6
      parses both amber10 and amber12 formats
    """
    format = ['TOP', 'PRMTOP', 'PARM7']

    def parse(self):
        """Parse Amber PRMTOP topology file *filename*.

        Returns
        -------
        A MDAnalysis Topology object
        """
        # Sections that we grab as we parse the file
        sections = {
            "ATOM_NAME": (1, 20, self.parse_names, "name", 0),
            "CHARGE": (1, 5, self.parse_charges, "charge", 0),
            "ATOMIC_NUMBER": (1, 10, self.parse_elements, "elements", 0),
            "MASS": (1, 5, self.parse_masses, "mass", 0),
            "ATOM_TYPE_INDEX": (1, 10, self.parse_type_indices, "type_indices", 0),
            "AMBER_ATOM_TYPE": (1, 20, self.parse_types, "types", 0),
            "RESIDUE_LABEL": (1, 20, self.parse_resnames, "resname", 11),
            "RESIDUE_POINTER": (2, 10, self.parse_residx, "respoint", 11),
        }

        attrs = {}  # empty dict for attrs that we'll fill

        # Open and check top validity
        # Reading header info POINTERS
        with openany(self.filename) as self.topfile:
            header = self.topfile.next()
            if not header.startswith("%VE"):
                raise ValueError(
                    "{0} is not a valid TOP file. %VE Missing in header"
                    "".format(self.filename))
            title = self.topfile.next().split()
            if not (title[1] == "TITLE"):
                raise ValueError(
                    "{0} is not a valid TOP file. 'TITLE' missing in header"
                    "".format(self.filename))
            while not header.startswith('%FLAG POINTERS'):
                header = self.topfile.next()
            self.topfile.next()

            topremarks = [self.topfile.next().strip() for i in range(4)]
            sys_info = [int(k) for i in topremarks for k in i.split()]

            header = self.topfile.next()
            # grab the next section title
            next_section = header.split("%FLAG")[1].strip()

            while next_section is not None:
                try:
                    (atoms_per, per_line,
                     func, name, sect_num) = sections[next_section]
                except KeyError:
                    next_getter = self.skipper
                else:
                    num = sys_info[sect_num]
                    numlines = (num // per_line)
                    if num % per_line != 0:
                        numlines += 1

                    attrs[name] = func(atoms_per, numlines)

                    next_getter = self.topfile.next

                try:
                    line = next_getter()
                except StopIteration:
                    next_section = None
                else:
                    next_section = line.split("%FLAG")[1].strip()

        # strip out a few values to play with them
        n_atoms = len(attrs['name'])

        resptrs = attrs.pop('respoint')
        resptrs.append(n_atoms)
        residx = np.zeros(n_atoms, dtype=np.int32)
        for i, (x, y) in enumerate(zip(resptrs[:-1], resptrs[1:])):
            residx[x:y] = i

        n_res = len(attrs['resname'])

        # Guess elements if not in topology
        if not 'elements' in attrs:
            attrs['elements'] = Elements(
                guessers.guess_types(attrs['types'].values),
                guessed=True)

        # atom ids are mandatory
        attrs['atomids'] = Atomids(np.arange(n_atoms) + 1)
        attrs['resids'] = Resids(np.arange(n_res) + 1)
        attrs['resnums'] = Resnums(np.arange(n_res) + 1)
        attrs['segids'] = Segids(np.array(['SYSTEM'], dtype=object))

        top = Topology(n_atoms, n_res, 1,
                       attrs=attrs.values(),
                       atom_resindex=residx,
                       residue_segindex=None)

        return top

    def skipper(self):
        """Skip until we find the next %FLAG entry and return that"""
        line = self.topfile.next()
        while not line.startswith("%FLAG"):
            line = self.topfile.next()
        return line

    def parse_names(self, atoms_per, numlines):
        vals = self.parsesection_mapper(
            atoms_per, numlines, lambda x: x)
        attr = Atomnames(np.array(vals, dtype=object))
        return attr

    def parse_resnames(self, atoms_per, numlines):
        vals = self.parsesection_mapper(
            atoms_per, numlines, lambda x: x)
        attr = Resnames(np.array(vals, dtype=object))
        return attr

    def parse_charges(self, atoms_per, numlines):
        vals = self.parsesection_mapper(
            atoms_per, numlines, lambda x: float(x))
        charges = np.array(vals, dtype=np.float32)
        charges /= 18.2223  # to electron charge units
        attr = Charges(charges)
        return attr

    def parse_masses(self, atoms_per, numlines):
        vals = self.parsesection_mapper(
            atoms_per, numlines, lambda x: float(x))
        attr = Masses(vals)
        return attr

    def parse_elements(self, atoms_per, numlines):
        vals = self.parsesection_mapper(
            atoms_per, numlines, lambda x: NUMBER_TO_ELEMENT[int(x)])
        attr = Elements(np.array(vals, dtype=object))
        return attr

    def parse_types(self, atoms_per, numlines):
        vals = self.parsesection_mapper(
            atoms_per, numlines, lambda x: x)
        attr = Atomtypes(np.array(vals, dtype=object))
        return attr

    def parse_type_indices(self, atoms_per, numlines):
        vals = self.parsesection_mapper(
            atoms_per, numlines, lambda x: int(x))
        attr = TypeIndices(np.array(vals, dtype=np.int32))
        return attr

    def parse_residx(self, atoms_per, numlines):
        vals = self.parsesection_mapper(
            atoms_per, numlines, lambda x: int(x) - 1)
        return vals

    def parsebond(self, atoms_per, numlines):
        y = self.topfile.next().strip("%FORMAT(")
        section = []
        for i in range(numlines):
            l = self.topfile.next()
            # Subtract 1 from each number to ensure zero-indexing for the atoms
            fields = np.int64(l.split()) - 1
            for j in range(0, len(fields), atoms_per):
                section.append(tuple(fields[j:j+atoms_per]))
        return section

    def parsesection_mapper(self, atoms_per, numlines, mapper):
        section = []
        y = self.topfile.next().strip("%FORMAT(")
        y.strip(")")
        x = FORTRANReader(y)
        for i in range(numlines):
            l = self.topfile.next()
            for j in range(len(x.entries)):
                val = l[x.entries[j].start:x.entries[j].stop].strip()
                if val:
                    section.append(mapper(val))
        return section
