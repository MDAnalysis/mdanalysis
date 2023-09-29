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

"""
AMBER PRMTOP topology parser
============================

Reads an AMBER top file to build the system.

Amber keywords are turned into the following attributes:

+----------------------------+----------------------+
| AMBER flag                 | MDAnalysis attribute |
+============================+======================+
| ATOM_NAME                  | names                |
+----------------------------+----------------------+
| CHARGE                     | charges              |
+----------------------------+----------------------+
| ATOMIC_NUMBER              | elements             |
+----------------------------+----------------------+
| MASS                       | masses               |
+----------------------------+----------------------+
| BONDS_INC_HYDROGEN         | bonds                |
| BONDS_WITHOUT_HYDROGEN     |                      |
+----------------------------+----------------------+
| ANGLES_INC_HYDROGEN        | angles               |
| ANGLES_WITHOUT_HYDROGEN    |                      |
+----------------------------+----------------------+
| DIHEDRALS_INC_HYDROGEN     | dihedrals / improper |
| DIHEDRALS_WITHOUT_HYDROGEN |                      |
+----------------------------+----------------------+
| ATOM_TYPE_INDEX            | type_indices         |
+----------------------------+----------------------+
| AMBER_ATOM_TYPE            | types                |
+----------------------------+----------------------+
| RESIDUE_LABEL              | resnames             |
+----------------------------+----------------------+
| RESIDUE_POINTER            | residues             |
+----------------------------+----------------------+

TODO:
  Add support for Chamber-style topologies
  More stringent tests

.. Note::

   The Amber charge is converted to electron charges as used in
   MDAnalysis and other packages. To get back Amber charges, multiply
   by 18.2223.

   Chamber-style Amber topologies (i.e. topologies generated via parmed
   conversion of a CHARMM topology to an AMBER one) are not currently
   supported. Support will likely be added in future MDAnalysis releases.

   As of version 2.0.0, elements are no longer guessed if ATOMIC_NUMBER records
   are missing. In those scenarios, if elements are necessary, users will have
   to invoke the element guessers after parsing the topology file. Please see
   :mod:`MDAnalysis.topology.guessers` for more details.

.. _`PARM parameter/topology file specification`:
   https://ambermd.org/FileFormats.php#topo.cntrl

Classes
-------

.. autoclass:: TOPParser
   :members:
   :inherited-members:

"""
import numpy as np
import itertools

from .tables import Z2SYMB
from ..lib.util import openany, FORTRANReader
from .base import TopologyReaderBase, change_squash
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
    ChainIDs,
    AtomAttr,
    Bonds,
    Angles,
    Dihedrals,
    Impropers
)

import warnings
import logging

logger = logging.getLogger('MDAnalysis.topology.TOPParser')


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
    - Bonds
    - Angles
    - Dihedrals (inc. impropers)
    - ChainIDs (from %RESIDUE_CHAINID)
    - Segids (from %RESIDUE_CHAINID)

    The format is defined in `PARM parameter/topology file
    specification`_.  The reader tries to detect if it is a newer
    (AMBER 12?) file format by looking for the flag "ATOMIC_NUMBER".

    .. _`PARM parameter/topology file specification`:
       https://ambermd.org/FileFormats.php#topo.cntrl

    Additionally, the RESIDUE_CHAINID non-standard flag is supported. This
    can be added with the `addPDB`_ command from parmed:

    .. _`addPDB`: https://parmed.github.io/ParmEd/html/parmed.html#addpdb

    Notes
    -----
    Elements are obtained from the atomic numbers (if present). If a given
    input atomic number does not belong to an element (usually either -1 or 0),
    the element will be assigned an empty record.

    .. versionchanged:: 0.7.6
      parses both amber10 and amber12 formats
    .. versionchanged:: 0.19.0
      parses bonds, angles, dihedrals, and impropers
    .. versionchanged:: 1.0.0
      warns users that chamber-style topologies are not currently supported
    .. versionchanged:: 2.0.0
      no longer guesses elements if missing
    .. versionchanged:: 2.7.0
      gets Segments and chainIDs from flag RESIDUE_CHAINID, when present
    """
    format = ['TOP', 'PRMTOP', 'PARM7']

    def parse(self, **kwargs):
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
            "ATOM_TYPE_INDEX": (1, 10, self.parse_type_indices, "type_indices",
                                0),
            "AMBER_ATOM_TYPE": (1, 20, self.parse_types, "types", 0),
            "RESIDUE_LABEL": (1, 20, self.parse_resnames, "resname", 11),
            "RESIDUE_POINTER": (1, 10, self.parse_residx, "respoint", 11),
            "BONDS_INC_HYDROGEN": (3, 10, self.parse_bonded, "bondh", 2),
            "BONDS_WITHOUT_HYDROGEN": (3, 10, self.parse_bonded, "bonda", 3),
            "ANGLES_INC_HYDROGEN": (4, 10, self.parse_bonded, "angh", 4),
            "ANGLES_WITHOUT_HYDROGEN": (4, 10, self.parse_bonded, "anga", 5),
            "DIHEDRALS_INC_HYDROGEN": (5, 10, self.parse_bonded, "dihh", 6),
            "DIHEDRALS_WITHOUT_HYDROGEN": (5, 10, self.parse_bonded, "diha", 7),
            "RESIDUE_CHAINID": (1, 20, self.parse_chainids, "segids", 11),
        }

        attrs = {}  # empty dict for attrs that we'll fill

        # Open and check top validity
        # Reading header info POINTERS
        with openany(self.filename) as self.topfile:
            header = next(self.topfile)
            if not header.startswith("%VE"):
                raise ValueError(
                    "{0} is not a valid TOP file. %VE Missing in header"
                    "".format(self.filename))
            title = next(self.topfile).split()

            if not (title[1] == "TITLE"):
                # Raise a separate warning if Chamber-style TOP is detected
                if title[1] == "CTITLE":
                    emsg = ("{0} is detected as a Chamber-style TOP file. "
                            "At this time MDAnalysis does not support such "
                            "topologies".format(self.filename))
                else:
                    emsg = ("{0} is not a valid TOP file. "
                            "'TITLE' missing in header".format(self.filename))
                raise ValueError(emsg)

            while not header.startswith('%FLAG POINTERS'):
                header = next(self.topfile)
            next(self.topfile)

            topremarks = [next(self.topfile).strip() for i in range(4)]
            sys_info = [int(k) for i in topremarks for k in i.split()]

            header = next(self.topfile)
            # grab the next section title
            next_section = header.split("%FLAG")[1].strip()

            while next_section is not None:
                try:
                    (num_per_record, per_line,
                     func, name, sect_num) = sections[next_section]
                except KeyError:
                    def next_getter():
                        return self.skipper()
                else:
                    num = sys_info[sect_num] * num_per_record
                    numlines = (num // per_line)
                    if num % per_line != 0:
                        numlines += 1

                    attrs[name] = func(num_per_record, numlines)

                    def next_getter():
                        return next(self.topfile)

                try:
                    line = next_getter()
                    # Capture case where section is empty w/ 1 empty line
                    if numlines == 0 and not line.strip():
                        line = next_getter()
                except StopIteration:
                    next_section = None
                else:
                    try:
                        next_section = line.split("%FLAG")[1].strip()
                    except IndexError:
                        errmsg = (f"%FLAG section not found, formatting error "
                                  f"for PARM7 file {self.filename} ")
                        raise IndexError(errmsg) from None

        # strip out a few values to play with them
        n_atoms = len(attrs['name'])

        resptrs = attrs.pop('respoint')
        resptrs.append(n_atoms)
        residx = np.zeros(n_atoms, dtype=np.int32)
        for i, (x, y) in enumerate(zip(resptrs[:-1], resptrs[1:])):
            residx[x:y] = i

        n_res = len(attrs['resname'])

        # Deal with recreating bonds and angle records here
        attrs['bonds'] = Bonds([i for i in itertools.chain(
                                attrs.pop('bonda'), attrs.pop('bondh'))])

        attrs['angles'] = Angles([i for i in itertools.chain(
                                attrs.pop('anga'), attrs.pop('angh'))])

        attrs['dihedrals'], attrs['impropers'] = self.parse_dihedrals(
                              attrs.pop('diha'), attrs.pop('dihh'))

        # Warn user if elements not in topology
        if 'elements' not in attrs:
            msg = ("ATOMIC_NUMBER record not found, elements attribute will "
                   "not be populated. If needed these can be guessed using "
                   "MDAnalysis.topology.guessers.")
            logger.warning(msg)
            warnings.warn(msg)
        elif np.any(attrs['elements'].values == ""):
            # only send out one warning that some elements are unknown
            msg = ("Unknown ATOMIC_NUMBER value found for some atoms, these "
                   "have been given an empty element record. If needed these "
                   "can be guessed using MDAnalysis.topology.guessers.")
            logger.warning(msg)
            warnings.warn(msg)

        # atom ids are mandatory
        attrs['atomids'] = Atomids(np.arange(n_atoms) + 1)
        attrs['resids'] = Resids(np.arange(n_res) + 1)
        attrs['resnums'] = Resnums(np.arange(n_res) + 1)

        # Amber's 'RESIDUE_CHAINID' is a by-residue attribute, turn it into
        # a by-atom attribute when present. See PR #4007.
        if "segids" in attrs and len(attrs["segids"]) == n_res:
            segidx, (segids,) = change_squash((attrs["segids"],), (attrs["segids"],))
            chainids = [attrs["segids"][r] for r in residx]

            attrs["segids"] = Segids(segids)
            attrs["ChainIDs"] = ChainIDs(chainids)
            n_segs = len(segids)
        else:
            if "segids" in attrs:
                msg = (
                    f"Number of residues ({n_res}) does not match number of "
                    f"%RESIDUE_CHAINID ({len(attrs['segids'])}). Skipping section."
                )
                logger.warning(msg)
                warnings.warn(msg)
            attrs["segids"] = Segids(np.array(["SYSTEM"], dtype=object))
            segidx = None
            n_segs = 1

        top = Topology(
            n_atoms,
            n_res,
            n_segs,
            attrs=list(attrs.values()),
            atom_resindex=residx,
            residue_segindex=segidx,
        )

        return top

    def skipper(self):
        """TOPParser :class: helper function, skips lines of input parm7 file
        until we find the next %FLAG entry and return that

        Returns
        -------
        line : string
            String containing the current line of the parm7 file
        """
        line = next(self.topfile)
        while not line.startswith("%FLAG"):
            line = next(self.topfile)
        return line

    def parse_names(self, num_per_record, numlines):
        """Extracts atoms names from parm7 file

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in the section (unused input)
        numlines : int
            The number of lines to be parsed in current section

        Returns
        -------
        attr : :class:`Atomnames`
            A :class:`Atomnames` instance containing the names of each atom as
            defined in the parm7 file
        """
        vals = self.parsesection_mapper(numlines, lambda x: x)
        attr = Atomnames(np.array(vals, dtype=object))
        return attr

    def parse_resnames(self, num_per_record, numlines):
        """Extracts the names of each residue

        Parameters
        ----------
        num_per_record : int
            The number of entries for each recrod in section (unused input)
        numlines : int
            The number of lines to be parsed in current section

        Returns
        -------
        attr : :class:`Resnames`
            A :class:`Resnames` instance containing the names of each residue
            as defined in the parm7 file
        """
        vals = self.parsesection_mapper(numlines, lambda x: x)
        attr = Resnames(np.array(vals, dtype=object))
        return attr

    def parse_charges(self, num_per_record, numlines):
        """Extracts the partial charges for each atom

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in section (unused input)
        numlines : int
            The number of lines to be parsed in current section

        Returns
        -------
        attr : :class:`Charges`
            A :class:`Charges` instance containing the partial charges of each
            atom as defined in the parm7 file
        """
        vals = self.parsesection_mapper(numlines, lambda x: float(x))
        charges = np.array(vals, dtype=np.float32)
        charges /= 18.2223  # to electron charge units
        attr = Charges(charges)
        return attr

    def parse_masses(self, num_per_record, numlines):
        """Extracts the mass of each atom

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in section (unused input)
        numlines : int
            The number of lines to be parsed in current section

        Returns
        -------
        attr : :class:`Masses`
            A :class:`Masses` instance containing the mass of each atom as
            defined in the parm7 file
        """
        vals = self.parsesection_mapper(numlines, lambda x: float(x))
        attr = Masses(vals)
        return attr

    def parse_elements(self, num_per_record, numlines):
        """Extracts the atomic numbers of each atom and converts to element type

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in section(unused input)
        numlines : int
            The number of lines to be pasred in current section

        Returns
        -------
        attr : :class:`Elements`
            A :class:`Elements` instance containing the element of each atom
            as defined in the parm7 file

        Note
        ----
        If the record contains unknown atomic numbers (e.g. <= 0), these will
        be treated as unknown elements and assigned an empty string value. See
        issues #2306 and #2651 for more details.

        .. versionchanged:: 2.0.0
           Unrecognised elements will now return a empty string. The parser
           will no longer attempt to guess the element by default.
        """

        vals = self.parsesection_mapper(
                numlines,
                lambda x: Z2SYMB[int(x)] if int(x) > 0 else "")
        attr = Elements(np.array(vals, dtype=object))
        return attr

    def parse_types(self, num_per_record, numlines):
        """Extracts the force field atom types of each atom

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in section (unused input)
        numlines : int
            The number of lines to be parsed in current section

        Returns
        -------
        attr : :class:`Atomtypes`
            A :class:`Atomtypes` instance containing the atom types for each
            atom as defined in the parm7 file
        """
        vals = self.parsesection_mapper(numlines, lambda x: x)
        attr = Atomtypes(np.array(vals, dtype=object))
        return attr

    def parse_type_indices(self, num_per_record, numlines):
        """Extracts the index of atom types of the each atom involved in Lennard
        Jones (6-12) interactions.

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in section (unused input)
        numlines : int
            The number of lines to be parsed in current section

        Returns
        -------
        attr :class:`TypeIndices`
            A :class:`TypeIndices` instance containing the LJ 6-12 atom type
            index for each atom
        """
        vals = self.parsesection_mapper(numlines, lambda x: int(x))
        attr = TypeIndices(np.array(vals, dtype=np.int32))
        return attr

    def parse_residx(self, num_per_record, numlines):
        """Extracts the residue pointers for each atom

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in section (unused input)
        numlines : int
            The number of lines to be parsed in current section

        Returns
        -------
        vals : list of int
            A list of zero-formatted residue pointers for each atom
        """
        vals = self.parsesection_mapper(numlines, lambda x: int(x) - 1)
        return vals

    def parse_chunks(self, data, chunksize):
        """Helper function to parse AMBER PRMTOP bonds/angles.

        Parameters
        ----------
        data : list of int
            Input list of the parm7 bond/angle section, zero-indexed
        num_per_record : int
            The number of entries for each record in the input list

        Returns
        -------
        vals : list of int tuples
            A list of tuples containing the atoms involved in a given bonded
            interaction

        Note
        ----
        In the parm7 format this information is structured in the following
        format: [ atoms 1:n, internal index ]
        Where 1:n represent the ids of the n atoms involved in the bond/angle
        and the internal index links to a given set of FF parameters.
        Therefore, to extract the required information, we split out the list
        into chunks of size num_per_record, and only extract the atom ids.
        """
        vals = [tuple(data[x:x+chunksize-1])
                for x in range(0, len(data), chunksize)]
        return vals

    def parse_bonded(self, num_per_record, numlines):
        """Extracts bond information from PARM7 format files

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in section
        numlines : int
            The number of lines to be parsed for this section

        Note
        ----
        For the bond/angle sections of parm7 files, the atom numbers are set to
        coordinate array index values. As detailed in `the specification`_,
        to recover the actual atom number, one
        should divide the values by 3 and add 1. Here, since we want to satisfy
        zero-indexing, we only divide by 3.

        .. _`the specification`: https://ambermd.org/FileFormats.php#topo.cntrl
        """
        fields = self.parsesection_mapper(numlines, lambda x: int(x) // 3)
        section = self.parse_chunks(fields, num_per_record)
        return section

    def parsesection_mapper(self, numlines, mapper):
        """Parses FORTRAN formatted section, and returns a list of all entries
        in each line

        Parameters
        ----------
        numlines : int
            The number of lines to be parsed in this section
        mapper : lambda operator
            Operator to format entries in current section

        Returns
        -------
        section : list
            A list of all entries in a given parm7 section
        """
        section = []

        def get_fmt(file):
            """ Skips '%COMMENT' lines until it gets the FORMAT specification
            for the section."""
            line = next(file)
            if line[:7] == "%FORMAT":
                return line[8:].split(")")[0]
            elif line[:8] == "%COMMENT":
                return get_fmt(file)
            else:
                raise ValueError(
                    "Invalid header line. Does not begin with either %FLAG, %FORMAT "
                    f"nor %COMMENT:\n{line}"
                )

        # There may be %COMMENT lines between %FLAG and %FORMAT statements. Skip them.
        fmt = get_fmt(self.topfile)
        x = FORTRANReader(fmt)
        for i in range(numlines):
            l = next(self.topfile)
            for j in range(len(x.entries)):
                val = l[x.entries[j].start:x.entries[j].stop].strip()
                if val:
                    section.append(mapper(val))
        return section

    def parse_dihedrals(self, diha, dihh):
        """Combines hydrogen and non-hydrogen containing AMBER dihedral lists
        and extracts sublists for conventional dihedrals and improper angles

        Parameters
        ----------
        diha : list of tuples
            The atom ids of dihedrals not involving hydrogens
        dihh : list of tuples
            The atom ids of dihedrals involving hydrogens

        Returns
        -------
        dihedrals : :class:`Dihedrals`
            A :class:`Dihedrals` instance containing a list of all unique
            dihedrals as defined by the parm7 file
        impropers : :class:`Impropers`
            A :class:`Impropers` instance containing a list of all unique
            improper dihedrals as defined by the parm7 file

        Note
        ----
        As detailed in `the specification`_, the dihedral sections
        of parm7 files contain information about both conventional dihedrals
        and impropers. The following must be accounted for:
        1) If the fourth atom in a dihedral entry is given a negative value,
        this indicates that it is an improper.
        2) If the third atom in a dihedral entry is given a negative value,
        this indicates that it 1-4 NB interactions are ignored for this
        dihedrals. This could be due to the dihedral within a ring, or if it is
        part of a multi-term dihedral definition or if it is an improper.

        .. _`the specification`: https://ambermd.org/FileFormats.php#topo.cntrl
        """
        improp = []
        dihed = []
        for i in itertools.chain(diha, dihh):
            if i[3] < 0:
                improp.append(i[:2]+(abs(i[2]),)+(abs(i[3]),))
            elif i[2] < 0:
                vals = i[:2] + (abs(i[2]),) + i[3:]
                dihed.append(vals)
            else:
                dihed.append(i)
        dihed = sorted(set(dihed))
        dihedrals = Dihedrals(dihed)
        impropers = Impropers(improp)
        return dihedrals, impropers

    def parse_chainids(self, num_per_record: int, numlines: int):
        """Extracts the chainID of each residue

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in section (unused input)
        numlines : int
            The number of lines to be parsed in current section

        Returns
        -------
        attr : numpy array
            A numpy array containing the chainID of each residue as defined in
            the parm7 file
        """
        vals = self.parsesection_mapper(numlines, lambda x: x)
        attr = np.array(vals)
        return attr
