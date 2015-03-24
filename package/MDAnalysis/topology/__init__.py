# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
Topology readers --- :mod:`MDAnalysis.topology`
===============================================

This submodule contains the topology readers. A topology file supplies the list
of atoms in the system, their connectivity and possibly additional information
such as B-factors, partial charges, etc. The details depend on the file format
and not every topology file provides all (or even any) additional data. As a
minimum, a topology file has to contain the names of atoms in the order of the
coordinate file and their residue names and numbers.

The following table lists the currently supported topology formats.

.. _`Supported topology formats`:

.. table:: Table of Supported topology formats

   ================ ==========  =====================================================
   Name             extension   remarks
   ================ ==========  =====================================================
   CHARMM/XPLOR     psf         reads either format, atoms, bonds, angles,
                                torsions/dihedrals information is all used;
                                :mod:`MDAnalysis.topology.PSFParser`

   CHARMM [#a]_     crd         "CARD" coordinate output from CHARMM; deals with
                                either standard or EXTended format;
                                :mod:`MDAnalysis.topology.CRDParser`

   Brookhaven [#a]_ pdb         a simplified PDB format (as used in MD simulations)
                                is read by default; the full format can be read by
                                supplying the `permissive=False` flag to
                                :class:`MDAnalysis.Universe`;
                                :mod:`MDAnalysis.topology.PrimitivePDBParser` and
                                :mod:`MDAnalysis.topology.PDBParser`

   XPDB             pdb         Extended PDB format (can use 5-digit residue
                                numbers). To use, specify the format "XPBD"
                                explicitly:
                                ``Universe(..., topology_format="XPDB")``.
                                Module :MDAnalysis.coordinates.PDB`

   PQR [#a]_        pqr         PDB-like but whitespace-separated files with charge
                                and radius information;
                                :mod:`MDAnalysis.topology.PQRParser`

   PDBQT [#a]_      pdbqt       file format used by AutoDock with atom types *t*
                                and partial charges *q*. Module:
                                :mod:`MDAnalysis.topology.PDBQTParser`

   GROMOS96 [#a]_   gro         GROMOS96 coordinate file;
                                :mod:`MDAnalysis.topology.GROParser`

   AMBER            top,        simple AMBER format reader (only supports a subset
                    prmtop      of flags);
                                :mod:`MDAnalysis.topology.TOPParser`

   DESRES [#a]_     dms         DESRES molecular sturcture reader (only supports
                                the atom and bond records);
                                :mod:`MDAnalysis.topology.DMSParser`

   TPR              tpr         Gromacs portable run input reader (limited
                                experimental support for some of the more recent
                                versions of the file format);
                                :mod:`MDAnalysis.topology.TPRParser`
   MOL2             mol2        Tripos MOL2 molecular structure format;
                                :mod:`MDAnalysis.topology.MOL2Parser`
   LAMMPS           data        LAMMPS Data file parser
                                :mod:`MDAnalysis.topology.LAMMPSParser`
   XYZ [#a]_        xyz         XYZ File Parser.  Reads only the labels from atoms and
                                constructs minimal topology data.
                                :mod:`MDAnalysis.topology.XYZParser`
   GMS              gms,        GAMESS output parser. Read only atoms of assembly 
                    log         section (atom, elems and coords) and construct topology.
                                :mod:`MDAnalysis.topology.GMSParser`
   ================ ==========  =====================================================

.. [#a] This format can also be used to provide *coordinates* so that
   it is possible to create a full
   :mod:`~MDAnalysis.core.AtomGroup.Universe` by simply providing a
   file of this format as the sole argument to
   :mod:`~MDAnalysis.core.AtomGroup.Universe`: ``u =
   Universe(filename)``

.. SeeAlso:: :ref:`Coordinates` with the :ref:`Supported coordinate formats`


Developer Notes
---------------

.. versionadded:: 0.8

Topology information consists of data that do not change over time,
i.e. information that is the same for all time steps of a
trajectory. This includes

* identity of atoms (name, type, number, partial charge, ...) and to
  which residue and segment they belong; atoms are identified in
  MDAnalysis by their :attr:`~MDAnalysis.core.AtomGroup.Atom.number`,
  an integer number starting at 0 and incremented in the order of
  atoms found in the topology.

* bonds (pairs of atoms)

* angles (triplets of atoms)

* dihedral angles (quadruplets of atoms) — proper and improper
  dihedrals should be treated separately

At the moment, only the identity of atoms is mandatory and at its most
basic, the topology is simply a list of atoms to be associated with a
list of coordinates.

The current implementation contains submodules for different topology
file types. Each submodule *must* contain a function :func:`parse`:

.. function: parse(filename)

   Read a topology from *filename* and return the structure dict.

The function returns the basic MDAnalysis representation of the
topology. At the moment, this is simply a dictionary with keys
*_atoms*, *_bonds*, *_angles*, *_dihe*, *_impr*. The dictionary is
stored as :attr:`MDAnalysis.AtomGroup.Universe._psf`.

.. warning::

   The internal dictionary representation is subject to change. User
   code should *not* access this dictionary directly. The information
   provided here is solely for developers who need to work with the
   existing parsers.

.. SeeAlso:: `Topology Data Structures Wiki page`_

.. _`Topology Data Structures Wiki page`:
   http://code.google.com/p/mdanalysis/wiki/TopologyDataStructures

The format of the individual keys is the following (see
:mod:`PSFParser` for a reference implementation):

_atoms
~~~~~~

The **atoms** are represented as a :class:`list` of
:class:`~MDAnalysis.core.AtomGroup.Atom` instances. The parser needs
to initialize the :class:`~MDAnalysis.core.AtomGroup.Atom` objects
with the data read from the topology file.

The order of atoms in the list must correspond to the sequence of
atoms in the topology file. The atom's
:attr:`~MDAnalysis.core.AtomGroup.Atom.number` corresponds to its
index in this list.


_bonds
~~~~~~

**Bonds** are represented as a :class:`tuple` of :class:`tuple`. Each tuple
contains two atom numbers, which indicate the atoms between which the
bond is formed. Only one of the two permutations is stored, typically
the one with the lower atom number first.


_bondorder
~~~~~~~~~~

Some **bonds** have additional information called **order**. When available
this is stored in a dictionary of format {bondtuple:order}.  This extra
information is then passed to Bond initialisation in u._init_bonds()


_angles
~~~~~~~

**Angles** are represented by a :class:`list` of :class:`tuple`. Each
tuple contains three atom numbers.

.. Note::

   At the moment, the order is not defined and depends on how the
   topology file defines angles.


_dihe
~~~~~

**Proper dihedral angles** are represented by a :class:`list` of :class:`tuple`. Each
tuple contains four atom numbers.

.. Note::

   At the moment, the order is not defined and depends on how the
   topology file defines proper dihedrals..



_impr
~~~~~

**Improper dihedral angles** are represented by a :class:`list` of :class:`tuple`. Each
tuple contains four atom numbers.

.. Note::

   At the moment, the order is not defined and depends on how the
   topology file defines improper dihedrals..

"""

__all__ = ['core', 'PSFParser', 'PDBParser', 'PQRParser', 'GROParser',
           'CRDParser', 'TOPParser', 'PDBQTParser', 'TPRParser',
           'LAMMPSParser', 'XYZParser', 'GMSParser']

import core
import PSFParser
import TOPParser
import PDBParser
import PrimitivePDBParser
import ExtendedPDBParser
import PQRParser
import GROParser
import CRDParser
import PDBQTParser
import DMSParser
import TPRParser
import MOL2Parser
import LAMMPSParser
import XYZParser
import GMSParser


# dictionary of known file formats and the corresponding file parser
# (all parser should essentially do the same thing; the PSFParser is
# the reference implementation). The keys in :data:`_topology_parsers`
# are the known topology formats.
_topology_parsers = {'PSF': PSFParser.PSFParser,
                     'PDB': PDBParser.PDBParser,
                     'Permissive_PDB': PrimitivePDBParser.PrimitivePDBParser,
                     'XPDB': ExtendedPDBParser.ExtendedPDBParser,
                     'PQR': PQRParser.PQRParser,
                     'GRO': GROParser.GROParser,
                     'CRD': CRDParser.CRDParser,
                     'TOP': TOPParser.TOPParser,
                     'PRMTOP': TOPParser.TOPParser,
                     'PDBQT': PDBQTParser.PDBQTParser,
                     'TPR': TPRParser.TPRParser,
                     'DMS': DMSParser.DMSParser,
                     'MOL2': MOL2Parser.MOL2Parser,
                     'DATA': LAMMPSParser.DATAParser,
                     'XYZ': XYZParser.XYZParser,
                     'GMS': GMSParser.GMSParser,
                     }
