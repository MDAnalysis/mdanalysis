# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
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
and not every topology file provides all (or even any) additional data. This
data is made accessible through AtomGroup properties.

As a minimum, all topology parsers will provide atom ids, atom types, masses,
resids, resnums and segids as well as assigning all atoms to residues and all
residues to segments.  For systems without residues and segments, this results
in there being a single residue and segment to which all atoms belong. Often
when data is not provided by a file, it will be guessed based on other data in
the file.  In the event that this happens, a UserWarning will always be issued.

The following table lists the currently supported topology formats along with
the attributes they provide.

.. _`Supported topology formats`:

.. table:: Table of Supported Topology Formats

   ================= ========= ================= ===================================================
   Name              extension attributes        remarks
   ================= ========= ================= ===================================================
   CHARMM/XPLOR PSF  psf       resnames,         :mod:`MDAnalysis.topology.PSFParser`
                               names, types,
                               charges,
                               bonds, angles,
                               dihedrals,
                               impropers

   CHARMM CARD [#a]_ crd       names,            "CARD" coordinate output from CHARMM; deals with
                               tempfactors,      either standard or EXTended format;
                               resnames,         :mod:`MDAnalysis.topology.CRDParser`

   Brookhaven [#a]_  pdb/ent   names, bonds,     a simplified PDB format (as used in MD simulations)
                               resids, resnums,  is read by default
                               types,
                               chainids,
                               occupancies,
                               bfactors,
                               resids, icodes,
                               resnames,
                               segids,

   XPDB [#a]_        pdb       As PDB except     Extended PDB format (can use 5-digit residue
                               icodes            numbers). To use, specify the format "XPBD"
                                                  explicitly:
                                                  ``Universe(..., topology_format="XPDB")``.
                                                  Module :mod:`MDAnalysis.coordinates.PDB`

   PQR [#a]_         pqr       names, charges,   PDB-like but whitespace-separated files with charge
                               types,            and radius information;
                               radii, resids,    :mod:`MDAnalysis.topology.PQRParser`
                               resnames, segids

   PDBQT [#a]_       pdbqt     names, types,     file format used by AutoDock with atom types and
                               altLocs, charges, partial charges. Module:
                               resnames,         :mod:`MDAnalysis.topology.PDBQTParser`
                               resids,
                               icodes,
                               occupancies,
                               tempfactors,
                               segids,

   GROMOS96 [#a]_    gro       names, resids,    GROMOS96 coordinate file;
                               resnames,         :mod:`MDAnalysis.topology.GROParser`

   AMBER             top,      names, charges    simple AMBER format reader (only supports a subset
                     prmtop,   type_indices,     of flags);
                     parm7     types,            :mod:`MDAnalysis.topology.TOPParser`
                               resnames,

   DESRES [#a]_      dms       names, numbers,   DESRES molecular sturcture reader (only supports
                               masses, charges,  the atom and bond records);
                               chainids, resids, :mod:`MDAnalysis.topology.DMSParser`
                               resnames, segids,
                               radii,   

   TPR [#b]_         tpr       names, types,     Gromacs portable run input reader (limited
                               resids, resnames, experimental support for some of the more recent
                               charges, bonds,   versions of the file format);
                               masses,           :mod:`MDAnalysis.topology.TPRParser`

   MOL2 [#a]_        mol2      ids, names,       Tripos MOL2 molecular structure format;
                               types, resids,    :mod:`MDAnalysis.topology.MOL2Parser`
                               charges, bonds,
                               resnames,

   LAMMPS [#a]_      data      ids, types,       LAMMPS Data file parser
                               masses, charges,  :mod:`MDAnalysis.topology.LAMMPSParser`
                               resids, bonds,
                               angles, dihedrals

   XYZ [#a]_         xyz       names             XYZ File Parser.  Reads only the labels from atoms
                                                 and constructs minimal topology data.
                                                 :mod:`MDAnalysis.topology.XYZParser`

   GAMESS [#a]_      gms,      names,            GAMESS output parser. Read only atoms of assembly
                     log       atomic charges,   section (atom, elems and coords) and construct
                                                 topology.
                                                 :mod:`MDAnalysis.topology.GMSParser`

   DL_Poly [#a]_     config,   ids, names        DL_Poly CONFIG or HISTORY file.  Reads only the
                     history                     atom names. If atoms are written out of order, will
                                                 correct the order.
                                                 :mod:`MDAnalysis.topology.DLPolyParser`

   Hoomd XML         xml       types, charges,   `HOOMD XML`_ topology file.  Reads atom types,
                               radii, masses     masses, and charges if possible. Also reads bonds,
                               bonds, angles,    angles, and dihedrals.
                               dihedrals         :mod:`MDAnalysis.topology.HoomdXMLParser`

   Macromolecular    mmtf      altLocs,          `Macromolecular Transmission Format (MMTF)`_.
   transmission                bfactors, bonds,  An efficient compact format for biomolecular 
   format                      charges, masses,  structures.
                               names,
                               occupancies,
                               types, icodes,
                               resnames, resids,
                               segids, models
   ================= ========= ================= ===================================================

.. [#a] This format can also be used to provide *coordinates* so that
   it is possible to create a full
   :mod:`~MDAnalysis.core.AtomGroup.Universe` by simply providing a
   file of this format as the sole argument to
   :mod:`~MDAnalysis.core.AtomGroup.Universe`: ``u =
   Universe(filename)``

.. [#b] The Gromacs TPR format contains coordinate information but
        parsing coordinates from a TPR file is currently not implemented
        in :mod:`~MDAnalysis.topology.TPRParser`.

.. SeeAlso:: :ref:`Coordinates` with the :ref:`Supported coordinate formats`

.. _HOOMD XML: http://codeblue.umich.edu/hoomd-blue/doc/page_xml_file_format.html
.. _Macromolecular Transmission Format (MMTF): https://mmtf.rcsb.org/

Developer Notes
---------------

.. versionadded:: 0.8

Topology information consists of data that do not change over time,
i.e. information that is the same for all time steps of a
trajectory. This includes

* identity of atoms (name, type, number, partial charge, ...) and to
  which residue and segment they belong; atoms are identified in
  MDAnalysis by their :attr:`~MDAnalysis.core.AtomGroup.Atom.index`,
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
*atoms*, *bonds*, *angles*, *dihedrals*, *impropers*. The dictionary is
stored as :attr:`MDAnalysis.AtomGroup.Universe._topology`.

.. warning::

   The internal dictionary representation is subject to change. User
   code should *not* access this dictionary directly. The information
   provided here is solely for developers who need to work with the
   existing parsers.

.. SeeAlso:: `Topology Data Structures Wiki page`_

.. _`Topology Data Structures Wiki page`:
   https://github.com/MDAnalysis/mdanalysis/wiki/TopologyDataStructures

The format of the individual keys is the following (see
:mod:`PSFParser` for a reference implementation):

atoms
~~~~~~

The **atoms** are represented as a :class:`list` of
:class:`~MDAnalysis.core.AtomGroup.Atom` instances. The parser needs
to initialize the :class:`~MDAnalysis.core.AtomGroup.Atom` objects
with the data read from the topology file.

The order of atoms in the list must correspond to the sequence of
atoms in the topology file. The atom's
:attr:`~MDAnalysis.core.AtomGroup.Atom.index` corresponds to its
index in this list.


bonds
~~~~~~

**Bonds** are represented as a :class:`tuple` of :class:`tuple`. Each tuple
contains two atom numbers, which indicate the atoms between which the
bond is formed. Only one of the two permutations is stored, typically
the one with the lower atom number first.


bondorder
~~~~~~~~~~

Some **bonds** have additional information called **order**. When available
this is stored in a dictionary of format {bondtuple:order}.  This extra
information is then passed to Bond initialisation in u._init_bonds()


angles
~~~~~~~

**Angles** are represented by a :class:`list` of :class:`tuple`. Each
tuple contains three atom numbers.  The second of these numbers
represents the apex of the angle.


dihedrals
~~~~~~~~~

**Proper dihedral angles** are represented by a :class:`list` of :class:`tuple`. Each
tuple contains four atom numbers. The angle of the torsion
is defined by the angle between the planes formed by atoms 1, 2, and 3,
and 2, 3, and 4.


impropers
~~~~~~~~~

**Improper dihedral angles** are represented by a :class:`list` of :class:`tuple`. Each
tuple contains four atom numbers.  The angle of the improper torsion
is again defined by the angle between the planes formed by atoms 1, 2, and 3,
and 2, 3, and 4.  Improper dihedrals differ from regular dihedrals as the
four atoms need not be sequentially bonded, and are instead often all bonded
to the second atom.

"""

__all__ = ['core', 'PSFParser', 'PDBParser', 'PQRParser', 'GROParser',
           'CRDParser', 'TOPParser', 'PDBQTParser', 'TPRParser',
           'LAMMPSParser', 'XYZParser', 'GMSParser', 'DLPolyParser',
           'HoomdXMLParser']
# Registry of all Parsers in MDAnalysis
# Gets filled on class definition by metaclass magic
_PARSERS = {}

from . import core
from . import PSFParser
from . import TOPParser
from . import PDBParser
from . import PrimitivePDBParser
from . import ExtendedPDBParser
from . import PQRParser
from . import GROParser
from . import CRDParser
from . import PDBQTParser
from . import DMSParser
from . import TPRParser
from . import MOL2Parser
from . import LAMMPSParser
from . import XYZParser
from . import GMSParser
from . import DLPolyParser
from . import HoomdXMLParser
from . import MMTFParser
 
