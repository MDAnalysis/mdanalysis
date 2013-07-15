# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
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
   CHARMM/XPLOR     psf         reads either format, but only atoms and bonds
                                information is used at the moment;
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
   ================ ==========  =====================================================

.. [#a] This format can also be used to provide *coordinates* so that
   it is possible to create a full
   :mod:`~MDAnalysis.core.AtomGroup.Universe` by simply providing a
   file of this format as the sole argument to
   :mod:`~MDAnalysis.core.AtomGroup.Universe`: ``u =
   Universe(filename)``

.. SeeAlso:: :ref:`Coordinates` with the :ref:`Supported coordinate formats`
"""

__all__ = ['core', 'PSFParser', 'PDBParser', 'PQRParser', 'GROParser', 'CRDParser','TOPParser', 'PDBQTParser']

import core
import PSFParser, TOPParser, \
    PDBParser, PrimitivePDBParser, PQRParser, GROParser, CRDParser, \
    PDBQTParser, DMSParser

# dictionary of known file formats and the corresponding file parser
# (all parser should essentially do the same thing; the PSFParser is
# the reference implementation). The keys in :data:`_topology_parsers`
# are the known topology formats.
_topology_parsers = {'PSF': PSFParser.parse,
                     'PDB': PDBParser.parse,
                     'PQR': PQRParser.parse,
                     'GRO': GROParser.parse,
                     'CRD': CRDParser.parse,
                     'TOP': TOPParser.parse,
                     'PRMTOP': TOPParser.parse,
                     'PDBQT': PDBQTParser.parse,
                     'DMS': DMSParser.parse,
                     }
_topology_parsers_permissive = _topology_parsers.copy()
_topology_parsers_permissive['PDB'] = PrimitivePDBParser.parse
