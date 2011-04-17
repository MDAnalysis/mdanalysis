# $Id$
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

=============  ==========  =====================================================
Name           extension   remarks
=============  ==========  =====================================================
CHARMM/XPLOR   psf         reads either format, but only atoms and bonds
                           information is used at the moment;
                           :mod:`MDAnalysis.topology.PSFParser`

CHARMM*        crd         "CARD" coordinate output from CHARMM; deals with
                           either standard or EXTended format;
                           :mod:`MDAnalysis.topology.CRDParser`

Brookhaven*    pdb         a simplified PDB format (as used in MD simulations)
                           is read by default; the full format can be read by
                           supplying the `permissive=False` flag to
                           :class:`MDAnalysis.Universe`;
                           :mod:`MDAnalysis.topology.PrimitivePDBParser` and
                           :mod:`MDAnalysis.topology.PDBParser`

PQR*           pqr         PDB-like but whitespace-separated files with charge
                           and radius information;
                           :mod:`MDAnalysis.topology.PQRParser`

GROMOS96*      gro         GROMOS96 coordinate file;
                           :mod:`MDAnalysis.topology.GROParser`

Amber          top         simple Amber format read (only supports a subset of
               prmtop      flags);
                           :mod:`MDAnalysis.topology.TOPParser`
=============  ==========  =====================================================

Formats marked with ans asterisk * also hold coordinates and thus can
be used as the sole argument to :class:`MDAnalysis.Universe` to set up
a system.

.. SeeAlso:: :ref:`Coordinates`
"""

__all__ = ['core', 'PSFParser', 'PDBParser', 'PQRParser', 'GROParser', 'CRDParser','TOPParser']

import core
import PSFParser, TOPParser, \
    PDBParser, PrimitivePDBParser, PQRParser, GROParser, CRDParser

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
                     }
_topology_parsers_permissive = _topology_parsers.copy()
_topology_parsers_permissive['PDB'] = PrimitivePDBParser.parse
