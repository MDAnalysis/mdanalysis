# $Id$
"""
:mod:`MDAnalysis.topology` -- topology readers
==============================================

This submodule contains topology readers.

A rudimentary PDB parser obtains the list of atoms from a PDB file but
it lacks bond information, charges, and masses at the moment.

Supporte topologies; formats marked with ans asterisk * also hold coordinates.

=============  ==========  =====================================================
Name           extension   remarks
=============  ==========  =====================================================
CHARMM/XPLOR   psf         reads either format, but only atoms and bonds 
                           information is used at the moment

CHARMM*        crd         "CARD" coordinate output from CHARMM; deals with
                           either standard or EXTended format

Brookhaven*    pdb         a simplified PDB format (as used in MD simulations)
                           is read by default; the full format can be read by
                           supplying the `permissive=False` flag to 
                           :class:`MDAnalysis.Universe`

PQR            pqr         PDB-like but whitespace-separated files with charge
                           and radius information

GROMOS96*      gro         

Amber          top         simple Amber format read (only supports a subset of 
               prmtop      flags)
=============  ==========  =====================================================
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
