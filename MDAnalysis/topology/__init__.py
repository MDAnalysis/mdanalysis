# $Id$
"""
:mod:`MDAnalysis.topology` -- topology readers
==============================================

This submodule contains topology readers; at the moment, only
Charmm/XPLOR psf is supported well.

A rudimentary PDB parser obtains the list of atoms from a PDB file but
it lacks bond information, charges, and masses at the moment.
"""

__all__ = ['core', 'PSFParser', 'PDBParser', 'GROParser']

import core
import PSFParser, PDBParser, GROParser

# dictionary of known file formats and the corresponding file parser
# (all parser should essentially do the same thing; the PSFParser is
# the reference implementation). The keys in :data:`_topology_parsers`
# are the known topology formats.
_topology_parsers = {'psf': PSFParser.parse,
                     'pdb': PDBParser.parse,
                     'gro': GROParser.parse,
                     }
