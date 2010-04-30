# $Id$
"""
:mod:`MDAnalysis.topology` -- topology readers
==============================================

This submodule contains topology readers; at the moment, only
Charmm/XPLOR psf is supported well.

A rudimentary PDB parser obtains the list of atoms from a PDB file but
it lacks bond information, charges, and masses at the moment.
"""

__all__ = ['core', 'PSFParser']

import core, PSFParser
