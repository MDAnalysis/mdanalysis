# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Primitive PDB topology parser
=============================

This topology parser uses a standard PDB file to build a minimum
internal structure representation (list of atoms).

The topology reader reads a PDB file line by line and ignores atom
numbers but only reads residue numbers up to 9,999 correctly. If you
have systems containing at least 10,000 residues then you need to use
a different file format (e.g. the "extended" PDB, *XPDB* format, see
:mod:`~MDAnalysis.topology.ExtendedPDBParser`) that can handle residue
numbers up to 99,999.

.. Note::

   The parser processes atoms and their names. Masses are guessed and set to 0
   if unknown. Partial charges are not set.

.. SeeAlso::

   * :mod:`MDAnalysis.topology.ExtendedPDBParser`
   * :class:`MDAnalysis.coordinates.PDB.PDBReader`
   * :class:`MDAnalysis.core.AtomGroup.Universe`

Classes
-------

.. autoclass:: PrimitivePDBParser
   :members:
   :inherited-members:

..deprecated:: 0.15.0
    PDBParser has been replaced with PrimitivePDBParser.
"""

from __future__ import absolute_import, print_function

import numpy as np
import warnings

from . import PDBParser
from ..core.AtomGroup import Atom
from .core import get_atom_mass, guess_atom_element
from ..lib.util import openany
from .base import TopologyReader



class PrimitivePDBParser(PDBParser.PDBParser):
    def __init__(self, *args, **kwargs):
        warnings.warn('PrimitivePDBParser is identical to the PDBParser,'
                    ' it is deprecated in favor of the shorter name',
                    category=DeprecationWarning)
        super(PDBParser.PDBParser, self).__init__(*args, **kwargs)

def _parse_conect(conect):

    """parse a CONECT record from pdbs

    Parameters
    ----------
    conect : str
        white space striped CONECT record

    Returns
    -------
    atom_id : int
        atom index of bond
    bonds : set
        atom ids of bonded atoms

    Raises
    ------
    RuntimeError
        Raised if ``conect`` is not a valid CONECT record
    """
    atom_id = np.int(conect[6:11])
    n_bond_atoms = len(conect[11:]) // 5
    if len(conect[11:]) % n_bond_atoms != 0:
        raise RuntimeError("Bond atoms aren't aligned proberly for CONECT "
                           "record: {}".format(conect))
    bond_atoms = (int(conect[11 + i * 5: 16 + i * 5]) for i in
                  range(n_bond_atoms))
    return atom_id, bond_atoms
