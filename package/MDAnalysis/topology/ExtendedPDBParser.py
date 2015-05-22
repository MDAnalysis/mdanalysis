# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
Extended PDB topology parser
============================

.. versionadded:: 0.8

This topology parser uses a PDB file to build a minimum internal structure
representation (list of atoms). The only difference from
:mod:`~MDAnalysis.topology.PrimitivePDBParser` is that this parser reads a
non-standard PDB-like format in which residue numbers can be five digits
instead of four.

The topology reader reads a PDB file line by line and ignores atom numbers but
reads residue numbers up to 99,999 correctly. If you have systems containing at
least 100,000 residues then you need to use a different file format that can
handle such residue numbers.

.. Note::

   The parser processes atoms and their names. Masses are guessed and set to 0
   if unknown. Partial charges are not set.

.. SeeAlso::

   * :mod:`MDAnalysis.topology.PrimitivePDBParser`
   * :class:`MDAnalysis.coordinates.PDB.ExtendedPDBReader`
   * :class:`MDAnalysis.core.AtomGroup.Universe`

Classes
-------

.. autoclass:: ExtendedPDBParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

from . import PrimitivePDBParser


class ExtendedPDBParser(PrimitivePDBParser.PrimitivePDBParser):
    """Parser that obtains a list of atoms from an non-standard "extended" PDB file.

    Extended PDB files (MDAnalysis format specifier *XPDB*) may contain residue
    sequence numbers up to 99,999 by utilizing the insertion character field of
    the PDB standard.

    .. SeeAlso:: :class:`MDAnalysis.coordinates.PDB.ExtendedPDBReader`

    .. versionadded:: 0.8
    """
    format = 'XPDB'
