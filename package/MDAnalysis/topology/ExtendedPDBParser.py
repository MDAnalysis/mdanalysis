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
   if unknown. Partial charges are not set. Bond connectivity can be guessed if
   the ``bonds=True`` keyword is set for
   :class:`~MDAnalysis.core.AtomGroup.Universe`.

.. SeeAlso::

   * :mod:`MDAnalysis.topology.PrimitivePDBParser`
   * :class:`MDAnalysis.coordinates.PDB.ExtendedPDBReader`
   * :class:`MDAnalysis.core.AtomGroup.Universe`

"""

import MDAnalysis.coordinates.PDB
import PrimitivePDBParser

class ExtendedPDBParser(PrimitivePDBParser.PrimitivePDBParser):
    """Parser that obtains a list of atoms from an non-standard "extended" PDB file.

    Extended PDB files (MDAnalysis format specifier *XPDB*) may contain residue
    sequence numbers up to 99,999 by utilizing the insertion character field of
    the PDB standard.

    .. SeeAlso:: :class:`MDAnalysis.coordinates.PDB.ExtendedPDBReader`

    .. versionadded:: 0.8
    """
    def __init__(self, *args, **kwargs):
        super(ExtendedPDBParser, self).__init__(*args, **kwargs)
        self.PDBReader = MDAnalysis.coordinates.PDB.ExtendedPDBReader

# function to keep compatible with the current API; should be cleaned up...
def parse(filename):
    """Parse atom information from extended PDB file *filename*.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                 :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.

    """
    return ExtendedPDBParser(filename).parse()

def parse_bonds(filename):
    """Parse atom information from extended PDB file *filename* and guesses bonds.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                 :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.
    """
    return ExtendedPDBParser(filename, guess_bonds_mode=True).parse()
