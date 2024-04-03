# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
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
:mod:`~MDAnalysis.topology.PDBParser` is that this parser reads a
non-standard PDB-like format in which residue numbers can be five digits
instead of four.

The topology reader reads a PDB file line by line and ignores atom numbers but
reads residue numbers up to 99,999 correctly. If you have systems containing at
least 100,000 residues then you need to use a different file format that can
handle such residue numbers.

See Also
--------
* :mod:`MDAnalysis.topology.PDBParser`
* :class:`MDAnalysis.coordinates.PDB.ExtendedPDBReader`
* :class:`MDAnalysis.core.universe.Universe`


Classes
-------

.. autoclass:: ExtendedPDBParser
   :members:
   :inherited-members:

"""

from . import PDBParser


class ExtendedPDBParser(PDBParser.PDBParser):
    """Parser that handles non-standard "extended" PDB file.

    Extended PDB files (MDAnalysis format specifier *XPDB*) may contain residue
    sequence numbers up to 99,999 by utilizing the insertion character field of
    the PDB standard.

    Creates a Topology with the following Attributes (if present):
     - serials
     - names
     - altLocs
     - chainids
     - tempfactors
     - occupancies
     - resids
     - resnames
     - segids
     - elements
     - bonds
     - formalcharges

    Guesses the following Attributes:
     - masses

    See Also
    --------
    :class:`MDAnalysis.coordinates.PDB.ExtendedPDBReader`

    .. versionadded:: 0.8
    """
    format = 'XPDB'
