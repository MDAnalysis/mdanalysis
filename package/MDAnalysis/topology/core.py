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
Common functions for topology building --- :mod:`MDAnalysis.topology.core`
==========================================================================

The various topology parsers make use of functions and classes in this
module. They are mostly of use to developers.

.. SeeAlso:: :mod:`MDAnalysis.topology.tables` for some hard-coded atom
   information that is used by functions such as :func:`guess_atom_type` and
   :func:`guess_atom_mass`.

"""

from __future__ import print_function
import six

# Global imports
import os.path
import numpy as np
from collections import defaultdict

# Local imports
from . import _PARSERS
from . import tables
from .guessers import (
    guess_atom_element, guess_atom_type,
    get_atom_mass, guess_atom_mass, guess_atom_charge,
    guess_bonds, guess_angles, guess_dihedrals, guess_improper_dihedrals,
)
from ..lib.util import cached
from ..lib import util


def get_parser_for(filename, format=None):
    """Return the appropriate topology parser for *filename*.

    Automatic detection is disabled when an explicit *format* is
    provided.

    :Raises:
      *ValueError*
        If no appropriate parser could be found.
    """
    if format is None:
        format = util.guess_format(filename)
    format = format.upper()
    try:
        return _PARSERS[format]
    except KeyError:
        raise ValueError(
            "Cannot autodetect topology type for file '{0}' "
            "(file extension could not be parsed).\n"
            "           You can use 'Universe(topology, ..., topology_format=FORMAT)' "
            "to explicitly specify the format and\n"
            "           override automatic detection. Known FORMATs are:\n"
            "           {1}\n"
            "           See http://docs.mdanalysis.org/documentation_pages/topology/init.html#supported-topology-formats\n"
            "           For missing formats, raise an issue at "
            "http://issues.mdanalysis.org".format(filename, _PARSERS.keys()))


