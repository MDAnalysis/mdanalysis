# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
TNG Topology Parser
=========================================================================

This topology parser uses a standard TNG file to build a minimum
internal structure representation (list of atoms).

See Also
--------
* :class:`MDAnalysis.coordinates.TNG.TNGReader`
"""

from __future__ import absolute_import, print_function

import numpy as np
import warnings
import pytng

from .guessers import guess_masses, guess_types
from ..lib import util
from .base import TopologyReaderBase, change_squash
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    AltLocs,
    Bonds,
    ChainIDs,
    Atomtypes,
    ICodes,
    Masses,
    Occupancies,
    Resids,
    Resnames,
    Resnums,
    Segids,
    Tempfactors,
)


class TNGParser(TopologyReaderBase):
    """Parser that obtains a list of atoms from a standard TNG file.

    Creates the following Attributes:
     - names
     - chainids
     - resids
     - resnames

    Guesses the following Attributes:
     - masses


    .. versionadded:: 0.16.1
    """
    format = ['tng']

    def parse(self):
        """Parse atom information from tng file

        Returns
        -------
        MDAnalysis Topology object
        """
        with pytng.TNGFile(self.filename) as tng:
            n_atoms = tng.n_atoms
            residx = tng.residues_ids
            n_residues = tng.n_residues
            n_segments = tng.n_chains


        top = Topology(n_atoms, n_residues, n_segments,
                       atom_resindex=residx)

        return top
