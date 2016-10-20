# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""DL Poly format Topology Readers --- :mod:`MDAnalysis.topology.DLPolyParser`
==============================================================================

Read DL Poly_ format topology files

DLPoly files have the following Attributes:
 - Atomnames
 - Atomids
Guesses the following attributes:
 - Elements
 - Masses

.. _Poly: http://www.stfc.ac.uk/SCD/research/app/ccg/software/DL_POLY/44516.aspx
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy as np

from . import guessers
from .base import TopologyReader
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Elements,
    Masses,
    Resids,
    Resnums,
    Segids,
)
from ..lib.util import openany


class ConfigParser(TopologyReader):
    """DL_Poly CONFIG file parser

    .. versionadded:: 0.10.1
    """
    format = 'CONFIG'

    def parse(self):
        with openany(self.filename, 'r') as inf:
            inf.readline()
            levcfg, imcon, megatm = map(int,
                                        inf.readline().split()[:3])
            if not imcon == 0:
                inf.readline()
                inf.readline()
                inf.readline()

            names = []
            ids = []

            line = inf.readline().strip()
            while line:
                name = line[:8].strip()
                names.append(name)
                try:
                    idx = int(line[8:])
                except ValueError:
                    pass
                else:
                    ids.append(idx)

                inf.readline()
                if levcfg > 0:
                    inf.readline()
                if levcfg == 2:
                    inf.readline()

                line = inf.readline()

        n_atoms = len(names)
        if ids:
            ids = np.array(ids)
            names = np.array(names, dtype=object)
            order = np.argsort(ids)
            ids = ids[order]
            names = names[order]
        else:
            ids = np.arange(n_atoms)

        elements = guessers.guess_types(names)
        masses = guessers.guess_masses(elements)

        attrs = [
            Atomnames(names),
            Atomids(ids),
            Elements(elements, guessed=True),
            Masses(masses, guessed=True),
            Resids(np.array([1])),
            Resnums(np.array([1])),
            Segids(np.array(['SYSTEM'], dtype=object)),
        ]
        top = Topology(n_atoms, 1, 1,
                       attrs=attrs)

        return top


class HistoryParser(TopologyReader):
    """DL_Poly History file parser

    .. versionadded:: 0.10.1
    """
    format = 'HISTORY'

    def parse(self):
        with openany(self.filename, 'r') as inf:
            inf.readline()
            levcfg, imcon, megatm = map(int, inf.readline().split()[:3])
            inf.readline()
            if not imcon == 0:
                inf.readline()
                inf.readline()
                inf.readline()

            names = []
            ids = []

            line = inf.readline()
            while line and not line.startswith('timestep'):
                name = line[:8].strip()
                names.append(name)
                try:
                    idx = int(line.split()[1])
                except ValueError:
                    pass
                else:
                    ids.append(idx)

                inf.readline()
                if levcfg > 0:
                    inf.readline()
                if levcfg == 2:
                    inf.readline()

                line = inf.readline()

        n_atoms = len(names)
        if ids:
            ids = np.array(ids)
            names = np.array(names, dtype=object)
            order = np.argsort(ids)
            ids = ids[order]
            names = names[order]
        else:
            ids = np.arange(n_atoms)

        elements = guessers.guess_types(names)
        masses = guessers.guess_masses(elements)
            
        attrs = [
            Atomnames(names),
            Atomids(ids),
            Elements(elements, guessed=True),
            Masses(masses, guessed=True),
            Resids(np.array([1])),
            Resnums(np.array([1])),
            Segids(np.array(['SYSTEM'], dtype=object)),
        ]
        top = Topology(n_atoms, 1, 1,
                       attrs=attrs)

        return top
