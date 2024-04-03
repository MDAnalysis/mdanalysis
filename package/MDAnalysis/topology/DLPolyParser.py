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

"""DL Poly format Topology Readers --- :mod:`MDAnalysis.topology.DLPolyParser`
==============================================================================

Read DL Poly_ format topology files

DLPoly files have the following Attributes:
 - Atomnames
 - Atomids
Guesses the following attributes:
 - Atomtypes
 - Masses

.. _Poly: http://www.stfc.ac.uk/SCD/research/app/ccg/software/DL_POLY/44516.aspx

Classes
-------

.. autoclass:: ConfigParser

.. autoclass:: HistoryParser


"""
import numpy as np

from . import guessers
from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Masses,
    Resids,
    Resnums,
    Segids,
)
from ..lib.util import openany


class ConfigParser(TopologyReaderBase):
    """DL_Poly CONFIG file parser

    .. versionadded:: 0.10.1
    """
    format = 'CONFIG'

    def parse(self, **kwargs):
        with openany(self.filename) as inf:
            inf.readline()
            levcfg, imcon, megatm = np.int64(inf.readline().split()[:3])
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

        atomtypes = guessers.guess_types(names)
        masses = guessers.guess_masses(atomtypes)

        attrs = [
            Atomnames(names),
            Atomids(ids),
            Atomtypes(atomtypes, guessed=True),
            Masses(masses, guessed=True),
            Resids(np.array([1])),
            Resnums(np.array([1])),
            Segids(np.array(['SYSTEM'], dtype=object)),
        ]
        top = Topology(n_atoms, 1, 1,
                       attrs=attrs)

        return top


class HistoryParser(TopologyReaderBase):
    """DL_Poly History file parser

    .. versionadded:: 0.10.1
    """
    format = 'HISTORY'

    def parse(self, **kwargs):
        with openany(self.filename) as inf:
            inf.readline()
            levcfg, imcon, megatm = np.int64(inf.readline().split()[:3])

            names = []
            ids = []

            line = inf.readline()
            while not (len(line.split()) == 4 or len(line.split()) == 5):
                line = inf.readline()
                if line == '':
                    raise EOFError("End of file reached when reading HISTORY.")

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

        atomtypes = guessers.guess_types(names)
        masses = guessers.guess_masses(atomtypes)

        attrs = [
            Atomnames(names),
            Atomids(ids),
            Atomtypes(atomtypes, guessed=True),
            Masses(masses, guessed=True),
            Resids(np.array([1])),
            Resnums(np.array([1])),
            Segids(np.array(['SYSTEM'], dtype=object)),
        ]
        top = Topology(n_atoms, 1, 1,
                       attrs=attrs)

        return top
