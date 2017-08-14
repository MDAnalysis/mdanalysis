# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

import cclib
import numpy as np

from . import base
from ..core.topology import Topology
from ..core.topologyattrs import (
    AtomAttr,
    Atomids,
    Resids,
    Segids,
)


class MullikenCharges(AtomAttr):
    attrname = 'mulliken_charges'
    singular = 'mulliken_charge'
    per_object = 'atom'


class LowdinCharges(AtomAttr):
    attrname = 'lowdin_charges'
    singular = 'lowdin_charge'
    per_object = 'atom'


class CHELPGCharges(AtomAttr):
    attrname = 'CHELPG_charges'
    singular = 'CHELPG_charge'
    per_object = 'atom'


class CCLibParser(base.TopologyReaderBase):
    format = ['ORCA', 'GAUSSIAN', 'QCHEM']

    def parse(self):
        # first, obtain a CCLib object
        if isinstance(self.filename, cclib.parser.data.ccData_optdone_bool):
            ccobj = self.filename
        else:
            ccobj = cclib.ccopen(self.filename).parse()

        natoms = len(ccobj.atomnos)

        attrs = [
            Atomids(np.arange(natoms) + 1),
            Resids(np.array([1])),
            Segids(np.array(['SYSTEM'], dtype=object)),
        ]
        if hasattr(ccobj, 'atomcharges'):
            if 'mulliken' in ccobj.atomcharges:
                attrs.append(
                    MullikenCharges(np.array(ccobj.atomcharges['mulliken'])))
            if 'lowdin' in ccobj.atomcharges:
                attrs.append(
                    LowdinCharges(np.array(ccobj.atomcharges['lowdin'])))
            if 'chelpg' in ccobj.atomcharges:
                attrs.append(
                    CHELPGCharges(np.array(ccobj.atomcharges['chelpg'])))

        return Topology(natoms, 1, 1, attrs=attrs)
