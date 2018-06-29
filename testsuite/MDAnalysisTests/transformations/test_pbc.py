# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from __future__ import absolute_import

import numpy as np
import pytest

import MDAnalysis as mda
from MDAnalysis.transformations import unwrap
from __builtin__ import AttributeError

class TestUnwrap(object):

    @pytest.fixture()
    def universe(self):
        universe = mda.Universe(Make_Whole, Make_Whole)
        bondlist = [(0, 1), (1, 2), (1, 3), (1, 4), (4, 5), (4, 6), (4, 7)]
        universe.add_TopologyAttr(Bonds(bondlist))
        return universe

    @staticmethod
    @pytest.fixture()
    def ag(universe):
        return universe.residues[0].atoms
    
    def test_not_atomgroup(self):
        ag = 1
        with pytest.raises(AttributeError):
            unwrap(ag)(universe.trajectory.ts)
