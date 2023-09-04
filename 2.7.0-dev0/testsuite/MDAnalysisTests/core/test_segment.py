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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from numpy.testing import (
    assert_equal,
)
import pytest
import pickle

import MDAnalysis as mda

from MDAnalysisTests import make_Universe
from MDAnalysis.tests.datafiles import PSF, DCD


class TestSegment(object):
    @pytest.fixture()
    def universe(self):
        return make_Universe(('segids',))

    @pytest.fixture()
    def sB(self, universe):
        return universe.segments[1]

    def test_type(self, sB):
        assert isinstance(sB, mda.core.groups.Segment)
        assert_equal(sB.segid, "SegB")

    def test_index(self, sB):
        s = sB
        res = s.residues[3]
        assert isinstance(res, mda.core.groups.Residue)

    def test_slicing(self, sB):
        res = sB.residues[:3]
        assert_equal(len(res), 3)
        assert isinstance(res, mda.core.groups.ResidueGroup)

    def test_advanced_slicing(self, sB):
        res = sB.residues[[2, 1, 0, 2]]
        assert_equal(len(res), 4)
        assert isinstance(res, mda.core.groups.ResidueGroup)

    def test_atom_order(self, universe):
        assert_equal(universe.segments[0].atoms.indices,
                     sorted(universe.segments[0].atoms.indices))

    @pytest.mark.parametrize("ix", (1, -1))
    def test_residue_pickle(self, universe, ix):
        seg_out = universe.segments[ix]
        seg_in = pickle.loads(pickle.dumps(seg_out))
        assert seg_in == seg_out
