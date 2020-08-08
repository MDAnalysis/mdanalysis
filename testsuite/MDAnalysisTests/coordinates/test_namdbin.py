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
import numpy as np
import os

import pytest
from numpy.testing import (assert_allclose,
                           assert_equal,
                           assert_almost_equal)

import MDAnalysis as mda
from MDAnalysisTests.datafiles import NAMDBIN, PDB_small
from MDAnalysisTests.coordinates.base import (_SingleFrameReader,
                                              BaseReference,
                                              BaseWriterTest)


class TestNAMDBINReader(_SingleFrameReader):
    """Test reading namd binary coordinate file"""
    __test__ = True

    def setUp(self):
        # can lead to race conditions when testing in parallel
        self.universe = mda.Universe(PDB_small, NAMDBIN)
        # 6 decimals in NAMDBIN spec
        self.prec = 6
        self.ref_n_atoms = 3341

    def test_parse_n_atoms(self):
        n_atoms = mda.coordinates.NAMDBIN.NAMDBINReader.parse_n_atoms(NAMDBIN)
        assert_equal(n_atoms, self.ref_n_atoms)

    def test_get_writer_from_reader(self):
        universe = mda.Universe(PDB_small, NAMDBIN)
        writer = universe.trajectory.Writer('NAMDBIN-test')
        assert isinstance(writer,
                          mda.coordinates.NAMDBIN.NAMDBINWriter)


class NAMDBINReference(BaseReference):
    def __init__(self):
        super(NAMDBINReference, self).__init__()
        self.trajectory = NAMDBIN
        self.topology = PDB_small
        self.reader = mda.coordinates.NAMDBIN.NAMDBINReader
        self.writer = mda.coordinates.NAMDBIN.NAMDBINWriter
        self.ext = 'coor'
        self.volume = 0
        self.dimensions = np.zeros(6)
        self.container_format = True


class NAMDBINWriter(BaseWriterTest):
    __test__ = True
    @staticmethod
    @pytest.fixture()
    def ref():
        return NAMDBINReference()
