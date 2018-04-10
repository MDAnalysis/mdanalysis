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

import pytest

import MDAnalysis as mda
import numpy as np
from numpy.testing import assert_almost_equal

from MDAnalysis.coordinates.CHEMFILES import ChemfilesReader

from MDAnalysisTests.datafiles import COORDINATES_XYZ
from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest, BaseWriterTest, BaseReference
)


class ChemfilesXYZReference(BaseReference):
    def __init__(self):
        super(ChemfilesXYZReference, self).__init__()
        self.trajectory = COORDINATES_XYZ
        self.topology = COORDINATES_XYZ
        self.reader = mda.coordinates.CHEMFILES.ChemfilesReader
        self.writer = mda.coordinates.CHEMFILES.ChemfilesWriter
        self.ext = 'xyz'
        self.volume = 0
        self.dimensions = np.zeros(6)
        self.dimensions[3:] = 90.0


class TestChemfilesReader(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return ChemfilesXYZReference()


class TestChemfilesWriter(BaseWriterTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return ChemfilesXYZReference()

    # Disable 'test_no_container' as it try to open a file for writing without
    # extension.
    def test_no_container(self, ref):
        pass
