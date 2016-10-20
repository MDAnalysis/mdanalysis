# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
from numpy.testing import (
    assert_,
    assert_array_equal,
)

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    DLP_CONFIG,
    DLP_CONFIG_order,
    DLP_CONFIG_minimal,
    DLP_HISTORY,
    DLP_HISTORY_order,
    DLP_HISTORY_minimal,
)


class DLPUniverse(ParserBase):
    def test_creates_universe(self):
        u = mda.Universe(self.filename, topology_format=self.format)
        assert_(isinstance(u, mda.Universe))
        
class DLPBase2(DLPUniverse):
    expected_attrs = ['ids', 'names']
    guessed_attrs = ['elements', 'masses']
    expected_n_atoms = 216
    expected_n_residues = 1
    expected_n_segments = 1

    def test_names(self):
        assert_(self.top.names.values[0] == 'K+')
        assert_(self.top.names.values[4] == 'Cl-')


class TestDLPHistoryParser(DLPBase2):
    parser = mda.topology.DLPolyParser.HistoryParser
    filename = DLP_HISTORY
    format='HISTORY'


class TestDLPConfigParser(DLPBase2):
    parser = mda.topology.DLPolyParser.ConfigParser
    filename = DLP_CONFIG
    format='CONFIG'


class DLPBase(DLPUniverse):
    expected_attrs = ['ids', 'names']
    guessed_attrs = ['elements', 'masses']
    expected_n_atoms = 3
    expected_n_residues = 1
    expected_n_segments = 1

    def test_dlp_names(self):
        assert_array_equal(self.top.names.values,
                           ['C', 'B', 'A'])


class TestDLPConfigMinimal(DLPBase):
    parser = mda.topology.DLPolyParser.ConfigParser
    filename = DLP_CONFIG_minimal
    format='CONFIG'


class TestDLPConfigOrder(DLPBase):
    parser = mda.topology.DLPolyParser.ConfigParser
    filename = DLP_CONFIG_order
    format='CONFIG'


class TestDLPHistoryMinimal(DLPBase):
    parser = mda.topology.DLPolyParser.HistoryParser
    filename = DLP_HISTORY_minimal
    format='HISTORY'


class TestDLPHistoryOrder(DLPBase):
    parser = mda.topology.DLPolyParser.HistoryParser
    filename = DLP_HISTORY_order
    format='HISTORY'
