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
from numpy.testing import assert_equal
import pytest

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    DLP_CONFIG,
    DLP_CONFIG_order,
    DLP_CONFIG_minimal,
    DLP_HISTORY,
    DLP_HISTORY_order,
    DLP_HISTORY_minimal,
    DLP_HISTORY_minimal_cell,
    DLP_HISTORY_classic
)


class DLPUniverse(ParserBase):
    def test_creates_universe(self, filename):
        u = mda.Universe(filename, topology_format=self.format)
        assert isinstance(u, mda.Universe)


class DLPBase2(DLPUniverse):
    expected_attrs = ['ids', 'names']
    guessed_attrs = ['types', 'masses']
    expected_n_atoms = 216
    expected_n_residues = 1
    expected_n_segments = 1

    def test_names(self, top):
        assert top.names.values[0] == 'K+'
        assert top.names.values[4] == 'Cl-'


class TestDLPHistoryParser(DLPBase2):
    parser = mda.topology.DLPolyParser.HistoryParser
    ref_filename = DLP_HISTORY
    format = 'HISTORY'


class TestDLPConfigParser(DLPBase2):
    parser = mda.topology.DLPolyParser.ConfigParser
    ref_filename = DLP_CONFIG
    format = 'CONFIG'


class DLPBase(DLPUniverse):
    expected_attrs = ['ids', 'names']
    guessed_attrs = ['types', 'masses']
    expected_n_atoms = 3
    expected_n_residues = 1
    expected_n_segments = 1

    def test_dlp_names(self, top):
        assert_equal(top.names.values,
                     ['C', 'B', 'A'])


class TestDLPConfigMinimal(DLPBase):
    parser = mda.topology.DLPolyParser.ConfigParser
    ref_filename = DLP_CONFIG_minimal
    format = 'CONFIG'


class TestDLPConfigOrder(DLPBase):
    parser = mda.topology.DLPolyParser.ConfigParser
    ref_filename = DLP_CONFIG_order
    format = 'CONFIG'


class TestDLPHistoryMinimal(DLPBase):
    parser = mda.topology.DLPolyParser.HistoryParser
    ref_filename = DLP_HISTORY_minimal
    format = 'HISTORY'


class TestDLPHistoryMinimal(DLPBase):
    parser = mda.topology.DLPolyParser.HistoryParser
    ref_filename = DLP_HISTORY_minimal_cell
    format = 'HISTORY'


class TestDLPHistoryOrder(DLPBase):
    parser = mda.topology.DLPolyParser.HistoryParser
    ref_filename = DLP_HISTORY_order
    format = 'HISTORY'


class TestDLPHistoryClassic(DLPBase):
    parser = mda.topology.DLPolyParser.HistoryParser
    ref_filename = DLP_HISTORY_classic
    format = 'HISTORY'


def test_HISTORY_EOFError():
    with pytest.raises(EOFError):
        mda.Universe(DLP_CONFIG, topology_format='HISTORY')
