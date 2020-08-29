from __future__ import division, absolute_import, print_function

import numpy as np
from MDAnalysis.coordinates.GRO import GROReader
from MDAnalysis.topology.GROParser import GROParser
from MDAnalysisTests.datafiles import GRO
import MDAnalysis as mda

class GROReadBench(object):
    def time_read_GRO_coordinates(self):
        """Benchmark reading of standard testsuite GRO file."""
        GROReader(GRO)

    def time_parse_GRO_file(self):
        """Time to create topology from GRO file"""
        p = GROParser(GRO)
        top = p.parse()

    def time_create_GRO_universe(self):
        """Time to create MDA Universe of GRO"""
        u = mda.Universe(GRO)
