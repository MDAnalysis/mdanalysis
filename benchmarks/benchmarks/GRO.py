from __future__ import division, absolute_import, print_function

import numpy as np
from MDAnalysis.coordinates.GRO import GROReader
from MDAnalysisTests.datafiles import GRO

class GROReadBench(object):

    def time_read_GRO_file(self):
        """Benchmark reading of standard test
        suite GRO file.
        """
        GROReader(GRO)
