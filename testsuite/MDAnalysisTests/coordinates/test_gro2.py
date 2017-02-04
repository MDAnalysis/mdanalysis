from MDAnalysis.coordinates.GRO import GROReader, GROWriter
from MDAnalysisTests.coordinates.base import BaseReference, BaseReaderTest, MultiframeReaderTest
from MDAnalysisTests.datafiles import COORDINATES_GRO
import numpy as np


class GROReference(BaseReference):
    def __init__(self):
        super(GROReference, self).__init__()
        self.trajectory = COORDINATES_GRO
        # self.topology = COORDINATES_GRO
        self.reader = GROReader
        self.writer = GROWriter
        self.ext = 'gro'
        self.n_frames = 1
        # self.volume = 0
        # self.dimensions = np.zeros(6)
        # self.container_format = True


class TestGROReader(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROReference()
        super(TestGROReader, self).__init__(reference)
