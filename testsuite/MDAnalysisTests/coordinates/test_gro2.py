import numpy as np
from MDAnalysis.coordinates.GRO import GROReader, GROWriter
from MDAnalysisTests.coordinates.base import BaseReference, BaseReaderTest, BaseWriterTest
from MDAnalysisTests.datafiles import COORDINATES_GRO


class GROReference(BaseReference):
    def __init__(self):
        super(GROReference, self).__init__()
        self.trajectory = COORDINATES_GRO
        self.topology = COORDINATES_GRO
        self.reader = GROReader
        self.writer = GROWriter
        self.ext = 'gro'
        self.n_frames = 1
        self.prec = 4
        self.first_frame.velocities = np.array(
            [[0.0000, 0.100, 0.200],
             [0.300, 0.400, 0.500],
             [0.600, 0.700, 0.800],
             [0.900, 1.000, 1.100],
             [1.200, 1.300, 1.400]],
            dtype=np.float32)
        self.totaltime = 0
        self.container_format = True


class TestGROReader(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROReference()
        super(TestGROReader, self).__init__(reference)


class TestGROWriter(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROReference()
        super(TestGROWriter, self).__init__(reference)
