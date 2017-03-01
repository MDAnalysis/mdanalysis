import numpy as np
from MDAnalysis.coordinates.GRO import GROReader, GROWriter
from MDAnalysisTests.coordinates.base import BaseReference, BaseReaderTest, BaseWriterTest
from MDAnalysisTests.datafiles import COORDINATES_GRO, COORDINATES_GRO_INCOMPLETE_VELOCITY, COORDINATES_GRO_BZ2


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


class GRONoConversionReference(GROReference):
    def __init__(self):
        super(GRONoConversionReference, self).__init__()
        self.first_frame.positions /= 10.0
        self.first_frame.velocities /= 10.0
        self.dimensions[:3] /= 10.0
        self.volume /= 1000


class TestGROReaderNoConversion(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GRONoConversionReference()
        super(TestGROReaderNoConversion, self).__init__(reference)
        self.reader = self.ref.reader(self.ref.trajectory, convert_units=False)


class TestGROWriterNoConversion(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GRONoConversionReference()
        super(TestGROWriterNoConversion, self).__init__(reference)
        self.writer = self.ref.writer(self.ref.trajectory, convert_units=False)


class GROReaderIncompleteVelocitiesReference(GROReference):
    def __init__(self):
        super(GROReaderIncompleteVelocitiesReference, self).__init__()
        self.trajectory = COORDINATES_GRO_INCOMPLETE_VELOCITY
        self.topology = COORDINATES_GRO_INCOMPLETE_VELOCITY
        self.first_frame.velocities = np.array(
            [[0.0000, 0.100, 0.200],
             [0.000, 0.000, 0.000],
             [0.600, 0.700, 0.800],
             [0.900, 1.000, 1.100],
             [1.200, 1.300, 1.400]],
            dtype=np.float32)


class TestGROReaderIncompleteVelocities(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROReaderIncompleteVelocitiesReference()
        super(TestGROReaderIncompleteVelocities, self).__init__(reference)


class TestGROWriterIncompleteVelocities(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROReaderIncompleteVelocitiesReference()
        super(TestGROWriterIncompleteVelocities, self).__init__(reference)


class GROBZReference(GROReference):
    def __init__(self):
        super(GROBZReference, self).__init__()
        self.trajectory = COORDINATES_GRO_BZ2
        self.topology = COORDINATES_GRO_BZ2


class TestGROBZ2Reader(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROBZReference()
        super(TestGROBZ2Reader, self).__init__(reference)


class TestGROBZ2Writer(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = GROBZReference()
        super(TestGROBZ2Writer, self).__init__(reference)
