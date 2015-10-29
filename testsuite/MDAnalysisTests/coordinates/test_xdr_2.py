import MDAnalysis as mda

from MDAnalysisTests.datafiles import (COORDINATES_XTC, COORDINATES_TOPOLOGY,
                                       COORDINATES_TRR)
from MDAnalysisTests.coordinates.base import (BaseReaderTest, BaseReference,
                                              BaseWriterTest)


class XTCReference(BaseReference):
    def __init__(self):
        super(XTCReference, self).__init__()
        self.trajectory = COORDINATES_XTC
        self.topology = COORDINATES_TOPOLOGY
        self.reader = mda.coordinates.XTC.XTCReader
        self.writer = mda.coordinates.XTC.XTCWriter
        self.ext = 'xtc'
        self.prec = 3


class TestXTCReader(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = XTCReference()
        super(TestXTCReader, self).__init__(reference)


class TestXTCWriter(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = XTCReference()
        super(TestXTCWriter, self).__init__(reference)


class TRRReference(BaseReference):
    def __init__(self):
        super(TRRReference, self).__init__()
        self.trajectory = COORDINATES_TRR
        self.topology = COORDINATES_TOPOLOGY
        self.reader = mda.coordinates.TRR.TRRReader
        self.writer = mda.coordinates.TRR.TRRWriter
        self.ext = 'xtc'
        self.prec = 3
        self.first_frame.velocities = self.first_frame.positions / 10
        self.first_frame.forces = self.first_frame.positions / 100

        self.second_frame.velocities = self.second_frame.positions / 10
        self.second_frame.forces = self.second_frame.positions / 100

        self.last_frame.velocities = self.last_frame.positions / 10
        self.last_frame.forces = self.last_frame.positions / 100

        self.jump_to_frame.velocities = self.jump_to_frame.positions / 10
        self.jump_to_frame.forces = self.jump_to_frame.positions / 100

    def iter_ts(self, i):
        ts = self.first_frame.copy()
        ts.positions = 2**i * self.first_frame.positions
        ts.velocities = ts.positions / 10
        ts.forces = ts.positions / 100
        ts.time = i
        ts.frame = i
        return ts


class TestTRRReader(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = TRRReference()
        super(TestTRRReader, self).__init__(reference)


class TestTRRWriter(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = TRRReference()
        super(TestTRRWriter, self).__init__(reference)
