from MDAnalysis.coordinates.PQR import PQRReader, PQRWriter
from MDAnalysisTests.coordinates.base import BaseReaderTest, BaseReference
from MDAnalysisTests.datafiles import COORDINATES_PQR


class PQRReference(BaseReference):
    def __init__(self):
        super(PQRReference, self).__init__()
        self.trajectory = COORDINATES_PQR
        self.topology = COORDINATES_PQR
        self.reader = PQRReader
        self.writer = PQRWriter
        self.ext = 'pqr'
        self.prec = 3
        self.n_frames = 1
        self.totaltime = 0


class TestPQRReader(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = PQRReference()
        super(TestPQRReader, self).__init__(reference)
