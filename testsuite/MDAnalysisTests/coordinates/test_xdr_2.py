import MDAnalysis as mda

from MDAnalysisTests.datafiles import (COORDINATES_XTC, COORDINATES_TOPOLOGY,
                                       COORDINATES_TRR)
from MDAnalysisTests.coordinates.base import (BaseReaderTest, BaseReference,
                                              BaseWriterTest)

from MDAnalysis.tests.datafiles import XTC, TRR
import MDAnalysis.core.AtomGroup
from MDAnalysis.coordinates import XDR

import numpy as np
from numpy.testing import (assert_array_almost_equal, TestCase,
                           assert_equal, dec)

import os
import shutil
import tempdir


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


class _GromacsReader_offsets(TestCase):
    # This base class assumes same lengths and dt for XTC and TRR test cases!
    filename = None
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.],
                            dtype=np.float32)
    # computed with Gromacs: 362.26999999999998 nm**3 * 1000 A**3/nm**3
    ref_volume = 362270.0
    ref_offsets = None
    _reader = None

    def setUp(self):
        # since offsets are automatically generated in the same directory
        # as the trajectory, we do everything from a temporary directory
        self.tmpdir = tempdir.TempDir()
        shutil.copy(self.filename, self.tmpdir.name)

        self.traj = os.path.join(self.tmpdir.name,
                                 os.path.basename(self.filename))

        self.trajectory = self._reader(self.traj)
        self.prec = 3
        self.ts = self.trajectory.ts

    def tearDown(self):
        del self.tmpdir
        del self.trajectory

    @dec.slow
    def test_offsets(self):
        self.trajectory._read_offsets(store=True)
        assert_array_almost_equal(self.trajectory._xdr.offsets,
                                  self.ref_offsets,
                                  err_msg="wrong frame offsets")

        outfile_offsets = XDR.offsets_filename(self.traj)
        with open(outfile_offsets) as f:
            saved_offsets = {k: v for k, v in np.load(f).iteritems()}

        assert_array_almost_equal(self.trajectory._xdr.offsets,
                                  saved_offsets['offsets'],
                                  err_msg="error saving frame offsets")
        assert_array_almost_equal(self.ref_offsets, saved_offsets['offsets'],
                                  err_msg="saved frame offsets don't match "
                                  "the known ones")

        self.trajectory._load_offsets()
        assert_array_almost_equal(self.trajectory._xdr.offsets,
                                  self.ref_offsets,
                                  err_msg="error loading frame offsets")
        assert_equal(saved_offsets['ctime'], os.path.getctime(self.traj))
        assert_equal(saved_offsets['size'], os.path.getsize(self.traj))

    # TODO: tests mismatchs
    # @dec.slow
    # def test_persistent_offsets_size_mismatch(self):
    #     # check that stored offsets are not loaded when trajectory
    #     # size differs from stored size
    #     with open(XDR.offsets_filename(self.traj), 'rb') as f:
    #         saved_offsets = {k: v for k, v in np.load(f).iteritems()}
    #     saved_offsets['size'] += 1
    #     with open(XDR.offsets_filename(self.traj), 'wb') as f:
    #         np.savez(f, **saved_offsets)

    # TODO: This doesn't test if the offsets work AT ALL. the old
    # implementation only checked if the offsets were ok to set back to the old
    # frame. But that doesn't check if any of the other offsets is potentially
    # wrong. Basically the only way to check that would be to scan through the
    # whole trajectory.
    # @dec.slow
    # def test_persistent_offsets_last_frame_wrong(self):
    #     # check that stored offsets are not loaded when the offsets
    #     # themselves appear to be wrong
    #     with open(XDR.offsets_filename(self.traj), 'rb') as f:
    #         saved_offsets = {k: v for k, v in np.load(f).iteritems()}
    #     saved_offsets['offsets'] += 1
    #     with open(XDR.offsets_filename(self.traj), 'wb') as f:
    #         np.savez(f, **saved_offsets)

    #     # with warnings.catch_warnings():
    #     #     u = MDAnalysis.Universe(self.top, self.traj)
    #     #     assert_equal((u.trajectory._xdr.offsets is None), True)

    @dec.slow
    def test_persistent_offsets_readonly(self):
        os.remove(XDR.offsets_filename(self.traj))
        assert_equal(os.path.exists(
            XDR.offsets_filename(self.trajectory.filename)), False)

        os.chmod(self.tmpdir.name, 0555)
        self.trajectory._read_offsets(store=True)
        assert_equal(os.path.exists(
            XDR.offsets_filename(self.trajectory.filename)), False)


class TestXTCReader_offsets(_GromacsReader_offsets):
    filename = XTC
    ref_offsets = np.array([0, 165188, 330364, 495520, 660708, 825872, 991044,
                            1156212, 1321384, 1486544])
    _reader = MDAnalysis.coordinates.XTC.XTCReader


class TestTRRReader_offsets(_GromacsReader_offsets):
    filename = TRR
    ref_offsets = np.array([0, 1144464, 2288928, 3433392, 4577856, 5722320,
                            6866784, 8011248, 9155712, 10300176])
    _reader = MDAnalysis.coordinates.TRR.TRRReader
