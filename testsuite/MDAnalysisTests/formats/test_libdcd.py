from __future__ import print_function

from nose.tools import raises
from numpy.testing import assert_equal, assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal

from MDAnalysis.lib.formats.libdcd import DCDFile
from MDAnalysisTests.datafiles import (PSF, DCD, DCD_NAMD_TRICLINIC,
                                       PSF_NAMD_TRICLINIC)

from unittest import TestCase
import MDAnalysis
from MDAnalysisTests.tempdir import run_in_tempdir
from MDAnalysisTests import tempdir
import numpy as np
import os

class DCDReadFrameTest(TestCase):

    def setUp(self):
        self.dcdfile = DCDFile(DCD)
        self.natoms = 3341
        self.traj_length = 98
        self.new_frame = 91
        self.context_frame = 22
        self.num_iters = 3
        self.expected_remarks = '''* DIMS ADK SEQUENCE FOR PORE PROGRAM                                            * WRITTEN BY LIZ DENNING (6.2008)                                               *  DATE:     6/ 6/ 8     17:23:56      CREATED BY USER: denniej0                '''

    def tearDown(self):
        del self.dcdfile

    def test_header_remarks(self):
        # confirm correct header remarks section reading
        with self.dcdfile as f:
            list_chars = []
            for element in f.remarks:
                list_chars.append(element)

            list_chars = []
            for element in self.expected_remarks:
                list_chars.append(element)
            self.assertEqual(len(f.remarks), len(self.expected_remarks))

    def test_read_coords(self):
        # confirm shape of coordinate data against result from previous
        # MDAnalysis implementation of DCD file handling
        dcd_frame = self.dcdfile.read()
        xyz = dcd_frame[0]
        assert_equal(xyz.shape, (self.natoms, 3))

    def test_read_unit_cell(self):
        # confirm unit cell read against result from previous
        # MDAnalysis implementation of DCD file handling
        dcd_frame = self.dcdfile.read()
        unitcell = dcd_frame[1]
        expected = np.array([  0.,   0.,   0.,  90.,  90.,  90.],
                            dtype=np.float32)
        assert_equal(unitcell, expected)

    def test_seek_over_max(self):
        # should raise IOError if beyond 98th frame
        with self.assertRaises(IOError):
            self.dcdfile.seek(102)

    def test_seek_normal(self):
        # frame seek within range is tested
        self.dcdfile.seek(self.new_frame)
        assert_equal(self.dcdfile.tell(), self.new_frame)

    def test_seek_negative(self):
        # frame seek with negative number
        with self.assertRaises(IOError):
            self.dcdfile.seek(-78)

    def test_iteration(self):
        expected = 0
        while self.num_iters > 0:
            self.dcdfile.__next__()
            self.num_iters -= 1
            expected += 1
        
        assert_equal(self.dcdfile.tell(), expected)

    def test_zero_based_frames(self):
        expected_frame = 0
        assert_equal(self.dcdfile.tell(), expected_frame)

    def test_length_traj(self):
        expected = self.traj_length
        assert_equal(len(self.dcdfile), expected)

    def test_context_manager(self):
        with self.dcdfile as f:
            f.seek(self.context_frame)
            assert_equal(f.tell(), self.context_frame)

    @raises(IOError)
    def test_open_wrong_mode(self):
        DCDFile('foo', 'e')

    @raises(IOError)
    def test_raise_not_existing(self):
        DCDFile('foo')

    def test_n_atoms(self):
        assert_equal(self.dcdfile.n_atoms, self.natoms)

    @raises(IOError)
    @run_in_tempdir()
    def test_read_write_mode_file(self):
        with DCDFile('foo', 'w') as f:
            f.read()

    @raises(IOError)
    def test_read_closed(self):
        self.dcdfile.close()
        self.dcdfile.read()

    def test_iteration_2(self):
        with self.dcdfile as f:
            for frame in f:
                pass

class DCDWriteHeaderTest(TestCase):

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.dcdfile = DCDFile(self.testfile, 'w')
        self.dcdfile_r = DCDFile(DCD, 'r')

    def tearDown(self):
        try: 
            os.unlink(self.testfile)
        except OSError:
            pass
        del self.tmpdir
    
    def test_write_header_crude(self):
        # test that _write_header() can produce a very crude
        # header for a new / empty file
        self.dcdfile._write_header(remarks='Crazy!', n_atoms=22,
                                   starting_step=12, ts_between_saves=10,
                                   time_step=0.02)
        self.dcdfile.close()

        # we're not actually asserting anything, yet
        # run with: nosetests test_libdcd.py --nocapture
        # to see printed output from nose
        with open(self.testfile, "rb") as f:
            for element in f:
                print(element)

    def test_write_header_mode_sensitivy(self):
        # an exception should be raised on any attempt to use
        # _write_header with a DCDFile object in 'r' mode
        with self.assertRaises(IOError):
            self.dcdfile_r._write_header(remarks='Crazy!', n_atoms=22,
                                         starting_step=12, ts_between_saves=10,
                                         time_step=0.02)



class DCDWriteTest(TestCase):

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.dcdfile = DCDFile(self.testfile, 'w')
        self.dcdfile_r = DCDFile(DCD, 'r')
        self.natoms = 3341
        self.expected_frames = 98
        self.seek_frame = 91
        self.expected_remarks = '''* DIMS ADK SEQUENCE FOR PORE PROGRAM                                            * WRITTEN BY LIZ DENNING (6.2008)                                               *  DATE:     6/ 6/ 8     17:23:56      CREATED BY USER: denniej0                '''

        with self.dcdfile_r as f_in, self.dcdfile as f_out:
            for frame in f_in:
                frame = frame._asdict()
                f_out.write(xyz=frame['x'],
                            box=frame['unitcell'].astype(np.float64),
                            step=f_in.istart,
                            natoms=frame['x'].shape[0],
                            charmm=0,
                            time_step=f_in.delta,
                            ts_between_saves=f_in.nsavc,
                            remarks=f_in.remarks)


    def tearDown(self):
        try: 
            os.unlink(self.testfile)
        except OSError:
            pass
        del self.tmpdir

    def test_write_mode(self):
        # ensure that writing of DCD files only occurs with properly
        # opened files
        with self.assertRaises(IOError):
            self.dcdfile_r.write(xyz=np.zeros((3,3)),
                                 box=np.zeros(6, dtype=np.float64),
                                 step=0,
                                 natoms=330,
                                 charmm=0,
                                 time_step=22.2,
                                 ts_between_saves=3,
                                 remarks='')

    def test_written_dcd_coordinate_data_shape(self):
        # written coord shape should match for all frames
        expected = (self.natoms, 3)
        with DCDFile(self.testfile) as f:
            if f.n_frames > 1:
                for frame in f:
                    xyz = f.read()[0]
                    assert_equal(xyz.shape, expected)
            else:
                xyz = f.read()[0]
                assert_equal(xyz.shape, expected)

    def test_written_unit_cell(self):
        # written unit cell dimensions should match for all frames
        expected = np.array([  0.,   0.,   0.,  90.,  90.,  90.],
                            dtype=np.float32)
        with DCDFile(self.testfile) as f:
            if f.n_frames > 1:
                for frame in f:
                    unitcell = f.read()[1]
                    assert_equal(unitcell, expected)
            else:
                unitcell = f.read()[1]
                assert_equal(unitcell, expected)

    def test_written_num_frames(self):
        with DCDFile(self.testfile) as f:
            assert_equal(len(f), self.expected_frames)

    def test_written_seek(self):
        # ensure that we can seek properly on written DCD file
        with DCDFile(self.testfile) as f:
            f.seek(self.seek_frame)
            assert_equal(f.tell(), self.seek_frame)

    def test_written_zero_based_frames(self):
        # ensure that the first written DCD frame is 0
        expected_frame = 0
        with DCDFile(self.testfile) as f:
            assert_equal(f.tell(), expected_frame)

    def test_written_remarks(self):
        # ensure that the REMARKS field *can be* preserved exactly
        # in the written DCD file
        with DCDFile(self.testfile) as f:
            assert_equal(f.remarks, self.expected_remarks)

    def test_written_nsavc(self):
        # ensure that nsavc, the timesteps between frames written
        # to file, is preserved in the written DCD file
        expected = self.dcdfile_r.nsavc
        actual = DCDFile(self.testfile).nsavc
        assert_equal(actual, expected)

    def test_written_istart(self):
        # ensure that istart, the starting timestep, is preserved
        # in the written DCD file
        expected = self.dcdfile_r.istart
        actual = DCDFile(self.testfile).istart
        assert_equal(actual, expected)

    def test_written_delta(self):
        # ensure that delta, the trajectory timestep, is preserved in
        # the written DCD file
        expected = self.dcdfile_r.delta
        actual = DCDFile(self.testfile).delta
        assert_equal(actual, expected)

class DCDByteArithmeticTest(TestCase):

    def setUp(self):
        self.dcdfile = DCDFile(DCD, 'r')
        self.filesize = os.path.getsize(DCD)

    def test_relative_frame_sizes(self):
        # the first frame of a DCD file should always be >= in size
        # to subsequent frames, as the first frame contains the same
        # atoms + (optional) fixed atoms
        first_frame_size = self.dcdfile.firstframesize
        general_frame_size = self.dcdfile.framesize
        self.assertGreaterEqual(first_frame_size, general_frame_size)

    def test_file_size_breakdown(self):
        # the size of a DCD file is equivalent to the sum of the header
        # size, first frame size, and (N - 1 frames) * size per general
        # frame
        expected = self.filesize
        actual = self.dcdfile.header_size + self.dcdfile.firstframesize + \
                 ((self.dcdfile.n_frames - 1) * self.dcdfile.framesize)
        assert_equal(actual, expected)

    def test_nframessize_int(self):
        # require that the (nframessize / framesize) value used by DCDFile 
        # is an integer (because nframessize / framesize + 1 = total frames,
        # which must also be an int)
        nframessize = self.filesize - self.dcdfile.header_size - \
                           self.dcdfile.firstframesize
        self.assertTrue(float(nframessize) % float(self.dcdfile.framesize) == 0)


class DCDByteArithmeticTestNAMD(DCDByteArithmeticTest, TestCase):
    # repeat byte arithmetic tests for NAMD format DCD

    def setUp(self):
        self.dcdfile = DCDFile(DCD_NAMD_TRICLINIC, 'r')
        self.filesize = os.path.getsize(DCD_NAMD_TRICLINIC)
    

class DCDWriteTestNAMD(DCDWriteTest, TestCase):
    # repeat writing tests for NAMD format DCD

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.dcdfile = DCDFile(self.testfile, 'w')
        self.dcdfile_r = DCDFile(DCD_NAMD_TRICLINIC, 'r')
        self.natoms = 5545
        self.expected_frames = 1
        self.seek_frame = 0
        self.expected_remarks = '''Created by DCD pluginREMARKS Created 06 July, 2014 at 17:29Y5~CORD,'''

        with self.dcdfile_r as f_in, self.dcdfile as f_out:
                for frame in f_in:
                    frame = frame._asdict()
                    f_out.write(xyz=frame['x'],
                                box=frame['unitcell'].astype(np.float64),
                                step=f_in.istart,
                                natoms=frame['x'].shape[0],
                                charmm=0,
                                time_step=f_in.delta,
                                ts_between_saves=f_in.nsavc,
                                remarks=f_in.remarks)

class DCDWriteHeaderTestNAMD(DCDWriteHeaderTest, TestCase):
    # repeat header writing tests for NAMD format DCD

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.dcdfile = DCDFile(self.testfile, 'w')
        self.dcdfile_r = DCDFile(DCD_NAMD_TRICLINIC, 'r')

class DCDReadFrameTestNAMD(DCDReadFrameTest, TestCase):
    # repeat frame reading tests for NAMD format DCD

    def setUp(self):
        self.dcdfile = DCDFile(DCD_NAMD_TRICLINIC)
        self.natoms = 5545
        self.traj_length = 1
        self.new_frame = 0
        self.context_frame = 0
        self.num_iters = 0
        self.expected_remarks = 'Created by DCD pluginREMARKS Created 06 July, 2014 at 17:29Y5~CORD,'
