from __future__ import print_function
from __future__ import absolute_import

from nose.tools import raises
from numpy.testing import assert_equal
from numpy.testing import assert_allclose

from MDAnalysis.lib.formats.libdcd import DCDFile
from MDAnalysisTests.datafiles import (PSF, DCD, DCD_NAMD_TRICLINIC,
                                       PSF_NAMD_TRICLINIC,
                                       legacy_DCD_ADK_coords,
                                       legacy_DCD_NAMD_coords,
                                       legacy_DCD_c36_coords,
                                       PSF_TRICLINIC,
                                       DCD_TRICLINIC)

from MDAnalysisTests.tempdir import run_in_tempdir
from MDAnalysisTests import tempdir
import numpy as np
import os
import math
from hypothesis import given
import hypothesis.strategies as st


class TestDCDReadFrame():

    def setUp(self):
        self.dcdfile = DCD
        self.natoms = 3341
        self.traj_length = 98
        self.new_frame = 91
        self.context_frame = 22
        self.num_iters = 3
        self.selected_legacy_frames = [5, 29]
        self.legacy_data = legacy_DCD_ADK_coords
        self.expected_remarks = '''* DIMS ADK SEQUENCE FOR PORE PROGRAM                                            * WRITTEN BY LIZ DENNING (6.2008)                                               *  DATE:     6/ 6/ 8     17:23:56      CREATED BY USER: denniej0                '''
        self.expected_unit_cell = np.array([  0.,   0.,   0.,  90.,  90.,  90.],
                            dtype=np.float32)

    def test_header_remarks(self):
        # confirm correct header remarks section reading
        with DCDFile(self.dcdfile) as f:
            assert_equal(len(f.remarks), len(self.expected_remarks))

    def test_read_coords(self):
        # confirm shape of coordinate data against result from previous
        # MDAnalysis implementation of DCD file handling
        with DCDFile(self.dcdfile) as dcd:
            dcd_frame = dcd.read()
        xyz = dcd_frame[0]
        assert_equal(xyz.shape, (self.natoms, 3))

    def test_read_coord_values(self):
        # test the actual values of coordinates read in versus
        # stored values read in by the legacy DCD handling framework

        # to reduce repo storage burden, we only compare for a few
        # randomly selected frames
        legacy_DCD_frame_data = np.load(self.legacy_data)

        with DCDFile(self.dcdfile) as dcd:
            for index, frame_num in enumerate(self.selected_legacy_frames):
                dcd.seek(frame_num)
                actual_coords = dcd.read()[0]
                desired_coords = legacy_DCD_frame_data[index]
                assert_equal(actual_coords,
                             desired_coords)


    def test_read_unit_cell(self):
        # confirm unit cell read against result from previous
        # MDAnalysis implementation of DCD file handling
        with DCDFile(self.dcdfile) as dcd:
            dcd_frame = dcd.read()
        unitcell = dcd_frame[1]
        assert_allclose(unitcell, self.expected_unit_cell, rtol=1e-05)

    @raises(IOError)
    def test_seek_over_max(self):
        # should raise IOError if beyond 98th frame
        with DCDFile(DCD) as dcd:
            dcd.seek(102)

    def test_seek_normal(self):
        # frame seek within range is tested
        with DCDFile(self.dcdfile) as dcd:
            dcd.seek(self.new_frame)
            assert_equal(dcd.tell(), self.new_frame)

    @raises(IOError)
    def test_seek_negative(self):
        # frame seek with negative number
        with DCDFile(self.dcdfile) as dcd:
            dcd.seek(-78)

    def test_iteration(self):
        with DCDFile(self.dcdfile) as dcd:
            for i in range(self.num_iters):
                dcd.__next__()
            assert_equal(dcd.tell(), self.num_iters)

    def test_zero_based_frames(self):
        expected_frame = 0
        with DCDFile(self.dcdfile) as dcd:
            assert_equal(dcd.tell(), expected_frame)

    def test_length_traj(self):
        expected = self.traj_length
        with DCDFile(self.dcdfile) as dcd:
            assert_equal(len(dcd), expected)

    @raises(IOError)
    def test_open_wrong_mode(self):
        DCDFile('foo', 'e')

    @raises(IOError)
    def test_raise_not_existing(self):
        DCDFile('foo')

    def test_n_atoms(self):
        with DCDFile(self.dcdfile) as dcd:
            assert_equal(dcd.n_atoms, self.natoms)

    @raises(IOError)
    @run_in_tempdir()
    def test_read_write_mode_file(self):
        with DCDFile('foo', 'w') as f:
            f.read()

    @raises(IOError)
    def test_read_closed(self):
        with DCDFile(self.dcdfile) as dcd:
            dcd.close()
            dcd.read()

    def test_iteration_2(self):
        with DCDFile(self.dcdfile) as dcd:
            with dcd as f:
                for frame in f:
                    pass
                # second iteration should work from start again
                for frame in f:
                    pass


class TestDCDWriteHeader():

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.dcdfile = DCD

    def tearDown(self):
        try:
            os.unlink(self.testfile)
        except OSError:
            pass
        del self.tmpdir

    def test_write_header_crude(self):
        # test that _write_header() can produce a very crude
        # header for a new / empty file
        with DCDFile(self.testfile, 'w') as dcd:
            dcd._write_header(remarks='Crazy!', n_atoms=22,
                              starting_step=12, ts_between_saves=10,
                              time_step=0.02,
                              charmm=1)

        # we're not actually asserting anything, yet
        # run with: nosetests test_libdcd.py --nocapture
        # to see printed output from nose
        with open(self.testfile, "rb") as f:
            for element in f:
                print(element)

    @raises(IOError)
    def test_write_header_mode_sensitivy(self):
        # an exception should be raised on any attempt to use
        # _write_header with a DCDFile object in 'r' mode
        with DCDFile(self.dcdfile) as dcd:
            dcd._write_header(remarks='Crazy!', n_atoms=22,
                              starting_step=12, ts_between_saves=10,
                              time_step=0.02,
                              charmm=1)



class TestDCDWrite():

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.testfile2 = self.tmpdir.name + '/test2.dcd'
        self.readfile = DCD
        self.natoms = 3341
        self.expected_frames = 98
        self.seek_frame = 91
        self.expected_remarks = '''* DIMS ADK SEQUENCE FOR PORE PROGRAM                                            * WRITTEN BY LIZ DENNING (6.2008)                                               *  DATE:     6/ 6/ 8     17:23:56      CREATED BY USER: denniej0                '''
        self._write_files(testfile=self.testfile,
                          remarks_setting='input')

    def _write_files(self, testfile, remarks_setting):

        with DCDFile(self.readfile) as f_in, DCDFile(self.testfile, 'w') as f_out:
            for frame in f_in:
                if remarks_setting == 'input':
                    remarks = f_in.remarks
                else: # accept the random remarks strings from hypothesis
                    remarks = remarks_setting
                box=frame.unitcell.astype(np.float64)
                f_out.write(xyz=frame.x,
                            box=box,
                            step=f_in.istart,
                            natoms=frame.x.shape[0],
                            charmm=1, # DCD should be CHARMM
                            time_step=f_in.delta,
                            ts_between_saves=f_in.nsavc,
                            remarks=remarks)

    def tearDown(self):
        try:
            os.unlink(self.testfile)
        except OSError:
            pass
        del self.tmpdir

    @raises(IOError)
    def test_write_mode(self):
        # ensure that writing of DCD files only occurs with properly
        # opened files
        with DCDFile(self.readfile) as dcd:
            dcd.write(xyz=np.zeros((3,3)),
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
        with DCDFile(self.testfile) as test, DCDFile(self.readfile) as ref:
            curr_frame = 0
            while curr_frame < test.n_frames:
                written_unitcell = test.read()[1]
                ref_unitcell = ref.read()[1]
                curr_frame += 1
                assert_equal(written_unitcell, ref_unitcell)

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

    @given(st.text()) # handle the full unicode range of strings
    def test_written_remarks_property(self, remarks_str):
        # property based testing for writing of a wide range of string
        # values to REMARKS field
        self._write_files(testfile=self.testfile2,
                          remarks_setting=remarks_str)
        expected_remarks = remarks_str
        with DCDFile(self.testfile2) as f:
            assert_equal(f.remarks, expected_remarks)

    def test_written_nsavc(self):
        # ensure that nsavc, the timesteps between frames written
        # to file, is preserved in the written DCD file
        with DCDFile(self.readfile) as dcd_r, DCDFile(self.testfile) as dcd:
            assert_equal(dcd.nsavc, dcd_r.nsavc)

    def test_written_istart(self):
        # ensure that istart, the starting timestep, is preserved
        # in the written DCD file
        with DCDFile(self.readfile) as dcd_r, DCDFile(self.testfile) as dcd:
            assert_equal(dcd.istart, dcd_r.istart)

    def test_written_delta(self):
        # ensure that delta, the trajectory timestep, is preserved in
        # the written DCD file
        with DCDFile(self.readfile) as dcd_r, DCDFile(self.testfile) as dcd:
            assert_equal(dcd.delta, dcd_r.delta)

    def test_coord_match(self):
        # ensure that all coordinates match in each frame for the
        # written DCD file relative to original
        with DCDFile(self.testfile) as test, DCDFile(self.readfile) as ref:
            curr_frame = 0
            while curr_frame < test.n_frames:
                written_coords = test.read()[0]
                ref_coords = ref.read()[0]
                curr_frame += 1
                assert_equal(written_coords, ref_coords)


class TestDCDByteArithmetic():

    def setUp(self):

        self.dcdfile = DCD
        self._filesize = os.path.getsize(DCD)

    def test_relative_frame_sizes(self):
        # the first frame of a DCD file should always be >= in size
        # to subsequent frames, as the first frame contains the same
        # atoms + (optional) fixed atoms
        with DCDFile(self.dcdfile) as dcd:
            first_frame_size = dcd._firstframesize
            general_frame_size = dcd._framesize

            # for frame in test:
            #     written_coords = test.read()[0]
            #     ref_coords = ref.read()[0]
            #     assert_equal(written_coords, ref_coords)


class TestDCDByteArithmetic():

    def setUp(self):
        self.dcdfile = DCD
        self._filesize = os.path.getsize(DCD)

    def test_relative_frame_sizes(self):
        # the first frame of a DCD file should always be >= in size
        # to subsequent frames, as the first frame contains the same
        # atoms + (optional) fixed atoms
        with DCDFile(self.dcdfile) as dcd:
            first_frame_size = dcd._firstframesize
            general_frame_size = dcd._framesize

        assert_equal(first_frame_size >= general_frame_size, True)

    def test_file_size_breakdown(self):
        # the size of a DCD file is equivalent to the sum of the header
        # size, first frame size, and (N - 1 frames) * size per general
        # frame
        expected = self._filesize
        with DCDFile(self.dcdfile) as dcd:
            actual = dcd._header_size + dcd._firstframesize + ((dcd.n_frames - 1) * dcd._framesize)
        assert_equal(actual, expected)

    def test_nframessize_int(self):
        # require that the (nframessize / framesize) value used by DCDFile
        # is an integer (because nframessize / framesize + 1 = total frames,
        # which must also be an int)
        with DCDFile(self.dcdfile) as dcd:
            nframessize = self._filesize - dcd._header_size - dcd._firstframesize
            assert_equal(float(nframessize) % float(dcd._framesize), 0)



class TestDCDByteArithmeticNAMD(TestDCDByteArithmetic):
    # repeat byte arithmetic tests for NAMD format DCD

    def setUp(self):
        self.dcdfile = DCD_NAMD_TRICLINIC
        self._filesize = os.path.getsize(DCD_NAMD_TRICLINIC)


class TestDCDByteArithmeticCharmm36(TestDCDByteArithmetic):
    # repeat byte arithmetic tests for Charmm36 format DCD

    def setUp(self):
        self.dcdfile = DCD_TRICLINIC
        self._filesize = os.path.getsize(DCD_TRICLINIC)


class  TestDCDWriteNAMD(TestDCDWrite):
    # repeat writing tests for NAMD format DCD

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.testfile2 = self.tmpdir.name + '/test2.dcd'
        self.readfile = DCD_NAMD_TRICLINIC
        self.natoms = 5545
        self.expected_frames = 1
        self.seek_frame = 0
        self.expected_remarks = '''Created by DCD pluginREMARKS Created 06 July, 2014 at 17:29Y5~CORD,'''
        self._write_files(testfile=self.testfile,
                          remarks_setting='input')

    def test_written_unit_cell(self):
        # there's no expectation that we can write unit cell
        # data in NAMD format at the moment
        pass


class TestDCDWriteCharmm36(TestDCDWrite):
    # repeat writing tests for Charmm36 format DCD
    # no expectation that we can write unit cell info though (yet)

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.testfile2 = self.tmpdir.name + '/test2.dcd'
        self.readfile = DCD_TRICLINIC
        self.natoms = 375
        self.expected_frames = 10
        self.seek_frame = 7
        self.expected_remarks = '* CHARMM TRICLINIC BOX TESTING                                                  * (OLIVER BECKSTEIN 2014)                                                       * BASED ON NPTDYN.INP : SCOTT FELLER, NIH, 7/15/95                              '
        self._write_files(testfile=self.testfile,
                          remarks_setting='input')

    def test_written_unit_cell(self):
        # there's no expectation that we can write unit cell
        # data in NAMD format at the moment
        pass


class TestDCDWriteHeaderNAMD(TestDCDWriteHeader):
    # repeat header writing tests for NAMD format DCD

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.dcdfile = DCD_NAMD_TRICLINIC


class TestDCDWriteHeaderCharmm36(TestDCDWriteHeader):
    # repeat header writing tests for Charmm36 format DCD

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.dcdfile = DCD_TRICLINIC


class TestDCDReadFrameTestNAMD(TestDCDReadFrame):
    # repeat frame reading tests for NAMD format DCD

    def setUp(self):
        self.dcdfile = DCD_NAMD_TRICLINIC
        self.natoms = 5545
        self.traj_length = 1
        self.new_frame = 0
        self.context_frame = 0
        self.num_iters = 0
        self.selected_legacy_frames = [0]
        self.legacy_data = legacy_DCD_NAMD_coords
        self.expected_remarks = 'Created by DCD pluginREMARKS Created 06 July, 2014 at 17:29Y5~CORD,'
        # expected unit cell based on previous DCD framework read in:
        self.expected_unit_cell = np.array([ 38.42659378,  38.39310074, 44.75979996,
                                             90.        ,  90.        , 60.02891541],
                                             dtype=np.float32)


class TestDCDReadFrameTestCharmm36(TestDCDReadFrame):
    # repeat frame reading tests for Charmm36 format DCD

    def setUp(self):
        self.dcdfile = DCD_TRICLINIC
        self.natoms = 375
        self.traj_length = 10
        self.new_frame = 2
        self.context_frame = 5
        self.num_iters = 7
        self.selected_legacy_frames = [1, 4]
        self.legacy_data = legacy_DCD_c36_coords
        self.expected_remarks = '* CHARMM TRICLINIC BOX TESTING                                                  * (OLIVER BECKSTEIN 2014)                                                       * BASED ON NPTDYN.INP : SCOTT FELLER, NIH, 7/15/95                              * TEST EXTENDED SYSTEM CONSTANT PRESSURE AND TEMPERATURE                        * DYNAMICS WITH WATER BOX.                                                      *  DATE:     7/ 7/14     13:59:46      CREATED BY USER: oliver                  '
        # expected unit cell based on previous DCD framework read in:
        self.expected_unit_cell = np.array([ 35.44603729,  35.06156158,  34.15850067,
                                             91.32801819,  61.73519516, 44.4070282],
                                             dtype=np.float32)


class TestDCDWriteRandom():
    # should only be supported for Charmm24 format writing (for now)

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.testfile = self.tmpdir.name + '/test.dcd'
        self.readfile = DCD
        self.natoms = 3341
        self.expected_frames = 98
        self.seek_frame = 91
        self.expected_remarks = '''* DIMS ADK SEQUENCE FOR PORE PROGRAM                                            * WRITTEN BY LIZ DENNING (6.2008)                                               *  DATE:     6/ 6/ 8     17:23:56      CREATED BY USER: denniej0                '''

        np.random.seed(1178083)
        self.random_unitcells = np.random.random((self.expected_frames, 6)).astype(np.float64)
        with DCDFile(self.readfile) as f_in, DCDFile(self.testfile, 'w') as f_out:
            for index, frame in enumerate(f_in):
                frame_dict = frame._asdict()
                box=frame_dict['unitcell'].astype(np.float64)
                f_out.write(xyz=frame_dict['x'],
                            box=self.random_unitcells[index],
                            step=f_in.istart,
                            natoms=frame_dict['x'].shape[0],
                            charmm=1, # DCD should be CHARMM
                            time_step=f_in.delta,
                            ts_between_saves=f_in.nsavc,
                            remarks=f_in.remarks)

    def tearDown(self):
        try: 
            os.unlink(self.testfile)
        except OSError:
            pass
        del self.tmpdir

    def test_written_unit_cell_random(self):
        # written unit cell dimensions should match for all frames
        # using randomly generated unit cells but some processing
        # of the cosine data stored in charmm format is needed
        # as well as shuffling of the orders in the unitcell
        # array based on the prcoessing performed by 
        # DCDFile read and more generally relating to Issue 187
        with DCDFile(self.testfile) as test:
            curr_frame = 0
            while curr_frame < test.n_frames:
                written_unitcell = test.read()[1]
                ref_unitcell = self.random_unitcells[curr_frame]
                ref_unitcell[1] = math.degrees(math.acos(ref_unitcell[1]))
                ref_unitcell[3] = math.degrees(math.acos(ref_unitcell[3]))
                ref_unitcell[4] = math.degrees(math.acos(ref_unitcell[4]))

                _ts_order = [0, 2, 5, 4, 3, 1]
                ref_unitcell = np.take(ref_unitcell, _ts_order)
                curr_frame += 1
                assert_allclose(written_unitcell, ref_unitcell,
                                rtol=1e-05)
