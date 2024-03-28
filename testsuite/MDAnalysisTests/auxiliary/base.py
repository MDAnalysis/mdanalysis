# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import MDAnalysis as mda
import numpy as np
import pytest
from numpy.testing import assert_almost_equal, assert_equal


def test_get_bad_auxreader_format_raises_ValueError():
    # should raise a ValueError when no AuxReaders with match the specified format
    with pytest.raises(ValueError):
        mda.auxiliary.core.get_auxreader_for(format='bad-format')


class BaseAuxReference(object):
    ## assumes the reference auxiliary data has 5 steps, with three values 
    ## for each step: i, 2*i and 2^i, where i is the step number.
    ## If a particular AuxReader is such that auxiliary data is read in a 
    ## format other than np.array([i, 2*i, 2**i]), format_data() should be 
    ## overwritten tp return the appropriate format

    def __init__(self):
        self.n_steps = 5
        self.dt = 1
        self.initial_time = 0
        self.name = 'test'

        # reference description of the (basic) auxiliary reader. Will
        # have to add 'format' and 'auxdata' when creating the reference
        # for each particular reader
        self.description= {'dt':self.dt, 'represent_ts_as':'closest', 
                           'initial_time':self.initial_time, 'time_selector':None,
                           'data_selector':None, 'constant_dt':True, 
                           'cutoff': None, 'auxname': self.name}

        def reference_auxstep(i):
            # create a reference AuxStep for step i
            auxstep = mda.auxiliary.base.AuxStep(dt=self.dt,
                                             initial_time=self.initial_time)
            auxstep.step = i
            auxstep._data = self.format_data([i, 2*i, 2**i])
            return auxstep

        self.auxsteps = [reference_auxstep(i) for i in range(self.n_steps)]

        ## testing __getitem__ with slice and list. Should allow us to iterate
        ## through the specified auxiliary steps...
        self.iter_list = [1, -2]
        self.iter_list_auxsteps = [self.auxsteps[1], self.auxsteps[3]]
        self.iter_slice = slice(None, None, 2) # every second step
        self.iter_slice_auxsteps = [self.auxsteps[0], self.auxsteps[2], 
                                    self.auxsteps[4]]

        def reference_timestep(dt=1, offset=0):
            # return a trajectory timestep with specified dt, offset + move to 
            # frame 1; for use in auxiliary reading of different timesteps 
            ts = mda.coordinates.base.Timestep(0, dt=dt,
                                               time_offset=offset)
            ts.frame = 1
            return ts

        ## test reading a timestep with lower frequency. Auxiliary steps with
        ## times between [1ps, 3ps) will be assigned to this timestep, i.e.
        ## step 1 (1 ps) and step 2 (2 ps).
        self.lower_freq_ts = reference_timestep(dt=2, offset=0)
        # 'closest' representative value will match step 2
        self.lowf_closest_rep = self.format_data([2, 2*2, 2**2])
        # 'average' representative value
        self.lowf_average_rep = self.format_data([1.5, 3, 3])

        ## test reading a timestep with higher frequency. Auxiliart steps with 
        ## times between [0.25ps, 0.75ps) will be assigned to this timestep, i.e.
        ## no auxiliary steps
        self.higher_freq_ts = reference_timestep(dt=0.5, offset=0)
        self.highf_rep = self.format_data([np.nan, np.nan, np.nan])

        ## test reading a timestep that is offset. Auxiliary steps with
        ## times between [0.75ps, 1.75ps) will be assigned to this timestep, i.e.
        ## step 1 (1 ps)
        self.offset_ts = reference_timestep(dt=1, offset=0.25)
        # 'closest' representative value will match step 1 data
        self.offset_closest_rep = self.format_data([1, 2*1, 2**1])

        ## testing cutoff for representative values
        self.cutoff = 0
        # for 'average': use low frequenct timestep, only step 2 within 0ps cutoff
        self.lowf_cutoff_average_rep = self.format_data([2, 2*2, 2**2])
        # for 'closest': use offset timestep; no timestep within 0ps cutoff
        self.offset_cutoff_closest_rep = self.format_data([np.nan, np.nan, np.nan])

        ## testing selection of time/data. Overload for each auxilairy format 
        ## as appropraite.
        # default None behavior set here so won't get errors when time/data 
        # selection not implemented. 
        self.time_selector = None 
        self.select_time_ref = np.arange(self.n_steps)
        self.data_selector = None 
        self.select_data_ref = [self.format_data([2*i, 2**i]) for i in range(self.n_steps)]


    def format_data(self, data):
        ## overload if auxiliary reader will read data with a format different  
        ## to e.g. np.array([0, 0, 1])
        return np.array(data)


class BaseAuxReaderTest(object):

    def test_raise_error_no_auxdata_provided(self, ref, ref_universe):
        with pytest.raises(ValueError, match="No input `auxdata`"):
            ref_universe.trajectory.add_auxiliary()

    def test_n_steps(self, ref, reader):
        assert len(reader) == ref.n_steps, "number of steps does not match"

    def test_dt(self, ref, reader):
        assert reader.dt == ref.dt, "dt does not match"

    def test_initial_time(self, ref, reader):
        assert reader.initial_time == ref.initial_time, "initial time does not match"

    def test_first_step(self, ref, reader):
        # on first loading we should start at step 0
        assert_auxstep_equal(reader.auxstep, ref.auxsteps[0])

    def test_next(self, ref, reader):
        # should take us to step 1
        next(reader)
        assert_auxstep_equal(reader.auxstep, ref.auxsteps[1])

    def test_rewind(self, ref, reader):
        # move to step 1...
        reader.next()
        # now rewind should read step 0
        reader.rewind()
        assert_auxstep_equal(reader.auxstep, ref.auxsteps[0])

    def test_move_to_step(self, ref, reader):
        # should take us to step 3
        reader[3]
        assert_auxstep_equal(reader.auxstep, ref.auxsteps[3])

    def test_last_step(self, ref, reader):
        # should take us to the last step
        reader[-1]
        assert_auxstep_equal(reader.auxstep, ref.auxsteps[-1])

    def test_next_past_last_step_raises_StopIteration(self, ref, reader):
        # should take us to the last step
        reader[-1]
        # if we try to move to next step from here, should raise StopIteration
        with pytest.raises(StopIteration):
            reader.next()

    def test_move_to_invalid_step_raises_IndexError(self, ref, reader):
        # last step is number n_steps -1 ; if we try move to step number
        # n_steps we should get a ValueError
        with pytest.raises(IndexError):
            reader[ref.n_steps]

    def test_invalid_step_to_time_raises_ValueError(self, reader):
        # last step is number n_steps-1; if we try to run step_to_time on
        # step n_steps we should get a ValueError
        with pytest.raises(ValueError):
            reader.step_to_time(reader.n_steps)

    def test_iter(self,ref, reader):
        for i, val in enumerate(reader):
            assert_auxstep_equal(val, ref.auxsteps[i])

    def test_iter_list(self, ref, reader):
        # test using __getitem__ with a list
        for i, val in enumerate(reader[ref.iter_list]):
            assert_auxstep_equal(val, ref.iter_list_auxsteps[i])

    def test_iter_slice(self, ref, reader):
        # test using __getitem__ with a slice
        for i, val in enumerate(reader[ref.iter_slice]):
            assert_auxstep_equal(val, ref.iter_slice_auxsteps[i])

    def test_slice_start_after_stop_raises_IndexError(self, reader):
        #should raise IndexError if start frame after end frame
        with pytest.raises(IndexError):
            reader[2:1]

    def test_slice_out_of_range_raises_IndexError(self, ref, reader):
        # should raise IndexError if indices our of range
        with pytest.raises(IndexError):
            reader[ref.n_steps:]

    def test_slice_non_int_raises_TypeError(self, reader):
        # should raise TypeError if try pass in non-integer to slice
        with pytest.raises(TypeError):
            reader['a':]

    def test_bad_represent_raises_ValueError(self, reader):
        # if we try to set represent_ts_as to something not listed as a
        # valid option, we should get a ValueError
        with pytest.raises(ValueError):
            reader.represent_ts_as = 'invalid-option'

    def test_time_selector(self, ref):
        # reload the reader, passing a time selector
        reader = ref.reader(ref.testdata,
                                      time_selector = ref.time_selector)
        # time should still match reference time for each step
        for i, val in enumerate(reader):
            assert val.time == ref.select_time_ref[i], "time for step {} does not match".format(i)

    def test_time_non_constant_dt(self, reader):
        reader.constant_dt = False
        with pytest.raises(ValueError, match="If dt is not constant, must have a valid time selector"):
            reader.time

    def test_time_selector_manual(self, ref):
        reader = ref.reader(ref.testdata,
                            time_selector = ref.time_selector)
        # Manually set time selector
        reader.time_selector = ref.time_selector
        for i, val in enumerate(reader):
            assert val.time == ref.select_time_ref[i], "time for step {} does not match".format(i)

    def test_data_selector(self, ref):
        # reload reader, passing in a data selector
        reader = ref.reader(ref.testdata,
                                      data_selector=ref.data_selector)
        # data should match reference data for each step
        for i, val in enumerate(reader):
            assert_equal(val.data, ref.select_data_ref[i], "data for step {0} does not match".format(i))

    def test_no_constant_dt(self, ref):
        ## assume we can select time...
        # reload reader, without assuming constant dt
        reader = ref.reader(ref.testdata,
                                      time_selector=ref.time_selector,
                                      constant_dt=False)
        # time should match reference for selecting time, for each step
        for i, val in enumerate(reader):
            assert val.time == ref.select_time_ref[i], "data for step {} does not match".format(i)

    def test_update_ts_without_auxname_raises_ValueError(self, ref):
        # reload reader without auxname
        with pytest.raises(ValueError):
            reader = ref.reader(ref.testdata)
            ts = ref.lower_freq_ts
            reader.update_ts(ts)

    def test_read_lower_freq_timestep(self, ref, reader):
        # test reading a timestep with lower frequency
        ts = ref.lower_freq_ts
        reader.update_ts(ts)
        # check the value set in ts is as we expect
        assert_almost_equal(ts.aux.test, ref.lowf_closest_rep,
                            err_msg="Representative value in ts.aux does not match")

    def test_represent_as_average(self, ref, reader):
        # test the 'average' option for 'represent_ts_as'
        # reset the represent method to 'average'...
        reader.represent_ts_as = 'average'
        # read timestep; use the low freq timestep
        ts = ref.lower_freq_ts
        reader.update_ts(ts)
        # check the representative value set in ts is as expected
        assert_almost_equal(ts.aux.test, ref.lowf_average_rep,
                            err_msg="Representative value does not match when "
                                    "using with option 'average'")

    def test_represent_as_average_with_cutoff(self, ref, reader):
        # test the 'represent_ts_as' 'average' option when we have a cutoff set
        # set the cutoff...
        reader.cutoff = ref.cutoff
        # read timestep; use the low frequency timestep
        ts = ref.lower_freq_ts
        reader.update_ts(ts)
        # check representative value set in ts is as expected
        assert_almost_equal(ts.aux.test, ref.lowf_cutoff_average_rep,
                            err_msg="Representative value does not match when "
                                    "applying cutoff")

    def test_read_offset_timestep(self, ref, reader):
        # try reading a timestep offset from auxiliary
        ts = ref.offset_ts
        reader.update_ts(ts)
        assert_almost_equal(ts.aux.test, ref.offset_closest_rep,
                            err_msg="Representative value in ts.aux does not match")

    def test_represent_as_closest_with_cutoff(self, ref, reader):
        # test the 'represent_ts_as' 'closest' option when we have a cutoff set
        # set the cutoff...
        reader.cutoff = ref.cutoff
        # read timestep; use the offset timestep
        ts = ref.offset_ts
        reader.update_ts(ts)
        # check representative value set in ts is as expected
        assert_almost_equal(ts.aux.test, ref.offset_cutoff_closest_rep,
                            err_msg="Representative value does not match when "
                                    "applying cutoff")

    def test_read_higher_freq_timestep(self, ref, reader):
        # try reading a timestep with higher frequency
        ts = ref.higher_freq_ts
        reader.update_ts(ts)
        assert_almost_equal(ts.aux.test, ref.highf_rep,
                            err_msg="Representative value in ts.aux does not match")

    def test_get_auxreader_for(self, ref, reader):
        # check guesser gives us right reader
        reader = mda.auxiliary.core.get_auxreader_for(ref.testdata)
        assert reader == ref.reader

    def test_iterate_through_trajectory(self, ref, ref_universe):
        # check the representative values of aux for each frame are as expected
        # trajectory here has same dt, offset; so there's a direct correspondence
        # between frames and steps
        for i, ts in enumerate(ref_universe.trajectory):
            assert_equal(ts.aux.test, ref.auxsteps[i].data,
                         "representative value does not match when "
                         "iterating through all trajectory timesteps")

    def test_iterate_as_auxiliary_from_trajectory(self, ref, ref_universe):
        # check representative values of aux for each frame are as expected
        # trajectory here has same dt, offset, so there's a direct correspondence
        # between frames and steps, and iter_as_aux will run through all frames
        for i, ts in enumerate(ref_universe.trajectory.iter_as_aux('test')):
            assert_equal(ts.aux.test, ref.auxsteps[i].data,
                         "representative value does not match when "
                         "iterating through all trajectory timesteps")

    def test_auxiliary_read_ts_rewind(self, ref_universe):
        # AuxiliaryBase.read_ts() should retrieve the correct step after
        # reading the last one. Issue #2674 describes a case in which the
        # object gets stuck on the last frame.
        aux_info_0 = ref_universe.trajectory[0].aux.test
        ref_universe.trajectory[-1]
        aux_info_0_rewind = ref_universe.trajectory[0].aux.test
        assert_equal(aux_info_0, aux_info_0_rewind,
                     "aux info was retrieved incorrectly "
                     "after reading the last step")

    def test_get_description(self, ref, reader):
        description = reader.get_description()
        for attr in ref.description:
            assert description[attr] == ref.description[attr], "'Description' does not match for {}".format(attr)

    def test_load_from_description(self, reader):
        description = reader.get_description()
        new = mda.auxiliary.core.auxreader(**description)
        assert new == reader, "AuxReader reloaded from description does not match"

    def test_step_to_frame_out_of_bounds(self, reader, ref):

        ts = mda.coordinates.base.Timestep(0, dt=ref.dt)

        assert reader.step_to_frame(-1, ts) is None
        assert reader.step_to_frame(reader.n_steps, ts) is None

    def test_step_to_frame_no_time_diff(self, reader, ref):

        ts = mda.coordinates.base.Timestep(0, dt=ref.dt)

        for idx in range(reader.n_steps):

            assert reader.step_to_frame(idx, ts) == idx

    def test_step_to_frame_time_diff(self, reader, ref):

        # Timestep is 0.1 longer than auxiliary data
        ts = mda.coordinates.base.Timestep(0, dt=ref.dt + 0.1)

        # Test all 5 frames
        for idx in range(5):

            frame, time_diff = reader.step_to_frame(idx, ts, return_time_diff=True)

            assert frame == idx
            np.testing.assert_almost_equal(time_diff, idx * 0.1)

    def test_go_to_step_fail(self, reader):

        with pytest.raises(ValueError, match="Step index [0-9]* is not valid for auxiliary"):
            reader._go_to_step(reader.n_steps)

    @pytest.mark.parametrize("constant", [True, False])
    def test_set_constant_dt(self, reader, constant):

        reader.constant_dt = constant

        assert reader.constant_dt == constant

    def test_copy(self, reader):
        new_reader = reader.copy()
        assert reader == new_reader
        assert reader is not new_reader


def assert_auxstep_equal(A, B):
    if not isinstance(A, mda.auxiliary.base.AuxStep):
        raise AssertionError('A is not of type AuxStep')
    if not isinstance(B, mda.auxiliary.base.AuxStep):
        raise AssertionError('B is not of type AuxStep')
    if A.step != B.step:
        raise AssertionError('A and B refer to different steps: A.step = {}, '
                             'B.step = {}'.format(A.step, B.step))
    if A.time != B.time:
        raise AssertionError('A and B have different times: A.time = {}, '
                             'B.time = {}'.format(A.time, B.time))
    if isinstance(A.data, dict):
        for term in A.data:
            assert_almost_equal(A.data[term], B.data[term])
    else:
        if any(A.data != B.data):
            # e.g. XVGReader
            raise AssertionError('A and B have different data: A.data = {}, '
                                 'B.data = {}'.format(A.data, B.data))
