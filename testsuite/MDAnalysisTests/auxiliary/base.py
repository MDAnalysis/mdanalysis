# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
from six.moves import range
import numpy as np
from numpy.testing import (assert_equal, assert_raises, assert_almost_equal,
                           raises)

import MDAnalysis as mda

from MDAnalysisTests.datafiles import (COORDINATES_XTC, COORDINATES_TOPOLOGY)

@raises(ValueError)
def test_get_bad_auxreader_format_raises_ValueError():
    # should raise a ValueError when no AuxReaders with match the specified format
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
                           'cutoff':-1, 'auxname':self.name}

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
    def __init__(self, reference):
        self.ref = reference
        self.reader = self.ref.reader(self.ref.testdata, initial_time=self.ref.initial_time,
                                      dt=self.ref.dt, auxname=self.ref.name,
                                      time_selector=None, data_selector=None)

    def tearDown(self):
        del self.reader

    def test_n_steps(self):
        assert_equal(len(self.reader), self.ref.n_steps,
                     "number of steps does not match")

    def test_dt(self):
        assert_equal(self.reader.dt, self.ref.dt,
                     "dt does not match")

    def test_initial_time(self):
        assert_equal(self.reader.initial_time, self.ref.initial_time,
                     "initial time does not match")

    def test_first_step(self):
        # on first loading we should start at step 0
        assert_auxstep_equal(self.reader.auxstep, self.ref.auxsteps[0])

    def test_next(self):
        # should take us to step 1
        next(self.reader)
        assert_auxstep_equal(self.reader.auxstep, self.ref.auxsteps[1])

    def test_rewind(self):
        # move to step 1...
        self.reader.next()
        # now rewind should read step 0
        self.reader.rewind()
        assert_auxstep_equal(self.reader.auxstep, self.ref.auxsteps[0])

    def test_move_to_step(self):
        # should take us to step 3
        self.reader[3]
        assert_auxstep_equal(self.reader.auxstep, self.ref.auxsteps[3])

    def test_last_step(self):
        # should take us to the last step
        self.reader[-1]
        assert_auxstep_equal(self.reader.auxstep, self.ref.auxsteps[-1])

    def test_next_past_last_step_raises_StopIteration(self):
        # should take us to the last step
        self.reader[-1]
        # if we try to move to next step from here, should raise StopIteration
        assert_raises(StopIteration, self.reader.next)

    @raises(IndexError)
    def test_move_to_invalid_step_raises_IndexError(self):
        # last step is number n_steps -1 ; if we try move to step number 
        # n_steps we should get a ValueError
        self.reader[self.ref.n_steps]

    @raises(ValueError)
    def test_invalid_step_to_time_raises_ValueError(self):
        # last step is number n_steps-1; if we try to run step_to_time on 
        # step n_steps we should get a ValueError
        self.reader.step_to_time(self.reader.n_steps)

    def test_iter(self):
        for i, val in enumerate(self.reader):
            assert_auxstep_equal(val, self.ref.auxsteps[i])

    def test_iter_list(self):
        # test using __getitem__ with a list
        for i, val in enumerate(self.reader[self.ref.iter_list]):
            assert_auxstep_equal(val, self.ref.iter_list_auxsteps[i])

        
    def test_iter_slice(self):
        # test using __getitem__ with a slice
        for i, val in enumerate(self.reader[self.ref.iter_slice]):
            assert_auxstep_equal(val, self.ref.iter_slice_auxsteps[i])

    @raises(IndexError)
    def test_slice_start_after_stop_raises_IndexError(self):
        #should raise IndexError if start frame after end frame
        self.reader[2:1]

    @raises(IndexError)
    def test_slice_out_of_range_raises_IndexError(self):
        # should raise IndexError if indices our of range
        self.reader[self.ref.n_steps:]

    @raises(TypeError)
    def test_slice_non_int_raises_TypeError(self):
        # should raise TypeError if try pass in non-integer to slice
        self.reader['a':]

    @raises(ValueError)
    def test_bad_represent_raises_ValueError(self):
        # if we try to set represent_ts_as to something not listed as a 
        # valid option, we should get a ValueError
        self.reader.represent_ts_as = 'invalid-option'

    def test_time_selector(self):
        # reload the reader, passing a time selector
        self.reader = self.ref.reader(self.ref.testdata,
                                      time_selector = self.ref.time_selector)
        # time should still match reference time for each step
        for i, val in enumerate(self.reader):
            assert_equal(val.time, self.ref.select_time_ref[i], 
                         "time for step {} does not match".format(i))

    def test_data_selector(self):
        # reload reader, passing in a data selector
        self.reader = self.ref.reader(self.ref.testdata, 
                                      data_selector=self.ref.data_selector)
        # data should match reference data for each step
        for i, val in enumerate(self.reader):
            assert_equal(val.data, self.ref.select_data_ref[i],
                         "data for step {0} does not match".format(i))

    def test_no_constant_dt(self):
        ## assume we can select time...
        # reload reader, without assuming constant dt
        self.reader = self.ref.reader(self.ref.testdata, 
                                      time_selector=self.ref.time_selector,
                                      constant_dt=False)
        # time should match reference for selecting time, for each step
        for i, val in enumerate(self.reader):
            assert_equal(val.time, self.ref.select_time_ref[i], 
                         "data for step {} does not match".format(i))

    @raises(ValueError)
    def test_update_ts_without_auxname_raises_ValueError(self):
        # reload reader without auxname
        self.reader = self.ref.reader(self.ref.testdata)
        ts = self.ref.lower_freq_ts       
        self.reader.update_ts(ts)

    def test_read_lower_freq_timestep(self):
        # test reading a timestep with lower frequency
        ts = self.ref.lower_freq_ts       
        self.reader.update_ts(ts)
        # check the value set in ts is as we expect
        assert_almost_equal(ts.aux.test, self.ref.lowf_closest_rep,
                            err_msg="Representative value in ts.aux does not match")

    def test_represent_as_average(self):
        # test the 'average' option for 'represent_ts_as'
        # reset the represent method to 'average'...
        self.reader.represent_ts_as = 'average'
        # read timestep; use the low freq timestep
        ts = self.ref.lower_freq_ts
        self.reader.update_ts(ts)
        # check the representative value set in ts is as expected
        assert_almost_equal(ts.aux.test, self.ref.lowf_average_rep, 
                            err_msg="Representative value does not match when "
                                    "using with option 'average'")

    def test_represent_as_average_with_cutoff(self):
        # test the 'represent_ts_as' 'average' option when we have a cutoff set
        # set the cutoff...
        self.reader.cutoff = self.ref.cutoff
        # read timestep; use the low frequency timestep
        ts = self.ref.lower_freq_ts
        self.reader.update_ts(ts)
        # check representative value set in ts is as expected
        assert_almost_equal(ts.aux.test, self.ref.lowf_cutoff_average_rep,
                            err_msg="Representative value does not match when "
                                    "applying cutoff")

    def test_read_offset_timestep(self):
        # try reading a timestep offset from auxiliary
        ts = self.ref.offset_ts
        self.reader.update_ts(ts)
        assert_almost_equal(ts.aux.test, self.ref.offset_closest_rep,
                            err_msg="Representative value in ts.aux does not match")

    def test_represent_as_closest_with_cutoff(self):
        # test the 'represent_ts_as' 'closest' option when we have a cutoff set
        # set the cutoff...
        self.reader.cutoff = self.ref.cutoff
        # read timestep; use the offset timestep
        ts = self.ref.offset_ts
        self.reader.update_ts(ts)
        # check representative value set in ts is as expected
        assert_almost_equal(ts.aux.test, self.ref.offset_cutoff_closest_rep,
                            err_msg="Representative value does not match when "
                                    "applying cutoff")

    def test_read_higher_freq_timestep(self):
        # try reading a timestep with higher frequency
        ts = self.ref.higher_freq_ts
        self.reader.update_ts(ts)
        assert_almost_equal(ts.aux.test, self.ref.highf_rep,
                            err_msg="Representative value in ts.aux does not match")

    def test_get_auxreader_for(self):
        # check guesser gives us right reader
        reader = mda.auxiliary.core.get_auxreader_for(self.ref.testdata)
        assert_equal(reader, self.ref.reader)

    def test_iterate_through_trajectory(self):
        # add to trajectory
        u = mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_XTC)
        u.trajectory.add_auxiliary('test', self.ref.testdata)
        # check the representative values of aux for each frame are as expected
        # trajectory here has same dt, offset; so there's a direct correspondence
        # between frames and steps
        for i, ts in enumerate(u.trajectory):
            assert_equal(ts.aux.test, self.ref.auxsteps[i].data,
                     "representative value does not match when iterating through "
                     "all trajectory timesteps")
        u.trajectory.close()

    def test_iterate_as_auxiliary_from_trajectory(self):
        # add to trajectory
        u = mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_XTC)
        u.trajectory.add_auxiliary('test', self.ref.testdata)
        # check representative values of aux for each frame are as expected
        # trahectory here has same dt, offset, so there's a direct correspondence
        # between frames and steps, and iter_as_aux will run through all frames
        for i, ts in enumerate(u.trajectory.iter_as_aux('test')):
            assert_equal(ts.aux.test, self.ref.auxsteps[i].data,
                     "representative value does not match when iterating through "
                     "all trajectory timesteps")
        u.trajectory.close()

    def test_get_description(self):
        description = self.reader.get_description()
        for attr in self.ref.description:
            assert_equal(description[attr], self.ref.description[attr],
                         "'Description' does not match for {}".format(attr))

    def test_load_from_description(self):
        description = self.reader.get_description()
        new = mda.auxiliary.core.auxreader(**description)
        assert_equal(new, self.reader,
                     "AuxReader reloaded from description does not match")


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
    if all(A.data != B.data):
        raise AssertionError('A and B have different data: A.data = {}, '
                             'B.data = {}'.format(A.data, B.data))
