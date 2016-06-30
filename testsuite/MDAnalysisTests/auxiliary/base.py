import numpy as np
from numpy.testing import (assert_equal, assert_raises, assert_almost_equal,
                           raises)

import MDAnalysis as mda

from MDAnalysisTests.datafiles import (COORDINATES_XTC, COORDINATES_TOPOLOGY)

@raises(ValueError)
def test_get_bad_auxreader_format_raises_ValueError():
    mda.auxiliary.core.get_auxreader_for('data', format='bad-format')

class BaseAuxReference(object):
    def __init__(self):
        self.all_data = [[0, 0, 1], [1, 2, 2], [2, 4, 4], [3, 6, 8], [4, 8, 16]]

        # make time first column
        self.time_col = 0
        self.all_times = [i[0] for i in self.all_data]
        self.all_step_data = [[i[1], i[2]] for i in self.all_data]

        # for testing data col selection
        self.set_data_cols = [1]
        self.set_data_cols_vals = [[i[1]] for i in self.all_data]

        self.n_steps = len(self.all_data)
        self.n_cols = len(self.all_data[0])
        self.dt = self.all_times[1] - self.all_times[0]
        self.initial_time = self.all_times[0]

        self.cutoff = 0
 
        # reference timestep with lower frequency (higher dt)
        lowf = 2
        self.lowf = {}
        self.lowf_ts = mda.coordinates.base.Timestep(0, dt=lowf*self.dt)
        self.lowf_ts.frame = 1
        # should correspond to steps 1 and 2...
        self.lowf['data'] = [self.all_step_data[1], self.all_step_data[2]]
        self.lowf['rep'] = self.all_step_data[2]
        self.lowf['last'] = 2    # should be on this step after reading ts
        # for testing average + cutoff calc_representative option...
        self.lowf['average'] = [3, 3]
        self.lowf['cutoff_average'] = [4, 4]


        # reference timestep, offset
        offset = 0.25
        self.offset = {}
        self.offset_ts = mda.coordinates.base.Timestep(0, dt=self.dt, 
                                                  time_offset=self.dt*offset)
        self.offset_ts.frame = 1
        # should correspond to step 1
        self.offset['data'] = [self.all_step_data[1]]
        self.offset['rep'] = self.all_step_data[1]
        self.offset['last'] = 1
        # for testing cutoff...
        self.offset['cutoff_closest'] = [np.nan, np.nan]


        # reference timestep with higher frequency (lower dt)
        highf = 0.5
        self.highf = {}
        self.highf_ts = mda.coordinates.base.Timestep(0, dt=highf*self.dt)
        self.highf_ts.frame = 1
        # should fall between step 0 and 1 - no steps assigned to this ts
        self.highf['data'] = []
        self.highf['rep'] = [np.nan, np.nan]
        self.highf['last'] = 0

        self.description= {'dt':self.dt, 'represent_ts_as':'closest', 
                           'initial_time':0, 'time_col':self.time_col,
                           'data_cols':[1,2], 'constant_dt':True, 
                           'cutoff':-1,}


class BaseAuxReaderTest(object):
    def __init__(self, reference):
        self.ref = reference
        # time assumed to be in first column!
        self.reader = self.ref.reader(self.ref.testdata, time_col=0, 
                                      auxname='test')
        self.ref.description['auxname'] = 'test'

    def tearDown(self):
        del self.reader

    def test_n_steps(self):
        assert_equal(len(self.reader), self.ref.n_steps,
                     "number of steps does not match")

    def test_n_cols(self):
        assert_equal(self.reader.n_cols, self.ref.n_cols,
                     "number of columns does not match")

    def test_dt(self):
        assert_equal(self.reader.dt, self.ref.dt,
                     "dt does not match")

    def test_initial_time(self):
        assert_equal(self.reader.initial_time, self.ref.initial_time,
                     "initial time does not match")

    def check_step(self, i):
        # check step_data and time match expected values for step i
        assert_equal(self.reader.step_data, self.ref.all_step_data[i],
                     "Step data for step {0} does not match".format(i))
        assert_equal(self.reader.time, self.ref.all_times[i],
                     "Time for step {0} does not match".format(i))

    def test_first_step(self):
        # on first loading we should start at step 0
        self.check_step(0)

    def test_go_to_first_step(self):
        self.reader.next()
        # go_to_first_step should read step 0
        self.reader.go_to_first_step()
        self.check_step(0)

    def test_next_to_second_frame(self):
        # should take us to step 1
        next(self.reader)
        self.check_step(1)

    def test_go_to_step(self):
        self.reader.go_to_step(3)
        self.check_step(3)

    def test_next_past_last_step_raises_StopIteration(self):
        self.reader.go_to_step(self.reader.n_steps-1)
        assert_raises(StopIteration, self.reader.next)

    def test_iter(self):
        for i, val in enumerate(self.reader):
            # check yielded val is as expected
            ref_data = {'time': self.ref.all_times[i], 
                        'data': self.ref.all_step_data[i]}
            assert_equal(val, ref_data)
            # also check values set in reader 
            self.check_step(i)

    def test_data_cols(self):
        # reload reader, imposing a selection for which columns to
        # include in step_data
        self.reader = self.ref.reader(self.ref.testdata, 
                                      data_cols=self.ref.set_data_cols)
        for i, val in enumerate(self.reader):
            assert_equal(val['data'], self.ref.set_data_cols_vals[i],
                         "step_data for step {0} does not match".format(i))

    def test_no_time_col(self):
        # reload reader, without setting time column; pass dt and initial_time
        # to compensate
        self.reader = self.ref.reader(self.ref.testdata, dt=self.ref.dt, 
                                      initial_time=self.ref.initial_time)
        for i, val in enumerate(self.reader):
            assert_equal(val['data'], self.ref.all_data[i],
                         "step_data for step {0} does not match".format(i))

    def test_no_constant_dt(self):
        # reload reader, without assuming constant dt
        self.reader = self.ref.reader(self.ref.testdata, time_col=0,
                                      constant_dt=False)
        for i in range(self.ref.n_steps):
            assert_almost_equal(self.reader.step_to_time(i), 
                                self.ref.all_times[i],
                              err_msg="step_data for step {0} does not match".format(i))

    @raises(ValueError)
    def test_bad_represent_raises_ValueError(self):
        self.reader.represent_ts_as = 'invalid-option'

    @raises(ValueError)
    def test_time_col_out_of_range_raises_ValueError(self):
        self.reader.time_col = self.reader.n_cols 

    @raises(ValueError)
    def test_data_col_out_of_range_raises_ValueError(self):
        self.reader.data_cols = [self.reader.n_cols]

    @raises(ValueError)
    def test_go_to_invalid_step_raises_ValueError(self):
        self.reader.go_to_step(self.reader.n_steps)

    @raises(ValueError)
    def test_invalid_step_to_time_raises_ValueError(self):
        self.reader.step_to_time(self.reader.n_steps)

    def test_not_setting_auxname(self):
        # reading a timestep without setting *auxname* should mean the 
        # timesteps aux namespace is unchanged
        self.reader = self.ref.reader(self.ref.testdata)
        aux_before = self.ref.lowf_ts.copy().aux
        self.reader.read_ts(self.ref.lowf_ts)
        assert_equal(self.ref.lowf_ts.aux.__dict__, aux_before.__dict__)

    def check_timestep(self, ref):
        assert_equal(self.reader.ts_data, ref['data'],
                    "List of step_data for steps assigned to ts does not match")
        assert_almost_equal(self.reader.ts_rep, ref['rep'],
                            err_msg="Representative value for ts does not math")
        assert_equal(self.reader.step, ref['last'],
                     "Auxiliary ends on wrong step after reading ts")

    def test_read_ts_lowf(self):
        # try reading a timestep with lower frequency/higher dt
        self.reader.read_ts(self.ref.lowf_ts)
        # check AuxReader values as expected
        self.check_timestep(self.ref.lowf)
        # check aux namespace in ts updated as expected
        assert_almost_equal(self.ref.lowf_ts.aux.test, self.ref.lowf['rep'],
                            err_msg="Representative value in ts.aux does not match")

    def test_read_ts_offset(self):
        # try reading a timeset offset from auxiliary
        self.reader.read_ts(self.ref.offset_ts)
        self.check_timestep(self.ref.offset)

    def test_read_ts_highf(self):
        # try reading a timestep with higher frequency/lower dt
        self.reader.read_ts(self.ref.highf_ts)
        self.check_timestep(self.ref.highf)

    def test_rep_as_average(self):
        self.reader.represent_ts_as = 'average'
        self.reader.read_ts(self.ref.lowf_ts)
        assert_almost_equal(self.reader.ts_rep, self.ref.lowf['average'], 
                            err_msg="Representative value does not match when "
                                    "using with option 'average'")

    def test_cutoff_closest(self):
        self.reader.cutoff = self.ref.cutoff
        self.reader.read_ts(self.ref.offset_ts)
        assert_almost_equal(self.reader.ts_rep, self.ref.offset['cutoff_closest'],
                            err_msg="Representative value does not match when "
                                    "applying cutoff")

    def test_cutoff_average(self):
        self.reader.cutoff = self.ref.cutoff
        self.reader.repreesnt_ts_as = 'average'
        self.reader.read_ts(self.ref.lowf_ts)
        assert_almost_equal(self.reader.ts_rep, self.ref.lowf['cutoff_average'],
                            err_msg="Representative value does not match when "
                                    "using option 'average' and applying cutoff")

    def test_get_auxreader_for(self):
        # check guesser gives us right reader
        reader = mda.auxiliary.core.get_auxreader_for(self.ref.testdata)
        assert_equal(reader, self.ref.reader)

    def test_read_with_trajectory(self):
        # check we can load with trajectory as expected
        u = mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_XTC)
        u.trajectory.add_auxiliary('test', self.ref.testdata, time_col=0)
        iter_values = [ts.aux.test for ts in u.trajectory]
        # ref trajectory has same ts as ref auxdata so expect the same values
        assert_equal(iter_values, self.ref.all_step_data,
                     "representative value does not match when iterating through "
                     "all trajectory timesteps")
        # and so should also be the same if we iter_as_aux
        aux_iter_values = [ts.aux.test for ts in u.trajectory.iter_as_aux('test')]
        assert_equal(aux_iter_values, self.ref.all_step_data, 
                     "representative value does not match when when iterating "
                     "using iter_as_aux")

    def test_get_description(self):
        description = self.reader.get_description()
        for attr in self.ref.description:
            assert_equal(description[attr], self.ref.description[attr],
                         "'Description' does not match for {0}".format(attr))

    def test_load_from_description(self):
        description = self.reader.get_description()
        auxdata = description.pop('auxdata')
        format = description.pop('format')
        reader = mda.auxiliary.core.get_auxreader_for(format=format)
        assert_equal(reader(auxdata, **description), self.reader,
                     "AuxReader reloaded from description does not match")
