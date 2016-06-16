import numpy as np
from numpy.testing import (assert_equal, assert_raises, assert_almost_equal,
                           assert_array_almost_equal, raises)
import MDAnalysis as mda

from MDAnalysisTests.datafiles import AUX_XVG

class BaseAuxReference(object):
    def __init__(self):
        self.all_data = [[0,1], [1,2], [2,4], [3,8], [4,16]]

        self.all_step_data = [[i[1]] for i in self.all_data]
        self.all_times = [i[0] for i in self.all_data]
        self.n_steps = len(self.all_data)
        self.n_cols = len(self.all_data[0])
        self.dt = self.all_times[1] - self.all_times[0]
        self.initial_time = self.all_times[0]

        lowf_dt = self.dt * 2
        self.ts_lowf = mda.coordinates.base.Timestep(0, dt=lowf_dt)
        self.ts_lowf.frame = 1
        self.ts_lowf_data = [[2], [4]]
        self.ts_lowf_rep = [4] 
        self.ts_lowf_rep_average = [3]
        self.ts_lowf_last_step = 2

        highf_dt = self.dt / 2.
        self.ts_highf = mda.coordinates.base.Timestep(0, dt=highf_dt)
        self.ts_highf.frame = 1
        self.ts_highf_data = []
        self.ts_highf_rep = [np.nan]
        self.ts_highf_last_step = 0


class BaseAuxReaderTest(object):
    def __init__(self, reference):
        self.ref = reference

    def setUp(self):
        self.reader = self.ref.reader(self.ref.testdata, name='test')

    def tearDown(self):
        del self.reader

    def test_n_steps(self):
        assert_equal(self.reader.n_steps, self.ref.n_steps)

    def test_n_cols(self):
        assert_equal(self.reader.n_cols, self.ref.n_cols)

    def test_dt(self):
        assert_equal(self.reader.dt, self.ref.dt)

    def test_initial_time(self):
        assert_equal(self.reader.initial_time, self.ref.initial_time)

    def test_first_step(self):
        self.reader.go_to_first_step()
        assert_equal(self.reader.step_data, self.ref.all_step_data[0])
        assert_equal(self.reader.time, self.ref.all_times[0])

    def test_next_to_second_frame(self):
        self.reader.next()
        assert_equal(self.reader.step_data, self.ref.all_step_data[1])
        assert_equal(self.reader.time, self.ref.all_times[1])

    def test_go_to_step(self):
        ## (not all aux readers might have go_to_step)
        self.reader.go_to_step(3)
        assert_equal(self.reader.step_data, self.ref.all_step_data[3])
        assert_equal(self.reader.time, self.ref.all_times[3])

    def test_next_past_last_step(self):
        self.reader.go_to_step(self.reader.n_steps-1)
        assert_raises(StopIteration, self.reader.next)

    def test_iter(self):
        for i, val in enumerate(self.reader):
            ref_data = {'time': self.ref.all_times[i], 
                        'data': self.ref.all_step_data[i]}
            assert_equal(val, ref_data)

    def test_read_ts_lowf(self):
        ## split up the following?
        self.reader.read_ts(self.ref.ts_lowf)
        assert_equal(self.reader.ts_data, self.ref.ts_lowf_data)
        assert_almost_equal(self.reader.ts_rep, self.ref.ts_lowf_rep)
        assert_almost_equal(self.ref.ts_lowf.aux.test, self.ref.ts_lowf_rep)
        assert_equal(self.reader.step, self.ref.ts_lowf_last_step)

    def test_read_ts_highf(self):
        ## split?
        self.reader.read_ts(self.ref.ts_highf)
        assert_equal(self.reader.ts_data, self.ref.ts_highf_data)
        assert_almost_equal(self.reader.ts_rep, self.ref.ts_highf_rep)
        assert_almost_equal(self.ref.ts_highf.aux.test, self.ref.ts_highf_rep)
        assert_equal(self.reader.step, self.ref.ts_highf_last_step)

    def test_ref_as_average(self):
        self.reader.represent_ts_as = 'average'
        self.reader.read_ts(self.ref.ts_lowf)
        assert_almost_equal(self.reader.ts_rep, self.ref.ts_lowf_rep_average)

    @raises(ValueError)
    def test_bad_represent_raises_ValueError(self):
        self.reader.represent_ts_as = 'invalid-option'

    @raises(ValueError)
    def test_time_col_out_of_range_raises_ValueError(self):
        self.reader.time_col = self.reader.n_cols 

    @raises(ValueError)
    def test_data_col_out_of_range_raises_ValueError(self):
        self.reader.data_cols = [self.reader.n_cols]

    ##TODO - file specific tests - opening/closing?

class XVGReference(BaseAuxReference):
    def __init__(self):
        super(XVGReference, self).__init__()
        self.testdata = AUX_XVG
        self.reader = mda.auxiliary.XVG.XVGReader

class TestXVGReader(BaseAuxReaderTest):
    def __init__(self):
        reference = XVGReference()
        super(TestXVGReader, self).__init__(reference)

