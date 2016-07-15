from numpy.testing import (assert_equal, assert_raises, assert_almost_equal,
                           raises)
import os

import MDAnalysis as mda

from MDAnalysisTests.datafiles import AUX_XVG, XVG_BAD_NCOL
from MDAnalysisTests.auxiliary.base import (BaseAuxReaderTest, BaseAuxReference)

class XVGReference(BaseAuxReference):
    def __init__(self):
        super(XVGReference, self).__init__()
        self.testdata = AUX_XVG
        self.reader = mda.auxiliary.XVG.XVGReader
        self.description['auxdata'] = os.path.abspath(self.testdata)
        self.description['format'] = self.reader.format

class TestXVGReader(BaseAuxReaderTest):
    def __init__(self):
        reference = XVGReference()
        super(TestXVGReader, self).__init__(reference)

    def test_no_time_selector(self):
        # XVGReader automatically sets time_select to 0 so we need to specifically
        # set it to none to test without
        self.reader = self.ref.reader(self.ref.testdata, dt=self.ref.dt, 
                                      initial_time=self.ref.initial_time,
                                      time_selector=None)
        for i, val in enumerate(self.reader):
            assert_equal(val.step_data, self.ref.all_data[i],
                         "step_data for step {0} does not match".format(i))

    @raises(ValueError)
    def test_wrong_n_col_raises_ValueError(self): 
        # encountering a different number of columns at a later step should 
        # raise ValueError
        self.reader = self.ref.reader(XVG_BAD_NCOL)
        next(self.reader)


class XVGFileReference(XVGReference):
    def __init__(self):
        super(XVGFileReference, self).__init__()
        self.reader = mda.auxiliary.XVG.XVGFileReader
        self.format = "XVG-F"
        self.description['format'] = self.format

class TestXVGFileReader(TestXVGReader):
    def __init__(self):
        reference = XVGFileReference()
        super(TestXVGReader, self).__init__(reference)

    def test_get_auxreader_for(self):
        # Default reader of .xvg files is intead XVGReader, not XVGFileReader
        # so test specifying format 
        reader = mda.auxiliary.core.get_auxreader_for(self.ref.testdata,
                                                      format=self.ref.format)
        assert_equal(reader, self.ref.reader)

    def test_reopen(self):
        self.reader._reopen()
        # should start us back at before step 0, so next takes us to step 0
        self.reader.next()
        assert_equal(self.reader.auxstep.step_data, self.ref.all_step_data[0])
