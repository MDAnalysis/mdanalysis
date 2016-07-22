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

        # add the auxdata and format for .xvg to the reference description
        self.description['auxdata'] = os.path.abspath(self.testdata)
        self.description['format'] = self.reader.format

        # for testing the selection of data/time
        self.time_selector = 0 # take time as first value in auxilairy
        self.select_time_ref = range(self.n_steps)
        self.data_selector = [1,2] # select the second/third columns from auxiliary
        self.select_data_ref = [self.format_data([2*i, 2**i]) for i in range(self.n_steps)]


class TestXVGReader(BaseAuxReaderTest):
    def __init__(self):
        reference = XVGReference()
        super(TestXVGReader, self).__init__(reference)

    @raises(ValueError)
    def test_changing_n_col_raises_ValueError(self): 
        # if number of columns in .xvg file is not consistent, a ValueError
        # should be raised
        self.reader = self.ref.reader(XVG_BAD_NCOL)
        next(self.reader)

    @raises(ValueError)
    def test_time_selector_out_of_range_raises_ValueError(self):
        # if time_selector is not a valid index of _data, a ValueError 
        # should be raised
        self.reader.time_selector = len(self.reader.auxstep._data) 

    @raises(ValueError)
    def test_data_selector_out_of_range_raises_ValueError(self):
        # if data_selector is not a valid index of _data, a ValueError 
        # should be raised
        self.reader.data_selector = [len(self.reader.auxstep._data)]


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
        assert_equal(self.reader.step, 0)
