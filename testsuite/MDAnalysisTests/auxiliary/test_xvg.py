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
from numpy.testing import (assert_equal, assert_raises, assert_almost_equal,
                           assert_array_equal, raises)
import numpy as np

import os

import MDAnalysis as mda

from MDAnalysisTests.datafiles import AUX_XVG, XVG_BAD_NCOL, XVG_BZ2
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
        self.select_time_ref = np.arange(self.n_steps)
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

def test_xvg_bz2():
    reader = mda.auxiliary.XVG.XVGReader(XVG_BZ2)
    assert_array_equal(reader.read_all_times(), np.array([0., 50., 100.]))

def test_xvg_bz2():
    reader = mda.auxiliary.XVG.XVGFileReader(XVG_BZ2)
    assert_array_equal(reader.read_all_times(), np.array([0., 50., 100.]))
