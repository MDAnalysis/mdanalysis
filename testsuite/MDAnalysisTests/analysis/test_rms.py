# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function

from six.moves import range

import MDAnalysis
import MDAnalysis.analysis.rms

from numpy.testing import TestCase, assert_almost_equal
import numpy as np

import os
import tempfile

from MDAnalysisTests.datafiles import GRO, XTC, rmsfArray

class TestRMSF(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(GRO, XTC)
        fd, self.outfile = tempfile.mkstemp(suffix=".xtc")  # output is always same as input (=XTC)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

        del self.universe

    def test_rmsf(self):
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(quiet=True)
        test_rmsfs = np.load(rmsfArray)

        assert_almost_equal(rmsfs.rmsf, test_rmsfs, 5,
                            err_msg="error: rmsf profile should match test " +
                            "values")

    def test_rmsf_single_frame(self):
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(start=5, stop=6, quiet=True)

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be zero")

    def test_rmsf_identical_frames(self):
        # write a dummy trajectory of all the same frame
        with MDAnalysis.Writer(self.outfile, self.universe.atoms.n_atoms) as W:
            for i in range(self.universe.trajectory.n_frames):
                W.write(self.universe)

        self.universe = MDAnalysis.Universe(GRO, self.outfile)
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(quiet=True)

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be 0")


