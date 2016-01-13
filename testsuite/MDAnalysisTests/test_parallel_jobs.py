import MDAnalysis as mda
import numpy as np
import MDAnalysis.analysis.parallel_jobs as pj
import MDAnalysis.analysis.electromagnetism as em

from numpy.testing import *
import time

from MDAnalysisTests.datafiles import (DCD, PSF)

class TestParallel(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.selection_string = 'all'

        # Single thread analysis
        start_time = time.time()
        single_analysis = em.TotalDipole(universe=self.universe, selection=self.selection_string)
        self.single_result = single_analysis.run()

    def test_parallel(self):
        start_time = time.time()
        jobs = [em.TotalDipole(selection=self.selection_string)]
        process = pj.ParallelProcessor(jobs,self.universe)
        assert_equal(self.single_result, process.parallel_run())

    def tearDown(self):
        del self.universe
