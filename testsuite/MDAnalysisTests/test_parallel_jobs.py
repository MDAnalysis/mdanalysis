import MDAnalysis as mda
import numpy as np
import MDAnalysis.analysis.parallel_jobs as pj
import MDAnalysis.analysis.electromagnetism as em

from numpy.testing import *

from MDAnalysisTests.datafiles import (DCD, PSF)

class TestParallel(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.selection_string = 'all'

        # Single thread analysis
        single_analysis = em.TotalDipole(universe=self.universe, selection=self.selection_string)
        self.single_result = single_analysis.run()

    def test_parallel(self):
        jobs = [em.TotalDipole(selection=self.selection_string)]
        process = pj.ParallelProcessor(jobs,self.universe)
        assert_equal(self.single_result, process.parallel_run())

    def test_parallel_base(self):
        single_analysis = em.TotalDipole(universe=self.universe, selection=self.selection_string)
        assert_equal(self.single_result, single_analysis.run(parallel=True))

    def tearDown(self):
        del self.universe
