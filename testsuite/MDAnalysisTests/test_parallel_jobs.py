import MDAnalysis as mda
import numpy as np
import MDAnalysis.analysis.parallel_jobs as pj
import MDAnalysis.analysis.electromagnetism as electromagnetism

from numpy.testing import *
import time

from MDAnalysisTests.datafiles import (DCD,PSF)

class TestParallel(TestCase):
    def setUp(self):
        self.system	= mda.Universe(PSF,DCD)
        self.selection  = self.system.select_atoms('all')
        self.traj       = self.system.trajectory

        # Single thread analysis
        print "Single thread analysis:"
        start_time = time.time()
        single_analysis = electromagnetism.TotalDipole(self.traj,selection=self.selection)
        self.single_result = single_analysis.run()
        print("--- %s seconds ---" % (time.time() - start_time))

    def test_parallel(self):
        print "Parallel analysis:"
        start_time = time.time()
        jobs = [ electromagnetism.TotalDipole(selection = self.selection) ]
        process = pj.ParallelProcessor(jobs,self.traj,progressbar=True)
        assert_equal(self.single_result, process.parallel_run())
        print("--- %s seconds ---" % (time.time() - start_time))

    def tearDown(self):
        del self.system
