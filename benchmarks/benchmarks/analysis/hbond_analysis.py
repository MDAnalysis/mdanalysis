"""
Benchmark for HydrogenBondAnalysis
"""

try:
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
    from MDAnalysisTests.datafiles import waterPSF , waterDCD
except ImportError:
    pass

import MDAnalysis

class HydrogenBondAnalysisBenchmark:
    """
    It tests performance of hbond.run() across different
    number of frames using waterPSF/waterDCD test files.
    """

    unit = "ms"
    timeout = 60.0
    params = [2, 5, 10]
    param_names = ["n_frames"]


    def setup(self, n_frames):
        u = MDAnalysis.Universe(waterPSF, waterDCD)
        self.hbonds = HBA(universe= u)
    
    def time_run(self, n_frames):
        self.hbonds.run(stop = n_frames)
        