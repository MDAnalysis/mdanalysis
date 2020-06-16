
try:
    from MDAnalysis.coordinates.DCD import DCDReader
    from MDAnalysisTests.datafiles import DCD
except ImportError:
    pass

try:
    from MDAnalysis.coordinates.XTC import XTCReader
    from MDAnalysisTests.datafiles import XTC
except ImportError:
    pass

try:
    from MDAnalysis.coordinates.TRR import TRRReader
    from MDAnalysisTests.datafiles import TRR
except ImportError:
    pass

try:
    from MDAnalysis.coordinates.TRJ import NCDFReader
    from MDAnalysisTests.datafiles import NCDF
except ImportError:
    pass

traj_dict = {'XTC': [XTC, XTCReader],
             'TRR': [TRR, TRRReader],
             'DCD': [DCD, DCDReader],
             'NCDF': [NCDF, NCDFReader]}

class TrajReaderCreation(object):
    """Benchmarks for trajectory file format reading."""
    params = (['XTC', 'TRR', 'DCD', 'NCDF'])
    param_names = ['traj_format']

    def setup(self, traj_format):
        self.traj_dict = traj_dict
        self.traj_file, self.traj_reader = self.traj_dict[traj_format]

    def time_reads(self, traj_format):
        """Simple benchmark for reading traj file formats
        from our standard test files.
        """
        self.traj_reader(self.traj_file)


class TrajReaderIteration(object):
    """Benchmarks for trajectory file format striding."""
    params = (['XTC', 'TRR', 'DCD', 'NCDF'])
    param_names = ['traj_format']

    def setup(self, traj_format):
        self.traj_dict = traj_dict
        self.traj_file, self.traj_reader = self.traj_dict[traj_format]
        self.reader_object = self.traj_reader(self.traj_file)

    def time_strides(self, traj_format):
        """Benchmark striding over full trajectory
        test files for each format.
        """
        for ts in self.reader_object:
            pass
