from __future__ import division, absolute_import, print_function

try:
    from MDAnalysis.coordinates.DCD import DCDReader
    from MDAnalysis.coordinates.XTC import XTCReader
    from MDAnalysis.coordinates.TRR import TRRReader
    from MDAnalysisTests.datafiles import XTC, TRR, DCD
except ImportError:
    pass

class TrajReaderBench(object):
    """Benchmarks for trajectory file format reading."""
    params = (['XTC', 'TRR', 'DCD'])
    param_names = ['traj_format']

    def setup(self, traj_format):
        self.traj_dict = {'XTC': [XTC, XTCReader],
                          'TRR': [TRR, TRRReader],
                          'DCD': [DCD, DCDReader]}
        self.traj_file, self.traj_reader = self.traj_dict[traj_format]

    def time_reads(self, traj_format):
        """Simple benchmark for reading traj file formats
        from our standard test files.
        """
        self.traj_reader(self.traj_file)

class TrajReaderStriding(object):
    """Benchmarks for trajectory file format striding."""
    params = (['XTC', 'TRR', 'DCD'])
    param_names = ['traj_format']

    def setup(self, traj_format):
        self.traj_dict = {'XTC': [XTC, XTCReader],
                          'TRR': [TRR, TRRReader],
                          'DCD': [DCD, DCDReader]}
        self.traj_file, self.traj_reader = self.traj_dict[traj_format]
        self.reader_object = self.traj_reader(self.traj_file)

    def time_strides(self, traj_format):
        """Benchmark striding over full trajectory
        test files for each format.
        """
        for ts in self.reader_object:
            pass
