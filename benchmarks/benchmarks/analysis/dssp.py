import MDAnalysis as mda

try:
    from MDAnalysisTests.datafiles import TPR, XTC, PDB_janin
except:
    pass

try:
    from MDAnalysis.analysis.dssp import DSSP
    from MDAnalysis.analysis.dssp.pydssp_numpy import assign
except:
    pass


class DsspBench(object):
    """Benchmarks for MDAnalysis.analysis.dssp.DSSP over a trajectory."""

    params = ([1, 5, 10], [True, False])
    param_names = ["n_frames", "guess_hydrogens"]

    def setup(self, n_frames, guess_hydrogens):
        self.u = mda.Universe(TPR, XTC)
        self.dssp = DSSP(self.u, guess_hydrogens=guess_hydrogens)
        self.n_frames = n_frames

    def time_dssp_run(self, n_frames, guess_hydrogens):
        """Benchmark DSSP.run() over n_frames of a trajectory."""
        self.dssp.run(stop=self.n_frames)


class DsspAssignBench(object):
    """Benchmarks for the core DSSP assign() function on a single frame.

    Uses PDB_janin as a large test system for number of residues
    """

    params = [50, 100, 214, 500]
    param_names = ["n_residues"]

    def setup(self, n_residues):
        self.u = mda.Universe(PDB_janin)
        dssp = DSSP(self.u)
        coords = dssp._get_coords()
        self.coords = coords[:n_residues]

    def time_assign(self, n_residues):
        """Benchmark single-frame assign() call."""
        assign(self.coords)
