import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD


class TimeTrajectoryIteration:
    def setup(self):
        self.u = mda.Universe(PSF, DCD)

    def time_iterate_trajectory(self):
        for ts in self.u.trajectory:
            pass
