"""
Benchmarks for the EinsteinMSD analysis module.
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import msd
from MDAnalysis.tests.datafiles import RANDOM_WALK_TOPO, RANDOM_WALK


class EinsteinMSDCustom:
    """
    Benchmark for Mean-Squared Displacement (EinsteinMSD)
    using custom value of n_atoms and range of n_frames.
    """

    unit = "ms"
    description = (
        "Performance of EinsteinMSD calculation using "
        "custom n_atoms and range of n_frames."
    )

    timeout = 300.0

    param_names = ["n_frames", "fft"]
    params = ([100, 1000, 10000], [True, False])

    def setup(self, n_frames, fft):
        """Setup method for MSD benchmark with custom number of frames."""
        np.random.seed(42)
        n_atoms = 100
        res = np.random.random((n_frames, n_atoms, 3)).astype(np.float32)

        self.u = mda.Universe.empty(n_atoms, n_frames=n_frames, trajectory=True)
        self.u.trajectory.set_array(res)

    def time_msd_run_custom(self, n_frames, fft):
        """
        Benchmark calculation of Mean-Squared Displacement
        using a custom trajectory.
        """
        msd_var = msd.EinsteinMSD(self.u, select="all", fft=fft)
        msd_var.run()


class EinsteinMSDExample:
    """
    Benchmark for Mean-Squared Displacement (EinsteinMSD)
    using the example trajectory RANDOM_WALK.
    """

    unit = "ms"
    description = (
        "Performance of EinsteinMSD calculation using RANDOM_WALK_TOPO and RANDOM_WALK"
    )

    timeout = 300.0

    param_names = ["fft"]
    params = [True, False]

    def setup(self, fft):
        """Setup method for MSD benchmark with RANDOM_WALK example."""
        self.u = mda.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)

    def time_msd_run_simple(self, fft):
        """
        Benchmark calculation of Mean-Squared Displacement
        using the example trajectory RANDOM_WALK.
        """
        msd_var = msd.EinsteinMSD(self.u, fft=fft)
        msd_var.run()
