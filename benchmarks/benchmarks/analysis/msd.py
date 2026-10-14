"""
Benchmarks for the EinsteinMSD analysis module.
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import msd
from MDAnalysis.tests.datafiles import RANDOM_WALK_TOPO, RANDOM_WALK


class EinsteinMSDCustom:
    """
    Benchmark for Mean-Squared Displacement (EinsteinMSD).

    The trajectory is created by using random forces
    that are generated from a normal distribution
    and then cumulatively summed to simulate a real random walk.
    """

    unit = "ms"
    description = (
        "Performance of EinsteinMSD calculation using "
        "custom n_atoms and range of n_frames."
    )

    timeout = 300.0

    param_names = ["n_frames", "fft"]
    params = ([100, 1000, 5000], [True, False])

    def setup(self, n_frames, fft):
        """Setup method for MSD benchmark with custom number of frames."""
        n_atoms = 100

        rng = np.random.default_rng(42)

        # Random forces are sampled from a normal distribution
        steps = rng.standard_normal((n_frames, n_atoms, 3), dtype=np.float32)

        # Cumulative sum is taken to simulate positions over time
        res = np.cumsum(steps, axis=0)

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
    using the RANDOM_WALK example trajectory provided in the MDAnalysis testsuite.

    The files for the trajectory can be found in the library at :
        RANDOM_WALK      - mdanalysis/testsuite/xyz_random_walk.xtc
        RANDOM_WALK_TOPO - mdanalysis/testsuite/RANDOM_WALK_TOPO.pdb
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
