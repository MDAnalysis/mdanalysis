from MDAnalysis.analysis.base import AnalysisBase
import numpy as np

class RMSFResidue(AnalysisBase):
    """
    Compute RMSF (Root Mean Square Fluctuation) per residue.
    """

    def __init__(self, universe, select="all", **kwargs):
        self.universe = universe
        self.select = select

        # Select atoms
        self.ag = universe.select_atoms(select)

        # Build residue -> atomgroup mapping manually (safe & correct)
        self.residue_groups = []
        for res in universe.residues:
            self.residue_groups.append(self.ag[self.ag.resids == res.resid])

        super().__init__(self.universe.trajectory, **kwargs)

    def _prepare(self):
        # Number of residues
        self.n_residues = len(self.residue_groups)

        # Residue metadata
        self.resids = np.array([res.resid for res in self.universe.residues])
        self.resnames = np.array([res.resname for res in self.universe.residues])

        # Running statistics
        self._sum = np.zeros((self.n_residues, 3))
        self._sum_sq = np.zeros((self.n_residues, 3))
        self._counts = 0

    def _single_frame(self):
        # For each residue compute mean atomic position
        for i, group in enumerate(self.residue_groups):
            if len(group) > 0:
                mean_pos = group.positions.mean(axis=0)
                self._sum[i] += mean_pos
                self._sum_sq[i] += mean_pos ** 2

        self._counts += 1

    def _conclude(self):
        # Compute RMSF
        mean = self._sum / self._counts
        mean_sq = self._sum_sq / self._counts
        self.rmsf = np.sqrt((mean_sq - mean**2).sum(axis=1))

        # Store results
        self.results.resids = self.resids
        self.results.resnames = self.resnames
        self.results.residue_rmsf = self.rmsf
