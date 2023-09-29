from .base import AnalysisBase
import numpy as np
from .. import Universe
import pydssp


class DSSP(AnalysisBase):
    def __init__(self, u: Universe):

        super().__init__(u.trajectory)

        selections = {
            t: u.select_atoms(f"protein and {t}")
            for t in (
                "backbone and name N",
                "backbone and name CA",
                "backbone and name C",
                "name O or name O1",
                # hydrogens don't work yet -- they give less atoms than there actually are
                # "name H1 or name H",
            )
        }
        self.selections = selections
        self.n_protein_atoms = sum(map(lambda s: s.n_atoms, selections.values()))
        self._trajectory = u.trajectory

    def _prepare(self):
        self.results.dssp = [None for _ in range(self.n_frames)]

    def _single_frame(self):
        coords = np.array(
            [group.positions for group in self.selections.values()]
        ).swapaxes(0, 1)
        dssp = pydssp.assign(coords)
        self.results.dssp[self._frame_index] = dssp

    def _conclude(self):
        pass
