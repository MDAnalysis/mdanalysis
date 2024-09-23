import logging
import warnings

import numpy as np

from . import base

try:
    import gemmi

    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False

logger = logging.getLogger("MDAnalysis.coordinates.MMCIF")


def get_Coordinates(model: gemmi.Model) -> np.ndarray:
    return np.array(
        [[*at.pos.tolist()] for chain in model for res in chain for at in res]
    )


class MMCIFReader(base.SingleFrameReaderBase):
    """Reads from an MMCIF file"""

    format = "MMCIF"
    units = {"time": None, "length": "Angstrom"}

    def _read_first_frame(self):
        structure = gemmi.read_structure(self.filename)
        if len(structure) > 1:
            warnings.warn(
                f"File {self.filename} has {len(structure)} models, but only the first one will be read"
            )
        model = structure[0]
        coords = get_Coordinates(model)
        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(coords, **self._ts_kwargs)
        self.ts.frame = 0

    def Writer(self, filename, n_atoms=None, **kwargs):
        raise NotImplementedError

    def close(self):
        pass
